"""
2D BdG self-consistent field solver for an s-wave vortex core.

Works directly on a Cartesian (x,y) grid without the gauge-rotation or
angular-separation trick, so no rotational symmetry is assumed.

The BdG eigenvalue problem (ℏ=1):

  [H_e,       Δ(x,y) ] [u]   [u]
  [Δ*(x,y),  -H_e    ] [v] = ε[v]

  H_e = -1/(2m) ∇² - μ + U(x,y)

Self-consistency (from mean-field decoupling):
  Δ(x,y) = V Σ_n v_n*(x,y) u_n(x,y) (1 - 2fₙ)
  U(x,y) = -V Σ_n (|u_n|² fₙ + |v_n|² (1−fₙ))    [if use_hartree]

GPU architecture: Delta and U live on the device for the entire SCF loop;
only the static He_base, boundary mask, and final results are ever on the CPU.
"""

using LinearAlgebra, SparseArrays, CUDA

fermi(ε, T) = T == 0.0 ? (ε < 0 ? 1.0 : 0.0) : 1.0 / (exp(ε / T) + 1.0)

"""
Build the 2D finite-difference Laplacian on an Ngrid×Ngrid grid with
spacing h, using Dirichlet (zero) boundary conditions.
Returns a sparse (Ngrid²)×(Ngrid²) matrix.
"""
function laplacian_2d(Ngrid::Int, h::Float64)
    N2 = Ngrid^2
    Is = Int[]; Js = Int[]; Vs = Float64[]

    idx(i, j) = (i-1)*Ngrid + j

    for i in 1:Ngrid, j in 1:Ngrid
        k = idx(i, j)
        push!(Is, k); push!(Js, k); push!(Vs, -4.0 / h^2)
        if i > 1;     push!(Is, k); push!(Js, idx(i-1, j)); push!(Vs, 1.0 / h^2); end
        if i < Ngrid; push!(Is, k); push!(Js, idx(i+1, j)); push!(Vs, 1.0 / h^2); end
        if j > 1;     push!(Is, k); push!(Js, idx(i, j-1)); push!(Vs, 1.0 / h^2); end
        if j < Ngrid; push!(Is, k); push!(Js, idx(i, j+1)); push!(Vs, 1.0 / h^2); end
    end
    return sparse(Is, Js, Vs, N2, N2)
end

"""
GL initial profile: Δ(x,y) = Δ₀ tanh(r/ξ) exp(−i n_vortex θ).
Returns a complex vector of length Ngrid².
"""
function gl_profile_2d(xs::Vector{Float64}, ys::Vector{Float64},
                        gap_inf::Float64, xi::Float64, n_vortex::Int)
    Ngrid = length(xs)
    Delta = zeros(ComplexF64, Ngrid^2)
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        r = sqrt(x^2 + y^2)
        theta = atan(y, x)
        Delta[(i-1)*Ngrid + j] = gap_inf * tanh(r / xi) * exp(-im * n_vortex * theta)
    end
    return Delta
end

"""
Assemble and diagonalize the 2D BdG Hamiltonian entirely on the GPU.
He_base_gpu is the static -∇²/(2m) - μ part (uploaded once before the SCF loop).
Delta_gpu and U_gpu are device arrays updated each iteration.
Returns vals_gpu, us_gpu, vs_gpu as CuArrays — no CPU transfer here.
"""
function solve_bdg(p::SolverParams, He_base_gpu::CuMatrix{ComplexF64},
                   Delta_gpu::CuVector{ComplexF64}, U_gpu::CuVector{Float64})
    N2 = p.Ngrid^2

    # Add the density-dependent Hartree shift on the GPU
    He_gpu = He_base_gpu .+ Diagonal(U_gpu)

    # Assemble 2N²×2N² BdG Hamiltonian on the GPU using block-view assignment.
    # @views avoids temporary copies; Diagonal(CuVector) broadcasts correctly into
    # CuArray views without scalar indexing.
    H_gpu = CUDA.zeros(ComplexF64, 2N2, 2N2)
    @views H_gpu[1:N2,     1:N2]      .= He_gpu
    @views H_gpu[N2+1:2N2, N2+1:2N2] .= -He_gpu
    @views H_gpu[1:N2,     N2+1:2N2] .= Diagonal(Delta_gpu)
    @views H_gpu[N2+1:2N2, 1:N2]     .= Diagonal(conj.(Delta_gpu))

    vals_gpu, vecs_gpu = eigen(Hermitian(H_gpu))

    us_gpu = vecs_gpu[1:N2, :]
    vs_gpu = vecs_gpu[N2+1:2N2, :]
    return vals_gpu, us_gpu, vs_gpu
end

"""
Vectorized gap and Hartree update on the GPU — no scalar indexing or per-state loops.

Gap update (vectorized over all 2N² states at once):
  Δ_new = (V/h²) · (conj(V_mat) .* U_mat) · w_gap
  where w_gap[n] = (1 − 2fₙ) · [εₙ > 0]

Hartree update:
  U_new = −(V/h²) · (|U_mat|² · f + |V_mat|² · (1−f))
"""
function update_fields(p::SolverParams, h::Float64,
                       vals_gpu::CuVector{Float64},
                       us_gpu::CuMatrix{ComplexF64},
                       vs_gpu::CuMatrix{ComplexF64})
    fn = fermi.(vals_gpu, p.temperature)

    # Only positive-energy states enter the gap equation
    pos_mask    = vals_gpu .> 0
    weights_gap = (1 .- 2 .* fn) .* pos_mask      # (2N²,) real

    # Delta_new[k] = (V/h²) Σ_n conj(v[k,n]) u[k,n] w_gap[n]
    #              = (V/h²) · row-wise dot of (conj(V) .* U) with w_gap
    Delta_new = (p.V / h^2) .* ((conj.(vs_gpu) .* us_gpu) * weights_gap)

    if p.use_hartree
        U_new = -(p.V / h^2) .* (abs2.(us_gpu) * fn .+ abs2.(vs_gpu) * (1 .- fn))
    else
        U_new = CUDA.zeros(Float64, p.Ngrid^2)
    end

    return Delta_new, U_new
end

"""
Run the SCF loop on a 2D Cartesian grid, keeping Delta and U on the GPU
throughout. Only final results are transferred to the CPU.
Returns (xs, ys, Delta_converged, U_converged, converged::Bool)
where Delta and U are Ngrid×Ngrid CPU matrices.
"""
function run_scf(p::SolverParams)
    @info "starting scf"
    h  = 2p.L / (p.Ngrid - 1)
    xs = range(-p.L, p.L; length=p.Ngrid) |> collect
    ys = xs
    N2 = p.Ngrid^2

    # --- One-time CPU work: build static He_base, upload to GPU ---
    lap     = laplacian_2d(p.Ngrid, h)
    He_base = ComplexF64.(Matrix(-lap ./ (2p.m)))
    He_base[diagind(He_base)] .-= p.mu
    He_base_gpu = CuArray(He_base)

    # --- Initial fields: GL profile on CPU, then push to GPU ---
    Delta_gpu = CuArray(gl_profile_2d(xs, ys, p.gap_inf, p.coherence_length, p.n_vortex))
    U_gpu     = CUDA.zeros(Float64, N2)

    # --- Precompute boundary mask and pinned GL values (uploaded once) ---
    # boundary_mask = zeros(Bool, N2)
    # Delta_bc      = zeros(ComplexF64, N2)
    # for i in 1:p.Ngrid, j in 1:p.Ngrid
    #     if i == 1 || i == p.Ngrid || j == 1 || j == p.Ngrid
    #         k           = (i-1)*p.Ngrid + j
    #         theta       = atan(ys[j], xs[i])
    #         boundary_mask[k] = true
    #         Delta_bc[k]      = p.gap_inf * exp(-im * p.n_vortex * theta)
    #     end
    # end
    # boundary_mask_gpu = CuArray(boundary_mask)
    # Delta_bc_gpu      = CuArray(Delta_bc)

    vals_gpu = us_gpu = vs_gpu = nothing
    converged = false
    for iter in 1:p.max_iter
        vals_gpu, us_gpu, vs_gpu = solve_bdg(p, He_base_gpu, Delta_gpu, U_gpu)
        Delta_new_gpu, U_new_gpu = update_fields(p, h, vals_gpu, us_gpu, vs_gpu)

        # Pin boundary to the bulk GL profile — fully vectorized, no scalar indexing
        # @. Delta_new_gpu = ifelse(boundary_mask_gpu, Delta_bc_gpu, Delta_new_gpu)

        # Convergence check: GPU reduction → CPU scalar
        err = maximum(abs.(Delta_new_gpu .- Delta_gpu)) / p.gap_inf

        # Linear mixing on GPU
        @. Delta_gpu = p.mixing * Delta_new_gpu + (1 - p.mixing) * Delta_gpu
        @. U_gpu     = p.mixing * U_new_gpu     + (1 - p.mixing) * U_gpu

        if err < p.tol
            converged = true
            @info "SCF converged in $iter iterations (err=$err)"
            break
        end
        @info "iter=$iter  err=$(round(err; sigdigits=4))"
    end
    converged || @warn "SCF did not converge within $(p.max_iter) iterations"

    # Transfer all final results to CPU only here
    Delta_2d = reshape(Array(Delta_gpu), p.Ngrid, p.Ngrid)
    U_2d     = reshape(Array(U_gpu),     p.Ngrid, p.Ngrid)
    vals     = Array(vals_gpu)
    us       = Array(us_gpu)
    vs       = Array(vs_gpu)

    return xs, ys, Delta_2d, U_2d, converged, vals, us, vs
end
