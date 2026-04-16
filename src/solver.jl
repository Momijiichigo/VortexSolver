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
"""

using LinearAlgebra, SparseArrays

fermi(ε, T) = T == 0.0 ? (ε < 0 ? 1.0 : 0.0) : 1.0 / (exp(ε / T) + 1.0)

"""
Build the 2D finite-difference Laplacian on an Ngrid×Ngrid grid with
spacing h, using Dirichlet (zero) boundary conditions.
Returns a sparse (Ngrid²)×(Ngrid²) matrix.
"""
function laplacian_2d(Ngrid::Int, h::Float64)
    N2 = Ngrid^2
    Is = Int[]; Js = Int[]; Vs = Float64[]

    function idx(i, j)
        return (i-1)*Ngrid + j
    end
    function add!(i, j, val)
        push!(Is, i); push!(Js, j); push!(Vs, val)
    end

    for i in 1:Ngrid, j in 1:Ngrid
        k = idx(i, j)
        boundary = (i == 1 || i == Ngrid || j == 1 || j == Ngrid)
        if boundary
            add!(k, k, 1.0)
        else
            add!(k, k,            -4.0 / h^2)
            add!(k, idx(i-1, j),   1.0 / h^2)
            add!(k, idx(i+1, j),   1.0 / h^2)
            add!(k, idx(i, j-1),   1.0 / h^2)
            add!(k, idx(i, j+1),   1.0 / h^2)
        end
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
Assemble and diagonalize the 2D BdG Hamiltonian for the given complex Δ(x,y)
and real U(x,y) (both as flat vectors of length Ngrid²).
Returns eigenvalues ε and (u, v) eigenvectors (each column is one eigenstate).
"""
function solve_bdg(p::SolverParams, lap::SparseMatrixCSC,
                   Delta::Vector{ComplexF64}, U::Vector{Float64})
    N2 = p.Ngrid^2

    He = Matrix(-lap ./ (2p.m) - p.mu * I + Diagonal(U))

    D  = Diagonal(Delta)
    Dc = Diagonal(conj.(Delta))

    H = ComplexF64[He  D;
                   Dc  -He]

    vals, vecs = eigen(Hermitian(H))
    us = vecs[1:N2, :]
    vs = vecs[N2+1:end, :]
    return vals, us, vs
end

"""
Update Δ(x,y) and U(x,y) from the BdG eigensolutions.
"""
function update_fields(p::SolverParams,
                       vals::Vector{Float64},
                       us::Matrix{ComplexF64},
                       vs::Matrix{ComplexF64})
    N2 = p.Ngrid^2
    Delta_new = zeros(ComplexF64, N2)
    U_new     = zeros(Float64, N2)
    for n in eachindex(vals)
        fn = fermi(vals[n], p.temperature)
        Delta_new .+= p.V .* conj.(vs[:, n]) .* us[:, n] .* (1 - 2fn)
        if p.use_hartree
            U_new .-= p.V .* (abs2.(us[:, n]) .* fn .+ abs2.(vs[:, n]) .* (1 - fn))
        end
    end
    return Delta_new, U_new
end

"""
Run the SCF loop on a 2D Cartesian grid.
Returns (xs, ys, Delta_converged, U_converged, converged::Bool)
where Delta is a complex Ngrid×Ngrid matrix.
"""
function run_scf(p::SolverParams)
    @info "starting scf"
    h  = 2p.L / (p.Ngrid - 1)
    xs = range(-p.L, p.L; length=p.Ngrid) |> collect
    ys = xs

    lap   = laplacian_2d(p.Ngrid, h)
    Delta = gl_profile_2d(xs, ys, p.gap_inf, p.coherence_length, p.n_vortex)
    U     = zeros(p.Ngrid^2)

    converged = false
    for iter in 1:p.max_iter
        vals, us, vs = solve_bdg(p, lap, Delta, U)
        Delta_new, U_new = update_fields(p, vals, us, vs)

        err = maximum(abs.(Delta_new .- Delta))
        Delta = p.mixing .* Delta_new .+ (1 - p.mixing) .* Delta
        U     = p.mixing .* U_new     .+ (1 - p.mixing) .* U

        if err < p.tol
            converged = true
            @info "SCF converged in $iter iterations (err=$err)"
            break
        end
        @info "iter=$iter  err=$(round(err; sigdigits=4))"
    end
    converged || @warn "SCF did not converge within $(p.max_iter) iterations"

    Delta_2d = reshape(Delta, p.Ngrid, p.Ngrid)
    U_2d     = reshape(U,     p.Ngrid, p.Ngrid)
    return xs, ys, Delta_2d, U_2d, converged
end
