"""
Parameters for the 2D s-wave vortex core BdG solver.

The solver works on a Cartesian (x,y) grid with the full complex BdG Hamiltonian,
without the gauge-rotation or angular-separation trick. This supports geometries
where continuous rotational symmetry does not hold.
"""
Base.@kwdef struct SolverParams
    gap_inf::Float64          # bulk gap Δ₀
    temperature::Float64      # temperature (same units as gap_inf; ℏ=1)
    V::Float64                # attractive pairing potential
    coherence_length::Float64 # GL coherence length ξ (sets initial profile)
    use_hartree::Bool = false
    mixing::Float64 = 0.3     # linear mixing: Δ_new = mixing*Δ_BdG + (1-mixing)*Δ_old
    tol::Float64 = 1e-6       # convergence threshold (max |ΔΔ|)
    Ngrid::Int = 40           # grid points per side (total grid: Ngrid × Ngrid)
    L::Float64 = 10.0         # half-width of square grid (grid spans [-L,L]×[-L,L])
    n_vortex::Int = 1         # vortex winding number for initial GL phase profile
    m::Float64 = 1.0          # effective mass (ℏ=1 units)
    mu::Float64 = 4.0         # chemical potential
    max_iter::Int = 200
end


function make_solver_params(;
    gap_inf::Float64 = 0.5,       # You now provide the target gap directly
    temperature::Float64,
    use_hartree::Bool = false,
    mixing::Float64 = 0.3,
    tol::Float64 = 1e-6,
    Ngrid::Int = 41,              # Strongly recommend an odd number to hit (0,0) exactly!
    L::Float64 = 10.0,
    n_vortex::Int = 1,
    m::Float64 = 1.0,
    mu::Float64 = 4.0,
    max_iter::Int = 200,
)
    # 1. Compute BCS coherence length analytically
    v_F = sqrt(2 * mu / m)
    coherence_length = v_F / (pi * gap_inf)
    @info "Coherence length calculated: ξ₀ = $(round(coherence_length, digits=4))"

    # 2. Create a temporary params object (V is a dummy value for now)
    p_tmp = SolverParams(
        gap_inf=gap_inf, temperature=temperature, V=1.0, coherence_length=coherence_length,
        use_hartree=use_hartree, mixing=mixing, tol=tol, Ngrid=Ngrid, L=L,
        n_vortex=n_vortex, m=m, mu=mu, max_iter=max_iter,
    )

    # 3. Compute exact V required to sustain your chosen gap_inf
    V_calc = calculate_V(p_tmp)

    # 4. Return final fully-consistent params
    return SolverParams(
        gap_inf=gap_inf, temperature=temperature, V=V_calc, coherence_length=coherence_length,
        use_hartree=use_hartree, mixing=mixing, tol=tol, Ngrid=Ngrid, L=L,
        n_vortex=n_vortex, m=m, mu=mu, max_iter=max_iter,
    )
end

function calculate_V(p::SolverParams)
    @info "Computing exact V to support Δ₀ = $(p.gap_inf) for an N=$(p.Ngrid) grid"
    h = 2 * p.L / (p.Ngrid - 1)
    N = p.Ngrid
    sum_val = 0.0
    
    # Sum exactly over the N^2 discrete momentum states of the grid
    for mx in 0:(N-1), my in 0:(N-1)
        kx = (2.0 * pi * mx) / (N * h)
        ky = (2.0 * pi * my) / (N * h)
        
        # Exact kinetic energy dispersion for a 2D finite-difference grid
        kinetic = (2.0 - cos(kx * h) - cos(ky * h)) / (p.m * h^2)
        xi_k = kinetic - p.mu
        
        Ek = sqrt(xi_k^2 + p.gap_inf^2)
        f_k = fermi(Ek, p.temperature)
        
        sum_val += (1.0 - 2.0 * f_k) / (2.0 * Ek)
    end
    
    # Real-space normalization: sum_val / (N^2 * h^2) perfectly matches the
    # (V / h^2) prefactor and the 1/N^2 eigenvector normalization in the BdG loop.
    V_required = 1.0 / (sum_val / (N^2 * h^2))
    
    @info "Pairing potential calculated: V = $(round(V_required, digits=4))"
    return V_required
end
