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
    mu::Float64 = 1.0         # chemical potential
    max_iter::Int = 200
end
