using Pkg
Pkg.activate(@__DIR__)

include("src/VortexSolver.jl")
using .VortexSolver

p = SolverParams(
    gap_inf          = 1.0,
    temperature      = 0.1,
    V                = 2.0,
    coherence_length = 1.0,
    use_hartree      = false,
    mixing           = 0.3,
    tol              = 1e-6,
    Ngrid            = 40,
    L                = 10.0,
    n_vortex         = 1,
)

cached = load_cache(p)
if cached !== nothing
    xs, ys, Delta, U = cached
    @info "Loaded from cache"
else
    xs, ys, Delta, U, _ = run_scf(p)
    save_cache(p, xs, ys, Delta, U)
end

plot_results(xs, ys, Delta, p; outdir=@__DIR__)
