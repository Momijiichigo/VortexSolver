using Pkg
Pkg.activate(@__DIR__)

# include("src/VortexSolver.jl")
using VortexSolver


p = make_solver_params(
    gap_inf          = 0.5,
    temperature      = 0.1,
    use_hartree      = false,
    mixing           = 0.3,
    tol              = 1e-2,
    Ngrid            = 35,
    L                = 14.0,
    n_vortex         = 1,
)

# cached = load_cache(p)
# if cached !== nothing
#     xs, ys, Delta, U = cached
#     @info "Loaded from cache"
# else
#     xs, ys, Delta, U, _ = run_scf(p)
#     save_cache(p, xs, ys, Delta, U)
# end

xs, ys, Delta, U, _ = run_scf(p)

plot_results(xs, ys, Delta, p; outdir=@__DIR__)
