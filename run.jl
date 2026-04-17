using Pkg
Pkg.activate(@__DIR__)

# include("src/VortexSolver.jl")
using VortexSolver


p = make_solver_params(
    gap_inf          = 0.5,
    temperature      = 0.1,
    use_hartree      = false,
    mixing           = 0.3,
    tol              = 4e-3,
    Ngrid            = 41,
    L                = 14.0,
    n_vortex         = 1,
    m                = 1.0,
)

cached = load_cache(p)
if cached !== nothing
    xs, ys, Delta, U, vals, us, vs = cached
    @info "Loaded from cache"
else
    xs, ys, Delta, U, _, vals, us, vs = run_scf(p)
    save_cache(p, xs, ys, Delta, U, vals, us, vs)
end

plot_results(xs, ys, Delta, p; vals, us, vs, outdir=joinpath(@__DIR__, "plots"))
