using JLD2

const CACHE_FILE = "vortex_cache.jld2"

function cache_key(p::SolverParams)
    return (p.gap_inf, p.temperature, p.V, p.coherence_length,
            p.use_hartree, p.Ngrid, p.L, p.n_vortex, p.m, p.mu)
end

function load_cache(p::SolverParams)
    isfile(CACHE_FILE) || return nothing
    jldopen(CACHE_FILE, "r") do f
        key = string(cache_key(p))
        haskey(f, key) ? (f[key]["xs"], f[key]["ys"], f[key]["Delta"], f[key]["U"]) : nothing
    end
end

function save_cache(p::SolverParams, xs, ys, Delta, U)
    jldopen(CACHE_FILE, "a+") do f
        key = string(cache_key(p))
        f[key] = Dict("xs" => xs, "ys" => ys, "Delta" => Delta, "U" => U)
    end
end
