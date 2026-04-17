using CairoMakie
using Colors

function plot_results(xs::Vector{Float64}, ys::Vector{Float64},
                      Delta::Matrix{ComplexF64}, params::SolverParams; outdir=".")
    gap_inf = params.gap_inf
    Ngrid   = params.Ngrid

    Delta_init = reshape(
        gl_profile_2d(xs, ys, params.gap_inf, params.coherence_length, params.n_vortex),
        Ngrid, Ngrid)

    # --- 1. Phase heatmap: hue = arg(Δ), saturation = |Δ|/Δ₀, value = 1 ---
    hue     = @. (angle(Delta) + π) / (2π)
    sat     = clamp.(abs.(Delta) ./ gap_inf, 0.0, 1.0)
    rgb_img = [convert(RGB, HSV(hue[i,j] * 360, sat[i,j], 1.0))
               for i in 1:Ngrid, j in 1:Ngrid]
    fig1 = Figure()
    ax1  = Axis(fig1[1,1]; xlabel="x", ylabel="y",
                title="Phase [hue]  ·  |Δ|/Δ₀ [saturation]", aspect=DataAspect())
    image!(ax1, (xs[1], xs[end]), (ys[1], ys[end]), rgb_img)
    save(joinpath(outdir, "delta_phase.png"), fig1)

    # --- 2. X-cut |Δ(x, 0)| comparing initial GL guess vs converged ---
    iy   = div(Ngrid, 2) + 1
    fig2 = Figure()
    ax2  = Axis(fig2[1,1]; xlabel="x", ylabel="|Δ|", title="|Δ(x, y=0)|")
    lines!(ax2, xs, abs.(Delta_init[:, iy]); label="initial (GL)", linestyle=:dash)
    lines!(ax2, xs, abs.(Delta[:, iy]);      label="converged")
    hlines!(ax2, [gap_inf]; linestyle=:dot, color=:gray, label="Δ₀")
    axislegend(ax2, position = :rb)
    save(joinpath(outdir, "delta_xcut.png"), fig2)

    @info "Saved delta_phase.png and delta_xcut.png to $outdir"
end
