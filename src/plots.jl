using CairoMakie
using Colors

function plot_results(xs::Vector{Float64}, ys::Vector{Float64},
                      Delta::Matrix{ComplexF64}, params::SolverParams;
                      vals::Union{Vector{Float64}, Nothing} = nothing,
                      us::Union{Matrix{ComplexF64}, Nothing} = nothing,
                      vs::Union{Matrix{ComplexF64}, Nothing} = nothing,
                      outdir=".")
    gap_inf = params.gap_inf
    Ngrid   = params.Ngrid

    Delta_init = reshape(
        gl_profile_2d(xs, ys, params.gap_inf, params.coherence_length, params.n_vortex),
        Ngrid, Ngrid)

    ξ    = params.coherence_length
    xs_ξ = xs ./ ξ
    ys_ξ = ys ./ ξ

    # --- 1. Phase heatmap: hue = arg(Δ), saturation = |Δ|/Δ₀, value = 1 ---
    hue     = @. (angle(Delta) + π) / (2π)
    sat     = clamp.(abs.(Delta) ./ gap_inf, 0.0, 1.0)
    rgb_img = [convert(RGB, HSV(hue[i,j] * 360, sat[i,j], 1.0))
               for i in 1:Ngrid, j in 1:Ngrid]
    fig1 = Figure()
    ax1  = Axis(fig1[1,1]; xlabel="x / ξ", ylabel="y / ξ",
                title="Phase [hue]  ·  |Δ|/Δ₀ [saturation]", aspect=DataAspect())
    image!(ax1, (xs_ξ[1], xs_ξ[end]), (ys_ξ[1], ys_ξ[end]), rgb_img)
    # Hidden heatmap anchors the Colorbar to the hsv phase mapping
    hm_phase = heatmap!(ax1, [0,1], [0,1], [-π π]; colormap=:hsv, visible=false)
    Colorbar(fig1[1,2], hm_phase; label="arg(Δ)",
             ticks=([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]))
    save(joinpath(outdir, "delta_phase.png"), fig1)

    # --- 2. X-cut |Δ(x, 0)|/Δ₀ comparing initial GL guess vs converged ---
    iy   = div(Ngrid, 2) + 1
    fig2 = Figure()
    ax2  = Axis(fig2[1,1]; xlabel="x / ξ", ylabel="|Δ| / Δ₀", title="|Δ(x, y=0)| / Δ₀")
    lines!(ax2, xs_ξ, abs.(Delta_init[:, iy]) ./ gap_inf; label="initial (GL)", linestyle=:dash)
    lines!(ax2, xs_ξ, abs.(Delta[:, iy])      ./ gap_inf; label="converged")
    hlines!(ax2, [1.0]; linestyle=:dot, color=:gray, label="Δ₀")
    axislegend(ax2, position = :rb)
    save(joinpath(outdir, "delta_xcut.png"), fig2)

    # --- 3. Quasiparticle spectrum (if eigenstates provided) ---
    if !isnothing(vals)
        plot_energy_levels(vals, params; outdir)
    end

    # --- 4. Lowest bound-state density (if eigenstates provided) ---
    if !isnothing(vals) && !isnothing(us) && !isnothing(vs)
        plot_lowest_bound_state(vals, us, vs, params, xs, ys; outdir)
    end

    @info "Saved plots to $outdir"
end

"""
Scatter plot of quasiparticle energies within ±1.5Δ₀ of the Fermi level,
marking the Fermi level (E=0) and the bulk gap edges (±Δ₀).
Saves spectrum.png.
"""
function plot_energy_levels(vals::Vector{Float64}, p::SolverParams; outdir=".")
    window    = 1.5 * p.gap_inf
    sub_vals  = filter(e -> abs(e) < window, vals)
    indices   = 1:length(sub_vals)

    sub_vals_scaled = sub_vals ./ p.gap_inf

    fig = Figure()
    ax  = Axis(fig[1,1];
               xlabel="State index (sorted)", ylabel="Energy (ε / Δ₀)",
               title="Quasiparticle Spectrum (CdGM States)")
    scatter!(ax, indices, sub_vals_scaled; color=:dodgerblue, markersize=4)
    hlines!(ax, [0.0];        color=:black, linestyle=:dot,  linewidth=2)
    hlines!(ax, [1.0, -1.0]; color=:red,   linestyle=:dash, linewidth=1.5)
    text!(ax, "+Δ₀ (continuum)"; position=(length(sub_vals)*0.05,  1.12), color=:red, fontsize=10)
    text!(ax, "-Δ₀ (continuum)"; position=(length(sub_vals)*0.05, -1.12), color=:red, fontsize=10)
    save(joinpath(outdir, "spectrum.png"), fig)
end

"""
Heatmap of |u|² + |v|² for the lowest positive-energy eigenstate.
Saves bound_state.png.
"""
function plot_lowest_bound_state(vals::Vector{Float64},
                                  us::Matrix{ComplexF64}, vs::Matrix{ComplexF64},
                                  p::SolverParams,
                                  xs::Vector{Float64}, ys::Vector{Float64}; outdir=".", up_to_n_level=4)
    pos_indices = findall(>(0.0), vals)
    n_plot = min(up_to_n_level, length(pos_indices))

    if n_plot == 0
        @warn "No positive energy states found!"
        return
    end
    @info "Plotting the first $n_plot bound states..."

    cols = ceil(Int, sqrt(n_plot))
    rows = ceil(Int, n_plot / cols)
    fig  = Figure(size=(300 * cols, 320 * rows))
    Label(fig[0, 1:cols], "Core Bound States (|u|² + |v|²)"; fontsize=14)

    local hm_last
    for i in 1:n_plot
        idx        = pos_indices[i]
        E_val      = vals[idx]
        density_2d = reshape(abs2.(us[:, idx]) .+ abs2.(vs[:, idx]), p.Ngrid, p.Ngrid)

        row = div(i-1, cols) + 1
        col = mod(i-1, cols) + 1
        ax  = Axis(fig[row, col];
                   title="n=$i  (E ≈ $(round(E_val/p.gap_inf, digits=3)) Δ₀)",
                   titlesize=10, aspect=DataAspect())
        hidedecorations!(ax)
        hm_last = heatmap!(ax, xs ./ p.coherence_length, ys ./ p.coherence_length,
                           density_2d; colormap=:viridis)
    end
    Colorbar(fig[1:rows, cols+1], hm_last)

    save(joinpath(outdir, "bound_state.png"), fig)
end
