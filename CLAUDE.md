# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A self-consistent Bogoliubov–de Gennes (BdG) solver for s-wave superconductor vortex cores, written in Julia with Makie for visualization. The physics roadmap has two stages:

1. **Default Hamiltonian** — cylindrically symmetric BdG with a single vortex (`n=1` winding); active focus.
2. **Bernal Bilayer Graphene** — replace the free-electron `H_0` with the BBG lattice Hamiltonian.

## Physics Summary (`draft.typ`)

The mean-field Hamiltonian is solved on a radial grid via the BdG eigenvalue problem. Key physical quantities:

- **Gap** `Δ(r)` — self-consistently determined; near a vortex takes the form `Δ(r) = Δ̂(r) e^{-iθ}`.
- **Hartree term** `U(r)` — optional density-dependent shift; toggled by `use_hartree`.
- **Quasi-particle spinor** — the solver works directly with `(u(x,y), v(x,y))` on a 2D Cartesian grid; no gauge rotation or angular separation is applied, so no rotational symmetry is assumed.

Self-consistency loop:
1. Initialize `Δ(r)` from Ginzburg-Landau profile.
2. Solve the radial BdG eigenvalue problem to get `{uₙ(r), vₙ(r)}`.
3. Update `Δ(r)` and (optionally) `U(r)` from the new eigenvectors and Fermi weights.
4. Mix old/new solutions, check convergence against `tol`, repeat.

## Planned Solver Parameters

| Parameter | Meaning |
|---|---|
| `gap_inf` | Bulk gap `Δ₀` |
| `temperature` | Sets Fermi weights `fₙ` |
| `V` | Attractive pairing potential |
| `coherence_length` | Initial GL profile scale |
| `use_hartree` | Include `U(r)` in `H_e` |
| `mixing` | Linear mixing ratio for SCF update |
| `tol` | Convergence threshold |
| `L` / `Ngrid` | Half-width of square grid and points per side |

## Expected Output / Plots

- `|Δ(r)|` vs `r`
- Phase colormap of `Δ` (`colormap = :hsv`)
- Cache file: converged `Δ` (and optionally `U`) keyed by parameters — solver skips to plotting on cache hit.

## Commands

```bash
# Run the solver (from VortexSolver/)
julia --project=. run.jl

# Compile physics notes
typst compile draft.typ
```

## Code Architecture (`VortexSolver/`)

All source files live in `src/`:

- **`params.jl`** — `SolverParams` struct (all tunable inputs via `Base.@kwdef`).
- **`solver.jl`** — core physics: `laplacian_2d` (sparse 2D FD Laplacian), `gl_profile_2d` (complex GL initial profile), `solve_bdg` (dense complex Hermitian BdG matrix diagonalization), `update_fields`, `run_scf`. The SCF loop lives entirely here.
- **`cache.jl`** — JLD2-backed cache keyed by a tuple of all `SolverParams` fields; `load_cache`/`save_cache`.
- **`plots.jl`** — `plot_results`: radial profile and 2D phase colormap using CairoMakie.
- **`VortexSolver.jl`** — module root; includes all files above.

Entry point is `run.jl` (project root of `VortexSolver/`): constructs `SolverParams`, checks cache, runs SCF, saves cache, plots.

## Toolchain

- **Julia** — solver and main codebase
- **Makie** — visualization (currently `CairoMakie`; swap to `GLMakie` for interactive use)
- **JLD2** — parameter-keyed cache (`vortex_cache.jld2` written next to `run.jl`)
- **Typst** — physics notes/draft (`draft.typ`); compile with `typst compile draft.typ`

## Development cycles

As you will notice if there will be any modifications (additional parameters, additional options, etc.), carefully design the structures and update `CLAUDE.md` appropriately.
