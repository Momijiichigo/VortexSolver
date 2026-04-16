#import "@preview/physica:0.9.7": *
#let nop = $hat(n)$
#let Hop = $hat(H)$
#let up = $arrow.t$
#let down = $arrow.b$
#let exte = $epsilon.alt$
#let gap = $Delta$
#let nambuSpinorK = $vec(a_(k up), a_(-k down)^dagger)$
#let nambuSpinorKdag = $vecrow(a^dagger_(k up), a_(-k down))$
#let quasiOp = $vec(gamma_(k up),gamma^dagger_(-k down))$
#let gndphi = $tilde(phi.alt)$
#let psid = $psi^dagger$
#let ga = $gamma$
#let gad = $gamma^dagger$
#let dr = $d r$
#let avg(arg) = $expectationvalue(arg)$

== S-wave SC
$
H_"eff" = integral dr { sum_sigma psid_sigma overbrace((cal(H)_0 +U(r)), cal(H)_e) psi_sigma + gap(r)psid_up (r)psid_down (r) + gap^*(r)psi_down (r)psi_up (r) + abs(gap)^2/V}
$

- Hartree term: $U(r)=-V avg(psid_up (r)psi_up (r))=-V avg(psid_down (r)psi_down (r))$
- Gap term: $gap(r)=-V avg(psi_down (r)psi_up (r))=+V avg(psi_up (r)psi_down (r))$

Substitution with $psi_sigma=sum_n gamma_(n, sigma) u_n - sigma gamma^dagger_(n, -sigma)v_n^*$ will result in:

$
U(r) = -V sum_n (|u_n (r)|^2f_n + |v_n (r)|^2 (1-f_n))\

gap(r) =+V sum_n (v_n^* (r) u_n (r)(1-2f_n))
$

where $f_n$ is Fermi-Dirac distribution with quasi-particle excitation energy $exte_n$ plugged in.

== Vortex Core
#let gapR = $hat(gap)$
#let uvvec = $vec(u (r), v (r))$

Around the vortex, the $gap$ becomes $gap(r) = gapR(r) e^(-i theta)$.
    
$exte_n uvvec = mat(cal(H)_e, gapR(r) e^(-i theta); gapR(r) e^(i theta), -cal(H)^*_e) uvvec$ 



== Computational solver for s-wave SC vortex core

Since the above vortex core equation is self-consistent equation:
- the equation requires $gap(r)$ which depends on $uvvec$
- $uvvec$ is obtained from $f(r)$
- etc.

We would like to implement a self-consistent field solver through the given diff-eqn, and iterate until it converges.
The brief steps will be
- set the condition $gap_infinity$, temperature, attractive potential, coherence length, use_hartree (bool), mixing, tol, grid side length, etc.
- Start by the initial $gap(r)$ given by Ginzburg-Landau
- iteration of plug the result values into the diff-eq until it converges (error is within `tol`)
- Plot:
  - |Δ(r)| vs r
  - phase colormap of Δ (colormap = :hsv)
- Save the converged $gap$ (and optionally $U$) with the used parameters into a cache file. 
  - the program will first read the cache file and if the parameter is the same, it will skip the solver step and directly plots using the saved data.

=== Roadmap
1. the default Hamiltonian (given above); *We will focus on this for now.*
2. Include Bernal Bilayer Graphene lattice structure for $cal(H)_0$

== Tools

- Julia
- Makie
