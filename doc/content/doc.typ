#import "../lib.typ": *

This doc is written for description of the TBSK code. The original motivation was to simulate the ARPES spectrum, noticed the chinook package, so I made my own SKTB model of the code base, more convenient to continue on this.
= What are Tight Binding models?
== Tight Binding Model
Tight binding models describe the electronic structure of solids and molecules by expanding electronic wave functions in a set of *localized basis functions*, typically *atomic orbitals*. The expansion coefficients (tight-binding parameters) are obtained by fitting to experimental data or first-principles calculations.

These models allow efficient calculation of electronic properties like *band structure, density of states, and conductivity*—especially effective when only nearest-neighbor interactions matter.

Construction methods include the Slater-Koster approach, which considers only nearest-neighbor interactions, and the extended tight binding approach, which includes longer-range interactions.

== Slater-Koster Tight Binding
The Slater-Koster method models electronic structures using a simplified Hamiltonian with *nearest-neighbor interactions*. It approximates the system with localized atomic orbitals connected by tight-binding parameters, fitted from experiments or first-principles.

This method efficiently computes properties such as band structure and conductivity, and is valued for its simplicity and accuracy in systems dominated by short-range interactions.
= Basic framework
qwerasdf
231
ewq

ewqe

#figure(image("fig/1.png",width:40%))



// #definition(
//   "Stokes' theorem",
//   footer: "Information extracted from a well-known public encyclopedia"
// )[
//   Let $Sigma$ be a smooth oriented surface in $RR^3$ with boundary $diff Sigma
//   equiv Gamma$. If a vector field $iboxed(bold(F)(x,y,z))=(F_x (x,y,z), F_y (x,y,z),
//   F_z (x,y,z))$ is defined and has continuous first order partial derivatives
//   in a region containing $Sigma$, then

//   $ integral.double_Sigma (bold(nabla) times bold(F)) dot bold(Sigma) =
//   dboxed(integral.cont_(diff Sigma) bold(F) dot dif bold(Gamma)) $ 
// ] <stokes>

