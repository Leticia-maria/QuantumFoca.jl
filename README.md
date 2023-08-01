# THIS REPO WAS TRANSFERRED TO: [Ohata.jl](https://github.com/HartreeFoca/Ohata.jl)

# QuantumFoca.jl

[![Issues](https://img.shields.io/github/issues-raw/Leticia-maria/Foca.jl?style=for-the-badge)](https://github.com/Leticia-maria/QuantumFoca.jl/)
[![Build Status](https://img.shields.io/github/workflow/status/Leticia-maria/Foca.jl/CI?style=for-the-badge)](https://github.com/Leticia-maria/QuantumFoca.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Commit Activity](https://img.shields.io/github/commit-activity/w/Leticia-maria/Foca.jl/main?style=for-the-badge)](https://github.com/Leticia-maria/QuantumFoca.jl/)
[![DOI](https://zenodo.org/badge/419452183.svg)](https://zenodo.org/badge/latestdoi/419452183)

![HartreeFoca jl](https://user-images.githubusercontent.com/60739184/170071106-68ba0e42-08a5-4923-b69a-d5db945bdf7b.svg)

## Overview

## Installation

To install the package, you will call the Julia Package Manager on your REPL:

```julia
]add QuantumFoca
```

Done! Now it is *ready to use*
## Package Features

- *Calculates overlap integrals*
- *Calculate kinetic integrals*
- *Calculate electron-nuclear attraction integrals*
- *Calculate electron-electron repulsion integrals*

## Quick Example

Consider that you want to calculate the electronic energy of a methane molecule. First, you will need the cartesian coordinates of the molecule of interest. This information is stored in our ```methane.xyz``` file, formatted as follows:

```julia
5

C  0.00001021434087  0.00001532972083 -0.00001493500137
H -0.19951695340554  0.87894179053067 -0.62713882127936
H  0.76712229809243  0.24863902907755  0.74526241504934
H  0.35580334399536 -0.82601803138729 -0.62993342769733
H -0.92343260142312 -0.30159515034176  0.51179839372872
```

```julia
methane = molecule("methane.xyz")
```

For any molecular calculations, you will need a basis set.

```julia
sto3g = buildbasis(methane)
```

With this information, we can build the molecular integrals.

```julia
S = overlap(sto3g, methane)
T = kinetic(sto3g, methane)
V = attraction(sto3g, methane)
G = repulsion(sto3g, methane)
```
