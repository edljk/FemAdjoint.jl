# FemAdjoint

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://edljk.github.io/FemAdjoint.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://edljk.github.io/FemAdjoint.jl/dev)


Designed to illustrate the use of `SparseDiffTools` approach in the context of simple finite element cost functions.

+ P1 triangular mesh, gradient computation with respect to the mesh points of $|| u || ^ 2$ where $K(p)u = M(p) \mathbb{1}$ by the adjoint method.


## Install FemAdjoint

```julia
>] add https://github.com/edljk/FemAdjoint.jl
```

## Benchmark

```julia
using FemAdjoint
>] test FemAdjoint
```