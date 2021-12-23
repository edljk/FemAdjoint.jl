# FemAdjoint

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://edljk.github.io/FemAdjoint.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://edljk.github.io/FemAdjoint.jl/dev)
[![Build Status](https://travis-ci.com/edljk/FemAdjoint.jl.svg?branch=main)](https://travis-ci.com/edljk/FemAdjoint.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/edljk/FemAdjoint.jl?svg=true)](https://ci.appveyor.com/project/edljk/FemAdjoint-jl)
[![Build Status](https://api.cirrus-ci.com/github/edljk/FemAdjoint.jl.svg)](https://cirrus-ci.com/github/edljk/FemAdjoint.jl)
[![Coverage](https://codecov.io/gh/edljk/FemAdjoint.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/edljk/FemAdjoint.jl)
[![Coverage](https://coveralls.io/repos/github/edljk/FemAdjoint.jl/badge.svg?branch=main)](https://coveralls.io/github/edljk/FemAdjoint.jl?branch=main)

Designed to illustrate the use of `SparseDiffTools` approach in the context of simple finite element cost functions.

+ P1 triangular mesh, gradient computation with respect to the mesh points of $|| u || ^ 2$ where $K(p)u = M(p) \mathbb{1}$ by the adjoint method.


## Install FemAdjoint

```julia
import Pkg
Pkg.clone("https://github.com/edljk/FemAdjoint.jl")
```

## Benchmark

```julia
using FemAdjoint
>] test FemAdjoint
```