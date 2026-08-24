# ParallelTemperingMonteCarlo.jl

[![CI](https://github.com/yeongjunk/ParallelTemperingMonteCarlo.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yeongjunk/ParallelTemperingMonteCarlo.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A lightweight, algorithm-independent implementation of parallel tempering
(replica-exchange Monte Carlo) for Julia.

This package provides the replica-exchange workflow while leaving the local
Monte Carlo update to the user. It can therefore be combined with algorithms
such as random-walk Metropolis, Langevin Monte Carlo (LMC), or Hamiltonian
Monte Carlo (HMC).

## Installation

Until the package is registered in Julia's General registry, install it
directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/yeongjunk/ParallelTemperingMonteCarlo.jl")
```

## Interface

Users must define the following methods for their own subtype of
`AbstractReplicas`:

```julia
step!(reps::AbstractReplicas) =
    error("step!(reps) is not implemented.")

step_slot!(reps::AbstractReplicas, slot::Int) =
    error("step_slot!(reps, slot) is not implemented.")

getenergy(reps::AbstractReplicas, slot::Int) =
    error("getenergy(reps, slot) is not implemented.")

getstate(reps::AbstractReplicas, slot::Int) =
    error("getstate(reps, slot) is not implemented.")

getbeta(reps::AbstractReplicas, slot::Int) =
    error("getbeta(reps, slot) is not implemented.")

getwalkerid(reps::AbstractReplicas, slot::Int) =
    error("getwalkerid(reps, slot) is not implemented.")

swapwalkers!(reps::AbstractReplicas, slot_i::Int, slot_j::Int) =
    error("swapwalkers!(reps, slot_i, slot_j) is not implemented.")

Base.length(reps::AbstractReplicas) =
    error("length(reps) is not implemented.")
```

Temperature parameters belong to fixed slots, while states belong to walkers.
Replica exchange is implemented by changing the mapping from slots to walkers.

## Example

A self-contained [Gaussian example](examples/gaussian.jl) implements a minimal
`AbstractReplicas` subtype using random-walk Metropolis updates. It runs
parallel tempering and compares the sampled variance with the exact result
`1 / beta`.

Run it from the repository root with:

```bash
julia --project=. examples/gaussian.jl
```

For an optimized Langevin Monte Carlo implementation, see
[LMC.jl](https://github.com/yeongjunk/LMC.jl).

## Testing

Run the test suite with:

```julia
using Pkg
Pkg.test()
```

## License

ParallelTemperingMonteCarlo.jl is released under the [MIT License](LICENSE).
