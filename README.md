## Interface

Users must define the following methods for their own subtype of `AbstractReplicas`:

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

This package provides only the abstract interface and replica-exchange workflow. Users are responsible for implementing the underlying Monte Carlo algorithm, such as Langevin Monte Carlo (LMC) or Hamiltonian Monte Carlo (HMC).

For an implementation example, see the parallel-tempering section of my LMC package.

