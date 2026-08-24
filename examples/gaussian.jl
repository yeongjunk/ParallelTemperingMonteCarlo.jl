using Random
using Statistics
using ParallelTemperingMonteCarlo

import ParallelTemperingMonteCarlo:
    AbstractReplicas,
    step!,
    step_slot!,
    getenergy,
    getstate,
    getbeta,
    getwalkerid,
    swapwalkers!

mutable struct GaussianReplicas{R<:AbstractRNG} <: AbstractReplicas
    states::Vector{Vector{Float64}}  # states[walker]
    betas::Vector{Float64}           # betas[slot]
    walkerids::Vector{Int}           # walkerids[slot]
    proposal_width::Float64
    rng::R
end

function GaussianReplicas(
    betas::AbstractVector;
    proposal_width=0.5,
    rng=Xoshiro(1234),
)
    K = length(betas)
    states = [[randn(rng)] for _ in 1:K]

    return GaussianReplicas(
        states,
        Float64.(betas),
        collect(1:K),
        proposal_width,
        rng,
    )
end

Base.length(reps::GaussianReplicas) = length(reps.betas)

getwalkerid(reps::GaussianReplicas, slot::Int) =
    reps.walkerids[slot]

getbeta(reps::GaussianReplicas, slot::Int) =
    reps.betas[slot]

getstate(reps::GaussianReplicas, slot::Int) =
    reps.states[getwalkerid(reps, slot)]

function getenergy(reps::GaussianReplicas, slot::Int)
    x = getstate(reps, slot)[1]
    return x^2 / 2
end

function swapwalkers!(
    reps::GaussianReplicas,
    slot_i::Int,
    slot_j::Int,
)
    reps.walkerids[slot_i], reps.walkerids[slot_j] =
        reps.walkerids[slot_j], reps.walkerids[slot_i]

    return nothing
end

function step_slot!(reps::GaussianReplicas, slot::Int)
    walker = getwalkerid(reps, slot)
    x = reps.states[walker][1]
    x_proposed = x + reps.proposal_width * randn(reps.rng)

    energy = x^2 / 2
    energy_proposed = x_proposed^2 / 2
    beta = getbeta(reps, slot)

    log_acceptance = -beta * (energy_proposed - energy)

    if log(rand(reps.rng)) < min(0.0, log_acceptance)
        reps.states[walker][1] = x_proposed
        return true
    end

    return false
end

function step!(reps::GaussianReplicas)
    return [step_slot!(reps, slot) for slot in 1:length(reps)]
end

# Replica setup
betas = [0.5, 1.0, 2.0, 4.0]
reps = GaussianReplicas(betas)

edge_groups = [
    [(1, 2), (3, 4)],
    [(2, 3)],
]

exchange_params = ExchangeParams(edge_groups, 10)

sampling_params = SamplingParams(
    100_000,        # number of sweeps
    10,             # sample every
    10_000,         # partition every
    [length(betas)], # sample the coldest slot
)

result = sample_replicas!(
    reps,
    sampling_params,
    exchange_params;
    rng=Xoshiro(5678),
    sample_eltype=Float64,
)

samples = vec(result.samples[1, 1, :])
beta = betas[end]

println("Sample mean:       ", mean(samples))
println("Sample variance:   ", var(samples))
println("Expected variance: ", 1 / beta)
