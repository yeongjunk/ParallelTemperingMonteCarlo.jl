module RR 

using Random, LinearAlgebra

export AbstractReplica, step!, energy, cross_energy, beta, replica_exchange!, try_exchange!,
    Walker, ExchangeParams, ExchangeStatus, WalkerStatus,
    make_walkers, equilibrate!, monitor_equilibration!,
    exchange_sweep!, attempt_exchange_group!, compute_exchange_rates,
    check_edge_groups, 
    getstate, SamplingParams, sample_replicas!
## AbstractReplica interface

abstract type AbstractReplica end

step!(r::AbstractReplica; rng=Random.GLOBAL_RNG) = error("step! is not implemented.")
energy(r::AbstractReplica) = error("energy is not implemented.")
cross_energy(ri::AbstractReplica, rj::AbstractReplica) = error("cross_energy is not implemented.")
beta(r::AbstractReplica) = error("beta is not implemented.")
replica_exchange!(ri::AbstractReplica, rj::AbstractReplica) = error("replica_exchange! is not implemented.")
replica_exchange!(ri::AbstractReplica, rj::AbstractReplica, Eij, Eji) = replica_exchange!(ri, rj) 
getstate(r::AbstractReplica) = error("getstate is not implemented.")

function try_exchange!(ri::AbstractReplica, rj::AbstractReplica; rng=Random.GLOBAL_RNG)
    βi = beta(ri)
    βj = beta(rj)

    Eii = energy(ri)
    Ejj = energy(rj)
    Eij = cross_energy(ri, rj)
    Eji = cross_energy(rj, ri)

    logα = -βi*Eij - βj*Eji + βi*Eii + βj*Ejj

    if isfinite(logα) && log(rand(rng)) < min(0.0, logα)
        replica_exchange!(ri, rj, Eij, Eji)
        return true
    end

    return false
end


## Exchange parameters and status updates

struct ExchangeParams
    edge_groups::Vector{Vector{Tuple{Int,Int}}}
    n_between_exchange::Int
    n_sweeps::Int
end

struct ExchangeStatus
    trials::Vector{Vector{Int}}
    hits::Vector{Vector{Int}}
    rates::Vector{Vector{Float64}}
end

function ExchangeStatus(params::ExchangeParams)
    n_groups = length(params.edge_groups)

    trials = Vector{Vector{Int}}(undef, n_groups)
    hits   = Vector{Vector{Int}}(undef, n_groups)
    rates  = Vector{Vector{Float64}}(undef, n_groups)

    for g in 1:n_groups
        n_edges = length(params.edge_groups[g])

        trials[g] = zeros(Int, n_edges)
        hits[g]   = zeros(Int, n_edges)
        rates[g]  = zeros(Float64, n_edges)
    end

    return ExchangeStatus(trials, hits, rates)
end

function compute_exchange_rates(status::ExchangeStatus)
    n_groups = length(status.trials)

    exchange_rates = Vector{Vector{Float64}}(undef, n_groups)

    for g in 1:n_groups
        n_edges = length(status.trials[g])

        exchange_rates[g] = Vector{Float64}(undef, n_edges)

        for e in 1:n_edges
            n_try = status.trials[g][e]
            n_hit = status.hits[g][e]

            if n_try == 0
                exchange_rates[g][e] = 0.0
            else
                exchange_rates[g][e] = n_hit / n_try
            end
        end
    end

    return exchange_rates
end

function update_exchange_rates!(status::ExchangeStatus)
    n_groups = length(status.trials)

    for g in 1:n_groups
        n_edges = length(status.trials[g])

        length(status.hits[g])  == n_edges || error("hits[$g] has inconsistent length.")
        length(status.rates[g]) == n_edges || error("rates[$g] has inconsistent length.")

        for e in 1:n_edges
            n_try = status.trials[g][e]
            n_hit = status.hits[g][e]

            if n_try == 0
                status.rates[g][e] = 0.0
            else
                status.rates[g][e] = n_hit / n_try
            end
        end
    end

    return status
end

function check_edge_groups(edge_groups::Vector{Vector{Tuple{Int,Int}}}, len_replica::Int)
    K = len_replica

    length(edge_groups) > 0 || error("edge_groups must be nonempty.")

    for edges in edge_groups
        used = falses(K)

        for (i, j) in edges
            1 <= i <= K || error("edge index out of range.")
            1 <= j <= K || error("edge index out of range.")
            i != j || error("self edge is not allowed.")
            !used[i] || error("edge group has duplicated vertex.")
            !used[j] || error("edge group has duplicated vertex.")

            used[i] = true
            used[j] = true
        end
    end

    return nothing
end


## Walker interface

mutable struct Walker{TR <: AbstractReplica}
    id::Int
    replica::TR
end

mutable struct WalkerStatus
    slot_counts::Matrix{Int}
    last_endpoint::Vector{Int}
    endpoint_crossings::Vector{Int}
    hot_visits::Vector{Int}
    cold_visits::Vector{Int}
    k_hot::Int
    k_cold::Int
end

function WalkerStatus(K::Int; k_hot::Int=1, k_cold::Int=K)
    1 <= k_hot <= K || error("Require 1 <= k_hot <= K.")
    1 <= k_cold <= K || error("Require 1 <= k_cold <= K.")
    k_hot != k_cold || error("k_hot and k_cold must be different.")

    return WalkerStatus(zeros(Int, K, K), zeros(Int, K), zeros(Int, K), zeros(Int, K), zeros(Int, K), k_hot, k_cold)
end

function update_walker_status!(mon::WalkerStatus, walkers::Vector{<:Walker})
    K = length(walkers)

    size(mon.slot_counts) == (K, K) || error("slot_counts must have size K × K.")

    @inbounds for k in 1:K
        w = walkers[k].id

        mon.slot_counts[w, k] += 1

        is_hot  = (k == mon.k_hot)
        is_cold = (k == mon.k_cold)

        is_hot && (mon.hot_visits[w] += 1)
        is_hot && (mon.last_endpoint[w] == 2) && (mon.endpoint_crossings[w] += 1)
        is_hot && (mon.last_endpoint[w] = 1)

        is_cold && (mon.cold_visits[w] += 1)
        is_cold && (mon.last_endpoint[w] == 1) && (mon.endpoint_crossings[w] += 1)
        is_cold && (mon.last_endpoint[w] = 2)
    end

    return nothing
end

## Exchange attempt primitives

function attempt_exchange_group!(walkers::Vector{<:Walker}, edges::Vector{Tuple{Int,Int}}; rng=Random.GLOBAL_RNG)
    n_edges = length(edges)
    accepted = Vector{Bool}(undef, n_edges)

    @inbounds for e in 1:n_edges
        i, j = edges[e]
        accepted[e] = try_exchange!(walkers[i].replica, walkers[j].replica; rng=rng)

        if accepted[e]
            walkers[i].id, walkers[j].id = walkers[j].id, walkers[i].id
        end
    end

    return accepted
end

function attempt_exchange_group!(walkers::Vector{<:Walker}, stats::ExchangeStatus, group_index::Int, edges::Vector{Tuple{Int,Int}}; rng=Random.GLOBAL_RNG)
    n_edges = length(edges)
    g = group_index

    length(stats.trials[g]) == n_edges || error("trials[$g] must have length equal to number of edges.")
    length(stats.hits[g])   == n_edges || error("hits[$g] must have length equal to number of edges.")

    accepted = attempt_exchange_group!(walkers, edges; rng=rng)

    @inbounds for e in 1:n_edges
        stats.trials[g][e] += 1
        accepted[e] && (stats.hits[g][e] += 1)
    end

    return accepted
end

function exchange_sweep!(walkers::Vector{<:Walker}, params::ExchangeParams; rng=Random.GLOBAL_RNG)
    accepted_groups = Vector{Vector{Bool}}(undef, length(params.edge_groups))

    for g in eachindex(params.edge_groups)
        accepted_groups[g] = attempt_exchange_group!(walkers, params.edge_groups[g]; rng=rng)
    end

    return accepted_groups
end

function exchange_sweep!(walkers::Vector{<:Walker}, params::ExchangeParams, stats::ExchangeStatus; rng=Random.GLOBAL_RNG)
    accepted_groups = Vector{Vector{Bool}}(undef, length(params.edge_groups))

    for g in eachindex(params.edge_groups)
        accepted_groups[g] = attempt_exchange_group!(walkers, stats, g, params.edge_groups[g]; rng=rng)
    end

    update_exchange_rates!(stats)

    return accepted_groups
end


## Utilities

function make_replica_rngs(K::Int; rng=Random.GLOBAL_RNG)
    seeds = rand(rng, UInt, K)
    return [Xoshiro(seed) for seed in seeds]
end

function make_walkers(replicas::Vector{TR}) where {TR <: AbstractReplica}
    walkers = Vector{Walker{TR}}(undef, length(replicas))

    for k in eachindex(replicas)
        walkers[k] = Walker{TR}(k, replicas[k])
    end

    return walkers
end


## Equilibration

function equilibrate!(walkers::Vector{<:Walker}, params::ExchangeParams; rng=Random.GLOBAL_RNG, exchange_stats=ExchangeStatus(params), walker_status=nothing)
    K = length(walkers)

    check_edge_groups(params.edge_groups, K)

    rngs = make_replica_rngs(K; rng=rng)

    for sweep in 1:params.n_sweeps
        @Threads.threads for k in 1:K
            for _ in 1:params.n_between_exchange
                step!(walkers[k].replica; rng=rngs[k])
            end
        end

        exchange_sweep!(walkers, params, exchange_stats; rng=rng)

        walker_status !== nothing && update_walker_status!(walker_status, walkers)
    end

    update_exchange_rates!(exchange_stats)

    return exchange_stats, walker_status
end

function equilibrate!(replicas::Vector{TR}, params::ExchangeParams; rng=Random.GLOBAL_RNG, walker_status=nothing) where {TR <: AbstractReplica}
    walkers = make_walkers(replicas)
    exchange_stats = ExchangeStatus(params)

    return equilibrate!(walkers, params; rng=rng, exchange_stats=exchange_stats, walker_status=walker_status)
end

function monitor_equilibration!(walkers::Vector{<:Walker}, params::ExchangeParams, n_blocks::Int; rng=Random.GLOBAL_RNG, k_hot::Int=1, k_cold::Int=length(walkers))
    K = length(walkers)

    exchange_history = Vector{ExchangeStatus}(undef, n_blocks)
    walker_status = WalkerStatus(K; k_hot=k_hot, k_cold=k_cold)

    energies = Matrix{Float64}(undef, K, n_blocks)
    
    for block in 1:n_blocks
        exchange_stats = ExchangeStatus(params)

        equilibrate!(walkers, params; rng=rng, exchange_stats=exchange_stats, walker_status=walker_status)

        exchange_history[block] = exchange_stats

        for r in 1:K
            energies[r, block] = energy(walkers[r].replica)
        end
    end

    return energies, exchange_history, walker_status
end

function monitor_equilibration!(replicas::Vector{TR}, params::ExchangeParams, n_blocks::Int; rng=Random.GLOBAL_RNG, k_hot::Int=1, k_cold::Int=length(replicas)) where {TR <: AbstractReplica}
    walkers = make_walkers(replicas)

    return monitor_equilibration!(walkers, params, n_blocks; rng=rng, k_hot=k_hot, k_cold=k_cold)
end

struct SamplingParams
    n_samples::Int
    beta_indices::Vector{Int}
end

SamplingParams(n_samples::Int, K::Int) = SamplingParams(n_samples, collect(1:K))

function sample_replicas!(replicas::Vector{TR}, sampling_params::SamplingParams, exchange_params::ExchangeParams; rng=Random.GLOBAL_RNG) where {TR <: AbstractReplica}
    K = length(replicas)

    sampling_params.n_samples > 0 || error("n_samples must be positive.")
    check_edge_groups(exchange_params.edge_groups, K)

    sample_indices = sampling_params.beta_indices
    all(1 .<= sample_indices .<= K) || error("beta_indices out of range.")

    n_samples = sampling_params.n_samples
    n_slots = length(sample_indices)

    state0 = getstate(replicas[sample_indices[1]])
    TS = typeof(state0)

    sample_states = Matrix{TS}(undef, n_slots, n_samples)
    sampled_betas = Matrix{Float64}(undef, n_slots, n_samples)
    exchange_history = Vector{ExchangeStatus}(undef, n_samples)

    for s in 1:n_samples
        exchange_stats, _ = equilibrate!(replicas, exchange_params; rng=rng, walker_status=nothing)
        exchange_history[s] = exchange_stats

        @inbounds for (j, k) in enumerate(sample_indices)
            r = replicas[k]
            sample_states[j, s] = getstate(r)
            sampled_betas[j, s] = beta(r)
        end
    end

    return (
        states = sample_states,
        betas = sampled_betas,
        exchange_history = exchange_history,
        exchange_params = exchange_params,
        sampling_params = sampling_params,
    )
end


end # module
