module PlanckLite

using Artifacts
using ChainRulesCore
using LinearAlgebra
using NPZ
using SparseArrays

function __init__()
    global data = npzread(joinpath(artifact"planck_lite", "data.npy"))
    global cov = npzread(joinpath(artifact"planck_lite", "cov.npy"))
    global blmin_TT = npzread(joinpath(artifact"planck_lite", "blmin_TT.npy"))
    global blmax_TT = npzread(joinpath(artifact"planck_lite", "blmax_TT.npy"))
    global blmin = npzread(joinpath(artifact"planck_lite", "blmin.npy"))
    global blmax = npzread(joinpath(artifact"planck_lite", "blmax.npy"))
    global bin_w_TT = npzread(joinpath(artifact"planck_lite", "bin_w_TT.npy"))
    global bin_w = npzread(joinpath(artifact"planck_lite", "bin_w.npy"))
    global plmin_TT = 30
    global plmin = 30
    global nbintt = 215
    global nbinte = 199
    global nbinee = 199
    global ellmin = 2
    global matrix_W_TT = zeros(nbintt, 2507)

    for i in 1:nbintt
        matrix_W_TT[i,blmin_TT[i]+1+plmin_TT-ellmin:blmax_TT[i]+1+plmin_TT-ellmin] .= bin_w_TT[blmin_TT[i]+1:blmax_TT[i]+1]
    end
    matrix_W_TT = sparse(matrix_W_TT)
    global matrix_W = zeros(nbinte, 1995)
    for i in 1:nbinte
        matrix_W[i,blmin[i]+1+plmin-ellmin:blmax[i]+1+plmin-ellmin] .= bin_w[blmin[i]+1:blmax[i]+1]
    end
    matrix_W = sparse(matrix_W)


    return nothing
end

function bin_Cℓ(Cltt, Clte, Clee)
    Cltt_bin = PlanckLite.matrix_W_TT * Cltt
    Clte_bin = PlanckLite.matrix_W * Clte
    Clee_bin = PlanckLite.matrix_W * Clee

    return vcat(Cltt_bin, Clte_bin, Clee_bin)
end

function ChainRulesCore.rrule(::typeof(new_bin_Cℓ), Cltt, Clte, Clee)
    Y = new_bin_Cℓ(Cltt, Clte, Clee)
    function new_bin_Cℓ_pullback(Ȳ)
        ∂Cltt = @thunk((Ȳ[1:nbintt]' * matrix_W_TT)[1,:])#, zeros(nbinte), zeros(nbinee))
        ∂Clte = @thunk((Ȳ[nbintt+1:nbintt+nbinte]' * matrix_W)[1,:])# vcat(zeros(nbintt), Ȳ[nbintt+1:nbintt+nbinte]' *matrix_W, zeros(nbinee))
        ∂Clee = @thunk((Ȳ[nbintt+nbinte+1:end]' * matrix_W)[1,:])#@thunk vcat(zeros(nbintt), zeros(nbinte), Ȳ[nbintt+nbinte+1:end]' *matrix_W)
        return NoTangent(), ∂Cltt, ∂Clte, ∂Clee
    end
    return Y, new_bin_Cℓ_pullback
end

end # module PlanckLite
