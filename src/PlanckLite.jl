module PlanckLite

using Artifacts
using ChainRulesCore
using LinearAlgebra
using NPZ

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
    return nothing
end

function bin_Cℓ(Cltt::AbstractArray{T}, Clte, Clee) where {T}
    Cltt_bin = zeros(T, nbintt)
    Clte_bin = zeros(T, nbinte)
    Clee_bin = zeros(T, nbinee)

    for i in 1:nbintt
        Cltt_bin[i] = LinearAlgebra.dot(Cltt[blmin_TT[i]+plmin_TT-ellmin+1:blmax_TT[i]+plmin_TT+1-ellmin],bin_w_TT[blmin_TT[i]+1:blmax_TT[i]+1])
    end

    for i in 1:nbinte
        Clte_bin[i] = LinearAlgebra.dot(Clte[blmin[i]+plmin-ellmin+1:blmax[i]+plmin+1-ellmin],bin_w[blmin[i]+1:blmax[i]+1])
    end

    for i in 1:nbinee
        Clee_bin[i] = LinearAlgebra.dot(Clee[blmin[i]+plmin-ellmin+1:blmax[i]+plmin+1-ellmin],bin_w[blmin[i]+1:blmax[i]+1])
    end

    return vcat(Cltt_bin, Clte_bin, Clee_bin)
end

function bin_Cℓ(Cl::AbstractArray{T}, nbin, blmin, plmin, ellmin, blmax, bin_W) where {T}
    #Cl_bin = zeros(T, nbin)

    #for i in 1:nbin
    #    Cl_bin[i] = LinearAlgebra.dot(Cl[blmin[i]+plmin-ellmin+1:blmax[i]+plmin+1-ellmin],bin_W[blmin[i]+1:blmax[i]+1])
    #end

    return [LinearAlgebra.dot(Cl[blmin[i]+plmin-ellmin+1:blmax[i]+plmin+1-ellmin],bin_W[blmin[i]+1:blmax[i]+1]) for i in 1:nbin]
end

function new_bin_Cℓ(Cltt::AbstractArray{T}, Clte, Clee) where {T}
    Cltt_bin = bin_Cℓ(Cltt, nbintt, blmin_TT, plmin_TT, ellmin, blmax_TT, bin_w_TT)
    Clte_bin = bin_Cℓ(Clte, nbinte, blmin, plmin, ellmin, blmax, bin_w)
    Clee_bin = bin_Cℓ(Clee, nbinee, blmin, plmin, ellmin, blmax, bin_w)

    return vcat(Cltt_bin, Clte_bin, Clee_bin)
end


#function ChainRulesCore.rrule(::typeof(bin_Cℓ), Cl, nbin, blmin, plmin, ellmin, blmax, bin_W)
#    Y = bin_Cℓ(Cl, nbin, blmin, plmin, ellmin, blmax, bin_W)
#    function bin_Cℓ_pullback(Ȳ)
 #       ∂Cl_bin = zeros(nbin)#
#
#        for i in 1:nbin
#            ∂Cl_bin[i] = LinearAlgebra.dot(Ȳ[blmin[i]+plmin-ellmin+1:blmax[i]+plmin+1-ellmin],bin_W[blmin[i]+1:blmax[i]+1])
#        end
#        return NoTangent(), ∂Cl_bin, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent()
#    end#
#
#    return Y, bin_Cℓ_pullback
#end

end # module PlanckLite
