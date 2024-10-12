using nuPGCM, PyPlot, PyCall, Printf, SparseArrays, LinearAlgebra, HDF5, Dierckx

plt.style.use("../../plots.mplstyle")
close("all")
pygui(false)

include("params.jl")
include("evolution.jl")
include("plotting.jl")
include("utils.jl")

################################################################################
# run single integration
################################################################################

ũ, ṽ, b̃, P̃x̃ = evolve(5*τ_A)

################################################################################
# plots
################################################################################

path = ""
iSaves = 0:1:5
dfiles = string.(path, "checkpoint", iSaves, ".h5")
profilePlot(dfiles)