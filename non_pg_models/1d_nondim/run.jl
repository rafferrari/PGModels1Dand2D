using PGModels1Dand2D
using Printf

# output directory for checkpoints, plots, and log files
out_dir = joinpath(@__DIR__, "out")
if !isdir(out_dir)
    mkdir(out_dir)
end

include("evolution.jl")
include("plotting.jl")
include("IO.jl")

################################################################################
# set up model
################################################################################

# canonical or transport-constrained case?
canonical = false

τ_A = 2e0 # nondim arrest time
τ_S = 1e2 # nondim spindown time
Ek = 1/τ_S^2 # Ekman number
S = 1/τ_A # slope Burger number
H = τ_S # depth (z ∈ [0, H] ⟹ z̃ ∈ [0, H/δ = 1/sqrt(Ek) = τ_S])
ṽ_0 = 1 # initial far-field along-slope flow
N = 1 # background stratification

# timestep
Δt̃ = minimum([τ_S/100, τ_A/100])
tSave = τ_A
α = 0.5 # time-stepping parameter

# number of grid points
nz̃ = 2^11 # good for anything at or below τ_S = 1e4

# grid (chebyshev, z̃ = 0 is bottom)
z̃ = @. H*(1 - cos(pi*(0:nz̃-1)/(nz̃-1)))/2

# bottom enhanced:
# ν0 = 1e-1
# ν1 = 1 - 1e-1
# κ0 = 1e-1
# κ1 = 1 - 1e-1
# h = 10

# not bottom enhanced:
ν0 = 1
ν1 = 0
κ0 = 1
κ1 = 0
h = 1
ν = @. ν0 + ν1*exp(-z̃/h)
κ = @. κ0 + κ1*exp(-z̃/h)

# log properties
ofile = joinpath(out_dir, "out.txt")
open(ofile, "w") do f
    write(f, "Nondimensional 1D model with Parameters:\n\n")
    write(f, @sprintf("nz̃  = %1.5e\n", nz̃))
    write(f, @sprintf("τ_A = %1.5e\n", τ_A))
    write(f, @sprintf("τ_S = %1.5e\n", τ_S))
    write(f, @sprintf("H   = %1.5e\n", H))
    write(f, @sprintf("S   = %1.5e\n", S))
    write(f, @sprintf("κ0  = %1.5e\n", κ0))
    write(f, @sprintf("κ1  = %1.5e\n", κ1))
    write(f, @sprintf("ν0  = %1.5e\n", ν0))
    write(f, @sprintf("ν1  = %1.5e\n", ν1))
    write(f, @sprintf("ṽ_0 = %1.5e\n", ṽ_0))
    write(f, @sprintf("N   = %1.5e\n", N))
    write(f, @sprintf("h   = %1.5e\n", h))
    write(f, @sprintf("Δt  = %1.5e\n", Δt̃))
    write(f, @sprintf("α   = %1.5e\n", α))
    write(f, string("\nCanonical: ", canonical, "\n"))
    write(f, @sprintf("τ_A/τ_S  = %1.5e\n", τ_A/τ_S))
end
println("Wrote '$ofile' with contents:")
open(ofile, "r") do f
    while !eof(f)
        println(readline(f))
    end
end

################################################################################
# run single integration
################################################################################

ũ, ṽ, b̃, P̃x̃ = evolve(5*τ_A)

################################################################################
# plots
################################################################################

path = ""
iSaves = 0:1:5
dfiles = [joinpath(out_dir, @sprintf("checkpoint%03d.h5", i)) for i in iSaves]
profilePlot(dfiles)