using PGModels1Dand2D
using Printf
using LaTeXStrings  

# output directory for checkpoints, plots, and log files
out_dir = joinpath(@__DIR__, "out")
if !isdir(out_dir)
    mkdir(out_dir)
end

include("model.jl")
include("evolution.jl")
include("plotting.jl")

################################################################################
# set up model
################################################################################

# canonical or transport-constrained case?
canonical = true

τ_A = 2e0 # nondim arrest time
τ_S = 1e2 # nondim spindown time
Ek = 1/τ_S^2 # Ekman number
S = 1/τ_A # slope Burger number
H = τ_S # depth (z ∈ [0, H] ⟹ z̃ ∈ [0, H/δ = 1/sqrt(Ek) = τ_S])
v₀ = 20 # initial far-field along-slope flow
N = 1 # background stratification
r = 0 # buoyancy damping coefficient
z_max = 100 # maximum vertical coordinate for plotting

# timestep
Δt = minimum([τ_S/5e4, τ_A/5e4])

# number of grid points
nz = 2^11

# grid (chebyshev, z = 0 is bottom)
# z = @. H*(1 - cos(pi*(0:nz-1)/(nz-1)))/2
z = range(0, H, nz)

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
ν = @. ν0 + ν1*exp(-z/h)
κ = @. κ0 + κ1*exp(-z/h)

# for BT12 mixing scheme
BT12 = true
BT12_debug = true
BT12_interval = 2e5/τ_A
κ_b = 100*κ0

# store in model
model = Model(S, v₀, N, Δt, z, ν, κ; canonical)

################################################################################
# run single integration
################################################################################

u, v, b, Px = evolve(model; t_final=10*τ_A, t_save=1*τ_A)

################################################################################
# plots
################################################################################

path = ""
i_saves = 0:1:10
dfiles = [joinpath(out_dir, @sprintf("checkpoint%03d.jld2", i)) for i in i_saves]
profile_plot(dfiles)