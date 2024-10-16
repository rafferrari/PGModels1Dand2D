using PGModels1Dand2D
using PyPlot
using PyCall
using Printf

pc = 1/6

gridspec = pyimport("matplotlib.gridspec")
slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)

warnings = pyimport("warnings")
warnings.filterwarnings("ignore")

plt.style.use("../plots.mplstyle")
plt.close("all")
pygui(false)

set_out_folder("../output/")

function run1d(; tc)
    # parameters
    f = -5.5e-5
    N2 = 1e-6
    nz = 2^8
    H = 2e3
    θ = 2.5e-3         
    U = [0.0]

    # use boundary layer model?
    bl = tc

    # grid
    if bl
        z = H*collect(-1:1/(nz-1):0)
    else
        z = @. -H*(cos(pi*(0:nz-1)/(nz-1)) + 1)/2 # chebyshev 
    end
    
    # diffusivity
    κ0 = 6e-5
    κ1 = 2e-3
    h = 200
    κ_func(z) = κ0 + κ1*exp(-(z + H)/h)
    κ_z_func(z) = -κ1/h*exp(-(z + H)/h)

    # viscosity
    μ = 1e0
    ν_func(z) = μ*κ_func(z)
    
    # timestepping
    Δt = 0.01*secs_in_year
    t_save = Δt
    
    # create model struct
    m = ModelSetup1D(bl, f, nz, z, H, θ, ν_func, κ_func, κ_z_func, N2, Δt, tc, U)

    # save and log params
    save_setup(m)

    # set initial state
    b = zeros(nz)
    χ = zeros(nz)
    χ, u, v = invert(m, b)
    i = [1]
    s = ModelState1D(b, χ, u, v, i)

    # solve transient
    evolve!(m, s, 3*secs_in_year, t_save) 
    
    return m, s
end

function run2d(topo; bl = false)
    # parameters
    f = -5.5e-5
    L = 2e6
    nξ = 2^8 + 1 
    nσ = 2^8
    coords = "cartesian"
    periodic = true

    # grids: even spacing in ξ and chebyshev in σ (unless bl)
    ξ = collect(0:L/nξ:(L - L/nξ))
    if bl
        σ = collect(-1:1/(nσ-1):0)
    else
        σ = @. -(cos(pi*(0:nσ-1)/(nσ-1)) + 1)/2  
    end
    
    # topography
    no_net_transport = true
    H0 = 2e3
    amp = 0.4*H0
    function H_func(x)
        if topo == "ridge"
            return H0 + amp*cos(2*π*x/L)
        elseif topo == "flat"
            return H0 + 0*x
        end
    end
    function Hx_func(x)
        if topo == "ridge"
            return -2*π/L*amp*sin(2*π*x/L)
        elseif topo == "flat"
            return 0*x
        end
    end

    # diffusivity
    κ0 = 6e-5
    κ1 = 2e-3
    h = 200
    κ_func(ξ, σ) = κ0 + κ1*exp(-H_func(ξ)*(σ + 1)/h)

    # viscosity
    μ = 1e0
    ν_func(ξ, σ) = μ*κ_func(ξ, σ)

    # stratification
    N2 = 1e-6
    N2_func(ξ, σ) = N2
    
    # timestepping
    # Δt = 10*secs_in_day
    Δt = 0.01*secs_in_year
    t_plot = 15*secs_in_year
    # t_save = 3*secs_in_year
    t_save = Δt
    
    # create model struct
    m = ModelSetup2D(bl, f, no_net_transport, L, nξ, nσ, coords, periodic, ξ, σ, H_func, Hx_func, ν_func, κ_func, N2_func, Δt)

    # save and log params
    save_setup(m)

    # set initial state
    b = zeros(nξ, nσ)
    for i=1:nξ
        b[i, :] = cumtrapz(m.N2[i, :], m.z[i, :]) .- trapz(m.N2[i, :], m.z[i, :])
    end
    χ = zeros(nξ, nσ)
    uξ = zeros(nξ, nσ)
    uη = zeros(nξ, nσ)
    uσ = zeros(nξ, nσ)
    i = [1]
    s = ModelState2D(b, χ, uξ, uη, uσ, i)

    # solve
    evolve!(m, s, 3*secs_in_year, t_plot, t_save) 

    return m, s
end

function plot_anim(folder2D, states; folder1Dcan="", folder1Dtc="")
    m2D = load_setup_2D("$folder2D/setup.h5")
    ix = argmin(abs.(m2D.ξ .- m2D.L/4))
    if folder1Dcan != ""
        m1Dcan = load_setup_1D("$folder1Dcan/setup.h5")
    end
    if folder1Dtc != ""
        m1Dtc  = load_setup_1D("$folder1Dtc/setup.h5")
    end
    for i ∈ states
        # load
        s2D = load_state_2D("$folder2D/state$i.h5")
        if folder1Dcan != ""
            s1Dcan = load_state_1D("$folder1Dcan/state$i.h5")
        end
        if folder1Dtc != ""
            s1Dtc = load_state_1D("$folder1Dtc/state$i.h5")
        end

        # figure
        fig = figure(figsize=(25*pc, 26*pc))
        gs = gridspec.GridSpec(2, 6)
        ax1 = fig.add_subplot(get(gs, (0, slice(0, 3))))
        ax2 = fig.add_subplot(get(gs, (0, slice(3, 6))))
        ax3 = fig.add_subplot(get(gs, (1, slice(0, 2))))
        ax4 = fig.add_subplot(get(gs, (1, slice(2, 4))))
        ax5 = fig.add_subplot(get(gs, (1, slice(4, 6))))
        fig.suptitle(latexstring(L"$t = $", @sprintf("%.1f years", (s2D.i[1] - 1)*m2D.Δt/secs_in_year)))

        # ridge plots
        ridge_plot(m2D, s2D, 1e3*s2D.χ,  "", L"Streamfunction $\chi$ ($\times 10^{-3}$ m$^2$ s$^{-1}$)"; ax=ax1, vext=2.0, cb_orientation="horizontal", pad=0.3)
        ridge_plot(m2D, s2D, 1e2*s2D.uη, "", L"Along-slope flow $v$ ($\times 10^{-2}$ m s$^{-1}$)"; ax=ax2, style="pcolormesh", vext=1.5, cb_orientation="horizontal", pad=0.3)
        ax1.plot([m2D.L/1e3/4, m2D.L/1e3/4], [m2D.z[ix, 1]/1e3, 0], "r-", alpha=0.5)
        ax2.plot([m2D.L/1e3/4, m2D.L/1e3/4], [m2D.z[ix, 1]/1e3, 0], "r-", alpha=0.5)
        ax2.set_ylabel("")
        ax2.set_yticklabels([])

        # profile plots
        ax3.set_ylabel(L"Vertical coordinate $z$ (km)")
        ax4.set_yticklabels([])
        ax5.set_yticklabels([])
        ax3.set_xlabel(string(L"Streamfunction $\chi$", "\n", L"($\times10^{-3}$ m$^2$ s$^{-1}$)"))
        ax4.set_xlabel(string(L"Along-slope flow $v$", "\n", L"($\times10^{-2}$ m s$^{-1}$)"))
        ax5.set_xlabel(string(L"Stratification $\partial_z b$", "\n", L"($\times10^{-6}$ s$^{-2}$)"))
        ax3.spines["left"].set_visible(false)
        ax4.spines["left"].set_visible(false)
        ax3.axvline(0, color="k", ls="-", lw=0.5)
        ax4.axvline(0, color="k", ls="-", lw=0.5)

        # 2D
        m = m2D
        s = s2D
        bz = differentiate(s.b[ix, :], m.z[ix, :])
        ax3.plot(1e3*s.χ[ix, :],   m.z[ix, :]/1e3, label="2D")
        ax4.plot(1e2*s.uη[ix, :],  m.z[ix, :]/1e3, label="2D")
        ax5.plot(1e6*bz,           m.z[ix, :]/1e3, label="2D")

        # canonical 1D
        if folder1Dcan != ""
            m = m1Dcan
            s = s1Dcan
            bz = m.N2 .+ differentiate(s.b, m.z)
            ax3.plot(1e3*s.χ,  m.z/1e3, label="Canonical 1D", zorder=0)
            ax4.plot(1e2*s.v,  m.z/1e3, label="Canonical 1D", zorder=0)
            ax5.plot(1e6*bz,   m.z/1e3, label="Canonical 1D", zorder=0)
        end

        # tc 1D BL
        if folder1Dtc != ""
            m = m1Dtc
            s = s1Dtc
            z = @. -m1Dtc.H*(cos(pi*(0:m1Dtc.nz-1)/(m1Dtc.nz-1)) - 1)/2
            χ, b = get_full_soln(m, s, z)
            δ, μ, S, q = get_BL_params(m)
            v = @. -m.f*s.χ[1]/q/m.ν[1] - tan(m.θ)/m.f*(s.b - s.b[1])
            bz = m.N2 .+ differentiate(b, z)
            ax3.plot(1e3*χ,  (z .- z[end])/1e3, c="tab:olive", ls=(0, (5, 5)), lw=1.0, label="Transport-constrained 1D BL")
            ax4.plot(1e2*v,            m.z/1e3, c="tab:olive", ls=(0, (5, 5)), lw=1.0, label="Transport-constrained 1D BL")
            ax5.plot(1e6*bz, (z .- z[end])/1e3, c="tab:olive", ls=(0, (5, 5)), lw=1.0, label="Transport-constrained 1D BL")
        end

        # legend
        if folder1Dtc != "" || folder1Dcan != ""
            ax3.legend(ncol=3, loc=(0.1, 1.01))
        end

        # limits
        ax3.set_ylim(-2, 0)
        ax4.set_ylim(-2, 0)
        ax5.set_ylim(-2, 0)
        ax3.set_xlim(-0.1, 2.0)
        ax4.set_xlim(-1.5, 1.5)
        ax5.set_xlim(0, 1.5)
        ax5.set_xticks(0:0.5:1.5)

        # adjust
        subplots_adjust(hspace=0.4, wspace=1.0)

        # save
        savefig(@sprintf("images/anim%03d.png", i))
        println(@sprintf("images/anim%03d.png", i))
    end

    # then run something like:
    #   ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" anim.mp4
end

# example: flat topography
set_out_folder("flat/")
m, s = run2d("flat")
plot_anim("flat", 0:300)

# example: ridge topography
set_out_folder("ridge2D/")
m, s = run2d("ridge")
set_out_folder("ridge1Dtc/")
m, s = run1d(tc=true)
set_out_folder("ridge1Dcan/")
m, s = run1d(tc=false)
plot_anim("ridge2D", 0:300; folder1Dtc="ridge1Dtc", folder1Dcan="ridge1Dcan")

println("Done.")
