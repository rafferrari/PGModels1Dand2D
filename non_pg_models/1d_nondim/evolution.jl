using SparseArrays
using LinearAlgebra
using JLD2

@doc raw"""
    LHS, RHS, rhs = build_linear_system(model::Model)
   

Set up the matrices and vectors for the non-dimensional 1D evolution equations.

The equations are:
```math
u_t = v - Px + S*b + (ν*u_z)_z
v_t = -u + (ν*v_z)_z
b_t = -u + [κ*(N² + b_z)]_z
```
where ``u``, ``v``, ``b`` are the non-dimensional cross-slope flow, along-slope flow, and 
buoyancy perturbation, respectively, and ``Px`` is the non-dimensional barotropic
pressure gradient. 

Boundary conditions are
```math
u = v = 0 \quad \text{at} \quad z = 0,
b_z = -N² \quad \text{at} \quad z = 0,
dz(u) = dz(v) = dz(b) = 0 \quad \text{at} \quad z = z_{top}
```
If `canonical=true`, the the barotropic pressure gradient is constant in time, equal to `v₀`. 
Otherwise, ``Px`` is set by the transport constraint ``∫ u dz = 0``.

They dynamics are determined by two non-dimensional numbers: the slope Burger number 
``S`` and background stratification ``N²``. ``ν`` and ``κ`` are 
the non-dimensional turbulent viscosity and diffusivity, respectively.

These equations are discretized using finite differences in the vertical direction, leading
to the matrix equation
```math
x_t = A*x + rhs
```
where ``x`` is the vector of variables.
We use a Crank-Nicolson time-stepping scheme,
```math
(x^(n+1) - x^n)/Δt = 1/2*[A*x^n + A*x^(n+1)] + rhs
```    
so that 
```math
LHS = I/Δt - 1/2*A,
RHS = I/Δt + 1/2*A.
```   
"""
function build_linear_system(model::Model)
    # unpack
    S = model.S
    N = model.N
    ν = model.ν
    κ = model.κ
    z = model.z
    Δt = model.Δt
    canonical = model.canonical
    v₀ = model.v₀

    # setup
    nz, ndof, imap = get_dof(z)
    LHS = Tuple{Int64,Int64,Float64}[]  
    RHS = Tuple{Int64,Int64,Float64}[]  
    rhs = zeros(ndof)

    # Main loop, insert stencil in matrices for each node point
    for j=2:nz-1
        # dz stencil
        fd_z = mkfdstencil(z[j-1:j+1], z[j], 1)
        ν_z = sum(fd_z.*ν[j-1:j+1])
        κ_z = sum(fd_z.*κ[j-1:j+1])

        # dzz stencil
        fd_zz = mkfdstencil(z[j-1:j+1], z[j], 2)

        # 1st eqtn: u_t = v - Px + S*b + (ν*u_z)_z
        row = imap[1, j]
        # time derivative
        add_to_matrix!(LHS, row, row, 1/Δt)
        add_to_matrix!(RHS, row, row, 1/Δt)
        # first term RHS
        add_to_CN_matrices!(LHS, RHS, row, imap[2, j], 1)
        # second term RHS
        add_to_CN_matrices!(LHS, RHS, row, ndof, -1)
        # third term RHS
        add_to_CN_matrices!(LHS, RHS, row, imap[3, j], S)
        # fourth term: dz(ν*dz(u))) = dz(ν)*dz(u) + ν*dzz(u)
        add_to_CN_matrices!(LHS, RHS, row, imap[1, j-1], ν_z*fd_z[1] + ν[j]*fd_zz[1])
        add_to_CN_matrices!(LHS, RHS, row, imap[1, j],   ν_z*fd_z[2] + ν[j]*fd_zz[2])
        add_to_CN_matrices!(LHS, RHS, row, imap[1, j+1], ν_z*fd_z[3] + ν[j]*fd_zz[3])

        # 2nd eqtn: v_t = -u + (ν*v_z)_z
        row = imap[2, j]
        # time derivative
        add_to_matrix!(LHS, row, row, 1/Δt)
        add_to_matrix!(RHS, row, row, 1/Δt)
        # first term:
        add_to_CN_matrices!(LHS, RHS, row, imap[1, j], -1)
        # second term: dz(ν*dz(v))) = dz(ν)*dz(v) + ν*dzz(v)
        add_to_CN_matrices!(LHS, RHS, row, imap[2, j-1], ν_z*fd_z[1] + ν[j]*fd_zz[1])
        add_to_CN_matrices!(LHS, RHS, row, imap[2, j],   ν_z*fd_z[2] + ν[j]*fd_zz[2])
        add_to_CN_matrices!(LHS, RHS, row, imap[2, j+1], ν_z*fd_z[3] + ν[j]*fd_zz[3])

        # 3rd eqtn: b_t = -u + [κ*(1 + b_z)]_z
        row = imap[3, j]
        # time derivative
        add_to_matrix!(LHS, row, row, 1/Δt)
        add_to_matrix!(RHS, row, row, 1/Δt)
        # first term
        add_to_CN_matrices!(LHS, RHS, row, imap[1, j], -1)
        # second term: dz(κ(N^2 + dz(b))) = N^2*dz(κ) + dz(κ)*dz(b) + κ*dzz(b)
        add_to_CN_matrices!(LHS, RHS, row, imap[3, j-1], (κ_z*fd_z[1] + κ[j]*fd_zz[1]))
        add_to_CN_matrices!(LHS, RHS, row, imap[3, j],   (κ_z*fd_z[2] + κ[j]*fd_zz[2]))
        add_to_CN_matrices!(LHS, RHS, row, imap[3, j+1], (κ_z*fd_z[3] + κ[j]*fd_zz[3]))
        rhs[row] = N^2*κ_z
    end

    # Boundary Conditions: Bottom
    # u = 0
    row = imap[1, 1] 
    add_to_matrix!(LHS, row, row, 1.0)
    # v = 0
    row = imap[2, 1] 
    add_to_matrix!(LHS, row, row, 1.0)
    # b_z = -N^2
    row = imap[3, 1] 
    fd_z = mkfdstencil(z[1:3], z[1], 1)
    add_to_matrix!(LHS, row, imap[3, 1], fd_z[1])
    add_to_matrix!(LHS, row, imap[3, 2], fd_z[2])
    add_to_matrix!(LHS, row, imap[3, 3], fd_z[3])
    rhs[row] = -N^2

    # Boundary Conditions: Top
    fd_z = mkfdstencil(z[nz-2:nz], z[nz], 1)
    # dz(u) = 0
    row = imap[1, nz] 
    add_to_matrix!(LHS, row, imap[1, nz-2], fd_z[1])
    add_to_matrix!(LHS, row, imap[1, nz-1], fd_z[2])
    add_to_matrix!(LHS, row, imap[1, nz],   fd_z[3])
    # dz(v) = 0
    row = imap[2, nz] 
    add_to_matrix!(LHS, row, imap[2, nz-2], fd_z[1])
    add_to_matrix!(LHS, row, imap[2, nz-1], fd_z[2])
    add_to_matrix!(LHS, row, imap[2, nz],   fd_z[3])
    # dz(b) = 0
    row = imap[3, nz]
    add_to_matrix!(LHS, row, imap[3, nz-2], fd_z[1])
    add_to_matrix!(LHS, row, imap[3, nz-1], fd_z[2])
    add_to_matrix!(LHS, row, imap[3, nz],   fd_z[3])

    # transport constraint
    row = ndof
    if canonical == true
        # canonical 1D: Px = constant in time (equal to v₀)
        add_to_matrix!(LHS, row, ndof, 1)
        rhs[row] = v₀
    else
        # transport-constrained 1D: Px such that ∫ u dz = 0
        for j=1:nz-1
            # trapezoidal rule: (u_j+1 + u_j)/2 * Δz_j
            add_to_matrix!(LHS, row, imap[1, j],   (z[j+1] - z[j])/2)
            add_to_matrix!(LHS, row, imap[1, j+1], (z[j+1] - z[j])/2)
        end
    end

    # Create CSC sparse matrices from matrix elements
    LHS = sparse((x->x[1]).(LHS), (x->x[2]).(LHS), (x->x[3]).(LHS), ndof, ndof)
    RHS = sparse((x->x[1]).(RHS), (x->x[2]).(RHS), (x->x[3]).(RHS), ndof, ndof)

    return lu(LHS), RHS, rhs
end
function get_dof(z)
    nz = length(z)
    ndof = 3*nz + 1 # number of degrees of freedom (u, v, b, Px)
    imap = reshape(1:(ndof-1), 3, nz)    
    return nz, ndof, imap
end
function add_to_matrix!(A, row, col, val)
    push!(A, (row, col, val))
end
function add_to_CN_matrices!(LHS, RHS, row, col, val)
    add_to_matrix!(LHS, row, col, -1/2*val)
    add_to_matrix!(RHS, row, col, +1/2*val)
end

"""
    u, v, b, Px = evolve(model::Model; t_final, t_save)

Solve 1D equations until `t_final`, saving every `t_save`.
"""
function evolve(model::Model; t_final, t_save)
    # unpack
    Δt = model.Δt
    z = model.z
    v₀ = model.v₀

    nz, ndof, imap = get_dof(z)

    # timestep
    n_steps = Int64(ceil(t_final/Δt))
    n_save = Int64(floor(t_save/Δt))

    # get matrices and vectors
    LHS, RHS, rhs = build_linear_system(model)

    # initial condition
    t = 0
    sol = zeros(ndof)
    # start with far-field geostrophic flow
    sol[imap[2, :]] .= v₀
    sol[ndof] = v₀

    # save initial condition
    u = sol[imap[1, :]]
    v = sol[imap[2, :]]
    b = sol[imap[3, :]]
    Px = sol[ndof]
    ofile = joinpath(out_dir, @sprintf("checkpoint%03d.jld2", 0))
    jldsave(ofile; u, v, b, Px, t, model)
    @info "Saved checkpoint to '$ofile'"
    i_save = 1

    # main loop
    for i=1:n_steps
        t += Δt

        # if model.κ changes, rebuild the linear system
        # LHS, RHS, rhs = build_linear_system(model)

        # solve linear system
        ldiv!(sol, LHS, RHS*sol + rhs)

        # log
        if i % n_save == 0
            # gather solution
            u = sol[imap[1, :]]
            v = sol[imap[2, :]]
            b = sol[imap[3, :]]
            Px = sol[ndof]

            # save data
            ofile = joinpath(out_dir, @sprintf("checkpoint%03d.jld2", i_save))
            jldsave(ofile; u, v, b, Px, t, model)
            @info "Saved checkpoint to '$ofile'"
            i_save += 1
        end
    end

    u = sol[imap[1, :]]
    v = sol[imap[2, :]]
    b = sol[imap[3, :]]
    Px = sol[ndof]
    return u, v, b, Px
end
