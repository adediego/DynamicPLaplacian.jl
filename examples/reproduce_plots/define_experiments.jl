using Setfield
using Statistics
using DifferentialEquations
using DynamicPLaplacian 
using Memoization


# functions to generate parameters for
# experiments

include("experiment_utils.jl")

function make_std_map_pars()

    name = "standard_map"
    function std_map(x, T) 
        a = 0.971635
        return VectorValue(x[1] + x[2] + a * sin(x[1]), 
                           x[2] + a * sin(x[1]) )
    end
    
    # rudimentary orthogonal projection 
    function change_context(lmm_context)
        new_pg = make_zero_mean_pseudogradient( lmm_context, 
                                               std_map)
        return @set lmm_context.pseudogradient = new_pg
    end
    pars = make_parameter_collection(; 
                     name, 
                     flow = std_map, 
                     nev=2,
                     nrange = [100], 
                     make_partition = n -> (n,n),
                     make_eval_part = n -> (2n+1,2n+1), 
                     prange = range(2, 1.3, 8),
                     domain = (0,2π,0,2π),
                     isperiodic = (true, true),
                     context_change = change_context)
    return pars
end

function make_static_unit_square_pars()
    name = "static_unit_square"
    pars = make_parameter_collection(; 
                     name, 
                     flow = (x,T) -> x, 
                     nev=1,
                     nrange = [100], 
                     make_partition = n -> (n,n),
                     make_eval_part = n -> (2n+1,2n+1), 
                     prange = range(2, 1.3, 8),
                     domain = (0,1,0,1),
                     isperiodic = (false, false))
    return pars
end



function make_rot_gyre_pars()
    name = "rot_gyre"
    pars = make_parameter_collection(; 
                     name, 
                     flow = rot_gyre_flow, 
                     nev=1,
                     nrange = [100], 
                     make_partition = n -> (n,n),
                     make_eval_part = n -> (2n+1,2n+1), 
                     prange = range(2, 1.3, 8),
                     domain = (0,1,0,1),
                     isperiodic = (false, false))
    return pars
end

function make_cylinder_flow_pars()
    name = "cylinder_flow"
    pars = make_parameter_collection(; 
                     name, 
                     flow = cyl_flow, 
                     nev = 1,
                     nrange = [100], 
                     make_partition = n -> (2n,n),
                     make_eval_part = n -> (4n+1,2n+1), 
                     prange = range(2, 1.3, 8),
                     domain = (0,2π,0,π),
                     isperiodic = (true, false))
    return pars
end

function cyl_flow_field!(dz, z, p, t)
    x = z[1]
    y = z[2]

    c = 0.5
    ν = 0.25
    ε = 0.25

    A(t) = 1 + 0.125sin(2sqrt(5) * t)
    G(ψ) = 1/(ψ^2 + 1)^2
    g(x,y,t) = sin(x-ν*t) * sin(y) + y/2 - π/4

    dz[1] = c - A(t)*sin(x-ν*t)*cos(y) + ε*G(g(x,y,t))*sin(t/2)
    dz[2] = A(t) * cos(x-ν*t) * sin(y)
    return dz
end


@memoize function cyl_flow(u0, T)
    tfact = 40
    tspan = (0.0, T * tfact)
    u0_array =  [u0[1], u0[2]]
    prob = ODEProblem(cyl_flow_field!, u0_array, tspan)
    sol_array = DifferentialEquations.solve(prob, Tsit5(), abstol=1e-7, reltol=1e-7)(T * tfact)
end

function make_zero_mean_pseudogradient(lmm_context, flow)

    DTinv = x -> transpose(DynamicPLaplacian.dinv_flow(x, flow; diffmethod = :finite_diff))

    p = lmm_context.p
    gcontext = lmm_context.context
    order = 1 
    reffe = ReferenceFE(lagrangian, Float64, order)
    model = gcontext.Ω.model
    V = TestFESpace(model, reffe, conformity=:H1, 
                    dirichlet_tags = "boundary", 
                    constraint=:zeromean)

    U = TrialFESpace(V)
    Ω = gcontext.Ω
    dΩ = gcontext.dΩ

    lp_norm_grad = lmm_context.lp_norm_grad
    lp_norm = lmm_context.lp_norm

    function pseudogradient(x)
        normalize2(x) = iszero(norm(x)) ? zero(x) : x / norm(x)

        η1d = Operation(x -> abs(x) ^ (p-1) * sign(x))
        η2d = Operation(x -> norm(x) ^ (p-1) * normalize2(x))

        a = lp_norm_grad(x)^p
        b = lp_norm(x)^p

        u = FEFunction(lmm_context.context.U, x)

        mCG = x->DynamicPLaplacian.mean_CG_tensor(x, flow)

        # solve -Δ_dyn d = J'(u) 
        lhs(d, v) = ∫(∇(d) ⋅ (mCG ⋅ ∇(v))) * dΩ
        rhs(v) = ∫( (p / (2b^2)) * 
                   (- b *  (    η2d(∇(u))     ⋅      ∇(v) )  
                    - b * ( η2d(DTinv ⋅ ∇(u)) ⋅ (DTinv ⋅∇(v)) )
                    + a * η1d(u) * v)) * dΩ

        op = AffineFEOperator(lhs, rhs, U, V)
        ls = LUSolver()
        solver = LinearFESolver(ls)
        d = Gridap.solve(solver, op)
        if_d = Interpolable(d)
        d_interp = interpolate(d, lmm_context.context.U)
        return d_interp.free_values
    end
    return pseudogradient
end
