using Gridap
using GridapGmsh
using Optim 

struct GridapContext
    U
    V
    Ω
    dΩ
end

function CartesianGridapContext(domain, partition; isperiodic=(false, false))
    model = CartesianDiscreteModel(domain, partition; isperiodic)
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)

    V = TestFESpace(model, reffe, conformity = :H1, dirichlet_tags="boundary")
    U = TrialFESpace(V)

    Ω = Triangulation(model)
    degree = 5
    dΩ = Measure(Ω, degree)
 
    return GridapContext(U, V, Ω, dΩ)
end


function GridapContext(mesh_filename)
    model = GmshDiscreteModel(mesh_filename)
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)

    V = TestFESpace(model, reffe, conformity = :H1, dirichlet_tags="border")
    U = TrialFESpace(V)

    Ω = Triangulation(model)
    degree = 5
    dΩ = Measure(Ω, degree)
 
    return GridapContext(U, V, Ω, dΩ)
end


"""
    struct LMMContext

Struct holding all necessary function definitions
for the local min-max algorithm by Yao & Zhou (2007). 

The function definitions should already be discretized, 
i.e. operate on raw arrays. 
"""
struct LMMContext{F1, F2, F3, F4, F5, F6}
    J::F1
    peak::F2
    pseudogradient::F3
    l2_norm_grad::F4
    lp_norm::F5
    lp_norm_grad::F6
    context::GridapContext
    p::Float64
end

function LMMContext(p, gcontext::GridapContext)
    normp = Operation(x-> norm(x)^p)
    absp = Operation(x -> abs(x)^p)

    U = gcontext.U
    V = gcontext.V
    Ω = gcontext.Ω
    dΩ = gcontext.dΩ

    function l2_norm_grad(x)
        u = FEFunction(U, x)
        return sqrt(sum( ∫(∇(u) ⋅ ∇(u))*dΩ ))
    end

    function lp_norm(x)
        u = FEFunction(U, x)
        return (sum( ∫(absp(u))*dΩ ))^(1/p)
    end

    function lp_norm_grad(x)
        u = FEFunction(U, x)
        return  (sum( ∫(normp(∇( u )))*dΩ ))^(1/p)
    end

    function J(x)
        F = lp_norm_grad(x)^p
        G = lp_norm(x)^p
        return F / G
    end

    function peak(xs, t)
        if length(xs) == 1
            return xs[1], [1.0]
        end

        objective(ts) = -J(sum(ts .* xs))
        manif = Optim.Sphere()
        res = optimize(objective, t, LBFGS(manifold=manif))
        ts = Optim.minimizer(res)
        return sum(ts .* xs), ts
    end

    function pseudogradient(x)
        normalize2(x) = iszero(norm(x)) ? zero(x) : x / norm(x)

        η1d = Operation(x -> abs(x) ^ (p-1) * sign(x))
        η2d = Operation(x -> norm(x) ^ (p-1) * normalize2(x))

        a = lp_norm_grad(x)^p
        b = lp_norm(x)^p

        u = FEFunction(U, x)

        # solve -Δd = J'(u) 
        lhs(d, v) = ∫(∇(d) ⋅ ∇(v)) * dΩ
        rhs(v) = ∫( (p / b^2) * (- b * (η2d(∇(u)) ⋅ ∇(v)) + a * η1d(u) * v)) * dΩ

        op = AffineFEOperator(lhs, rhs, U, V)
        ls = LUSolver()
        solver = LinearFESolver(ls)
        d = Gridap.solve(solver, op).free_values
        return d
    end
    return LMMContext(J, peak, pseudogradient, l2_norm_grad, lp_norm, lp_norm_grad, gcontext, p)
end


