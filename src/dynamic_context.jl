using Gridap
using GridapGmsh
using Optim 

function dynamic_LMMContext(p, flow, gcontext::GridapContext; diffmethod=:finite_diff)
    normp = Operation(x-> norm(x)^p)
    absp = Operation(x -> abs(x)^p)
    DTinv = x -> transpose(dinv_flow(x, flow; diffmethod = :finite_diff))

    U = gcontext.U
    V = gcontext.V
    Ω = gcontext.Ω
    dΩ = gcontext.dΩ

    function l2_norm_grad(x)
        u = FEFunction(U, x)
        mCG = x->mean_CG_tensor(x, flow)
        return sqrt(sum( ∫(∇(u) ⋅ (mCG ⋅ ∇(u)))*dΩ ))
    end

    function lp_norm(x)
        u = FEFunction(U, x)
        return (sum( ∫(absp(u))*dΩ ))^(1/p)
    end

    function lp_norm_grad(x)
        u = FEFunction(U, x)
        return  (sum( ∫(normp(∇(u)) + normp(DTinv ⋅∇(u)))*dΩ ))^(1/p)
    end

    function J(x)
        F = lp_norm_grad(x)^p
        G = 2lp_norm(x)^p
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

        mCG = x->mean_CG_tensor(x, flow)

        # solve -Δ_dyn d = J'(u) 
        lhs(d, v) = ∫(∇(d) ⋅ (mCG ⋅ ∇(v))) * dΩ
        rhs(v) = ∫( (p / (2b^2)) * 
                   (- b *  (    η2d(∇(u))     ⋅      ∇(v) )  
                    - b * ( η2d(DTinv ⋅ ∇(u)) ⋅ (DTinv ⋅∇(v)) )
                    + a * η1d(u) * v)) * dΩ

        op = AffineFEOperator(lhs, rhs, U, V)
        ls = LUSolver()
        solver = LinearFESolver(ls)
        d = Gridap.solve(solver, op).free_values
        return d
    end
    return LMMContext(J, peak, pseudogradient, l2_norm_grad, lp_norm, lp_norm_grad, gcontext, p)
end


