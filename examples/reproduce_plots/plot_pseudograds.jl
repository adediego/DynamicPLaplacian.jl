using DynamicPLaplacian
using Gridap
using ProgressMeter
using PGFPlotsX
using LaTeXStrings
using Statistics

@info "Pseudogradient plot"

name = "square_rot_gyre"
p = 1.83

domain = (0,1,0,1)
partition = (100,100)

context = CartesianGridapContext(domain, 
                                 partition)

lmm_context = dynamic_LMMContext(p, rot_gyre_flow, context)
U = lmm_context.context.U

dyn_egs = get_dyn_lap_eigs(rot_gyre_flow, 1, lmm_context)
start_u = dyn_egs[1][1]
start_x = start_u.free_values 
start_x /= lmm_context.lp_norm(start_x)


# hardcoded, just to illustrate that straightforward
# calculation of pseudogradient leads to irregular functions
function pseudogradient_classical(x)
    normalize2(x) = iszero(norm(x)) ? zero(x) : x / norm(x)
    dΩ = lmm_context.context.dΩ 
    V = lmm_context.context.V
    DTinv = x -> 
        transpose(DynamicPLaplacian.dinv_flow(x, rot_gyre_flow))

    function lp_norm_grad(x)
        u = FEFunction(U, x)
        return  (sum( ∫(normp(∇(u)) + normp(DTinv ⋅∇(u)))*dΩ ))^(1/p)
    end

    function lp_norm(x)
        u = FEFunction(U, x)
        return (sum( ∫(absp(u))*dΩ ))^(1/p)
    end


    normp = Operation(x-> norm(x)^p)
    absp = Operation(x -> abs(x)^p)
    η1d = Operation(x -> abs(x) ^ (p-1) * sign(x))
    η2d = Operation(x -> norm(x) ^ (p-1) * normalize2(x))

    a = lp_norm_grad(x)^p
    b = lp_norm(x)^p

    u = FEFunction(U, x)

    # solve -Δ_dyn d = J'(u) 
    lhs(d, v) = ∫(∇(d) ⋅ ∇(v)) * dΩ
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

w = lmm_context.pseudogradient(start_x)
w /= lmm_context.lp_norm(w)
w_func = FEFunction(lmm_context.context.U, w)

w2 = pseudogradient_classical(start_x)
w2 /= lmm_context.lp_norm(w2)
w_func2 = FEFunction(lmm_context.context.U, w2)

xx = range(domain[1], domain[2], 200)
yy = range(domain[3], domain[4], 200)
pts = [Point(x,y) for x in xx, y in yy]

cache = Gridap.FESpaces.return_cache(w_func, pts[1])
zz = map(pts) do pt
    evaluate!(cache, w_func, pt)
end

cache = Gridap.FESpaces.return_cache(w_func2, pts[1])
zz2 = map(pts) do pt
    evaluate!(cache, w_func2, pt)
end


# make sign deterministic (eigenfunction
# might flip sign.)
zz .*= sign((sum(zz)))


figure1 = @pgf Axis(
    {
      "axis equal", 
      width="207pt", 
      height="207pt", 
      "colormap name" = "viridis",
      title = L"\mathop{grad}^DJ^D(u)",
      view = (0,90), 
      xmin = minimum(xx), 
      xmax = maximum(xx), 
      ymin = minimum(yy), 
      ymax = maximum(yy), 
      colorbar
    }, 
    Plot3(
          { "contour lua" = {number = 34, labels = false}
          }, 
      Table(xx, yy, zz)
     )
   )

fig_name1 = "figures/pseudograd_dynamic.pdf"
pgfsave(fig_name1, figure1) 

figure2 = @pgf Axis(
    {
      "axis equal", 
      width="207pt", 
      height="207pt", 
      "colormap name" = "viridis",
      title = L"\log(|\mathop{grad}J^D(u)| + 0.1)",
      view = (0,90), 
      xmin = minimum(xx), 
      xmax = maximum(xx), 
      ymin = minimum(yy), 
      ymax = maximum(yy),
      colorbar
    }, 
    Plot3(
          { "contour lua" = {number = 34, labels = false}
          }, 
          Table(xx, yy, log.(abs.(zz2) .+ 0.1))
     )
   )


fig_name2 = "figures/pseudograd_classic.pdf"
pgfsave(fig_name2, figure2) 


