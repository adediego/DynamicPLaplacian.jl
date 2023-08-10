using SparseArrays
using Arpack

function get_dyn_lap_eigs(flow, N, lmm_context;
                          diffmethod=:finite_diff)
    U = lmm_context.context.U
    V = lmm_context.context.V
    dΩ= lmm_context.context.dΩ

    Cinv = x-> mean_CG_tensor(x, flow; diffmethod)

    a(u,v) = ∫( (Cinv ⋅ ∇(u)) ⋅ ∇(v)) * dΩ
    m(u,v) = ∫( u*v ) * dΩ

    b(u) = ∫( 0 * u) *dΩ

    op_A = AffineFEOperator(a, b, U, V)
    op_M = AffineFEOperator(m, b, U, V)

    A, M = get_matrix(op_A), get_matrix(op_M)

    egs = eigs(A, M, nev=N, which=:SM)
    return [(FEFunction(U, real.(egs[2][:,k])), real(egs[1][k])) for k =1:N]
end

function get_dyn_plap_eigs(nev, flow, lmm_context;
                           verbose = false, 
                           max_it=90, 
                           pgrad_stop = 1e-3, 
                           diffmethod = :finite_diff)

    U = lmm_context.context.U

    dyn_egs = get_dyn_lap_eigs(flow, nev, lmm_context;
                               diffmethod)
    lower_efs = []
    results = []

    for k = 1:nev
        start_u = dyn_egs[k][1]
        start_x = start_u.free_values 
        start_x /= lmm_context.lp_norm(start_x)
        
        result = local_min_max(lower_efs,  start_x, lmm_context; 
                               max_it, verbose, pgrad_stop)
        push!(lower_efs, result.u) 
        push!(results, result)
    end

    return results 
end
        
   



