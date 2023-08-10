using Gridap
using Gridap.CellData
using Statistics
using JLD2
using DynamicPLaplacian 
using Setfield
using ProgressMeter
import Base.@kwdef


@kwdef struct RunParameters 
    p::Float64
    n::Int64
    nev::Int64
    lmm_context
    domain::NTuple{4, Float64}
    isperiodic::NTuple{2, Bool}
    eval_partition::NTuple{2, Int64}
    flow
    name::String
end

function make_parameter_collection(; 
                     name, 
                     flow, 
                     nev=1,
                     nrange = 50 .* 2 .^(0:2), 
                     make_partition = n -> (n,n),
                     make_eval_part = n -> (2n,2n), 
                     prange = range(2, 1.05, 15),
                     domain = (0,1,0,1),
                     isperiodic = (false, false),
                     context_change = x -> x)
    pars = []
    for n in nrange 
        pth = "results/$name/data/n=$n"
        mkpath(pth)
        context = CartesianGridapContext(domain, 
                                         make_partition(n);
                                         isperiodic)
        for p in prange 
            lmm_context = dynamic_LMMContext(p, flow, context, 
                                       diffmethod=:forwarddiff)
            lmm_context = context_change(lmm_context)
            eval_partition = make_eval_part(n)
            par = RunParameters(; 
                                p, n, nev, lmm_context, 
                                domain, isperiodic, 
                                name, 
                                eval_partition,
                                flow)
            push!(pars, par)
        end
    end
    return pars
end

function run_with_parameters(par)
    @info "Running $(par.name) with p = $(par.p)"
    n = par.n
    p = par.p
    nev = par.nev
    flow = par.flow
    lmm_context = par.lmm_context
    eval_part = par.eval_partition
    domain = par.domain

    pth = "results/$(par.name)/data/n=$n"
    mkpath(pth)
    calc_and_save(;p, nev, 
                  flow, 
                  lmm_context, 
                  pth, 
                  eval_part, 
                  domain)
end

function calc_and_save(;p, 
                        nev, 
                        flow, 
                        lmm_context, 
                        pth, 
                        eval_part, 
                        domain)

    # calculate eigenfunctions of p-Laplacian 
    results = get_dyn_plap_eigs(nev, flow, 
                                lmm_context, 
                                verbose=true, 
                                diffmethod=:forwarddiff, 
                                max_it = 10000)

    # extract solutions
    us = [ 
          FEFunction(lmm_context.context.U,results[k].u)
          for k in eachindex(results)
         ]

    # save as vtk
    vtk_fields = [ "k=$(k)_p=$p" => us[k] for k in eachindex(us)]
    writevtk(lmm_context.context.Ω, 
             "$pth/p=$p.vtu", 
             cellfields = vtk_fields)

    # other stuff to save in jld
    Jvals = [lmm_context.J(res.u) for res in results]

    xx, yy, zzs = eval_functions(us; domain, eval_part)
    
    iterations = [r.iterations for r in results]
    l2_normGs = [r.l2_normGs for r in results]
    jldsave("$pth/p=$p.jld2";
            iterations,
            l2_normGs,
            Jvals, 
            p, xx, yy, zzs)
end

function eval_functions(us; domain, 
                        eval_part)

    xx = range(domain[1], domain[2], eval_part[1])
    yy = range(domain[3], domain[4], eval_part[2])
    pts = [Point(x,y) for x in xx, y in yy]
    zzs = []

    for u in us
        cache = Gridap.FESpaces.return_cache(u, pts[1])
        zz = map(pts) do pt
            evaluate!(cache, u, pt)
        end
        push!(zzs, zz)
    end

    return  xx, yy, zzs
end


