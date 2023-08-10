using ProgressMeter

mutable struct LMMResult
    Js::Vector{Float64}
    l2_normGs::Vector{Float64}
    ss::Vector{Float64}
    iterations::Int64
    context
    u
end

import Base.show
function show(io::IO, ::MIME"text/plain", res::LMMResult)
    if get(io, :compact, true)
        print(io, "LMMResult(its=$(res.iterations))")
    else
        println(io, "LMMResult:")
        println(io, "   Iterations: $(res.iterations)")
        println(io, "   J value:    $(res.Js[end])")
    end
end

"""
    local_min_max(us, v1, context; kwargs...)

Run the local min-max algorithm by Yao & Zhou (2007). 
- `us`: a vector of already computed critical points
- `vs`: initial value
"""
function local_min_max(us, v1, context::LMMContext; max_it=70, 
                       verbose=false, pgrad_stop=1e-3)
    # unpack functions from LMMContext
    J = context.J
    peak = context.peak
    pseudogradient = context.pseudogradient
    l2_norm_grad = context.l2_norm_grad
    lp_norm = context.lp_norm

    # logging containers
    Js, l2_normGs, ss = Float64[], Float64[], Float64[]
    prog = ProgressUnknown("Local min-max ($(length(us) + 1))...", 
                            color = :white, spinner=true)

    # arbitrary choice of λ > 0 (not specified in paper)
    λ = 1
    v = v1 / lp_norm(v1)
    t = [zeros(length(us)); 1]

    # step 1: find peak
    u, t = peak([us; [v]], t)

    for it = 1:max_it
        Ju = J(u)

        # step 2: compute descent direction
        G = pseudogradient(u)
        w = sign(t[end]) * G # "no projection needed"

        # step 3: stop if ∇J is small (we take G instead)
        if lp_norm(G) < pgrad_stop 
            verbose && finish!(prog)
            push!(Js, Ju); push!(l2_normGs, lp_norm(G)); push!(ss, NaN)
            return LMMResult(Js, l2_normGs, ss, it, context, u)
        end


        # step 4: Armijo-like step size control
        m = max(1, ceil(log2(lp_norm(G))))
        while true
            s = λ / 2^m
            new_v = (v + s * w) / lp_norm(v + s * w) 
            new_u, new_t = peak([us; [new_v]], t)
            ΔJ = J(new_u) - Ju
            if ΔJ <= - s * abs(new_t[end]) * l2_norm_grad(G) ^ 2 / 4
                # step 5
                v = new_v
                u, t = new_u, new_t

                # logging
                push!(Js, Ju); push!(l2_normGs, lp_norm(G)); push!(ss, s)
                break
            end
            m += 1
            # step became too small
            if m > 54
                verbose && finish!(prog)
                return LMMResult(Js, l2_normGs, ss, it, context, u)
            end
            verbose && update_progress_bar!(prog, length(us) + 1, 
                                            J(u), lp_norm(G), s, it)
         end
    end
    verbose && finish!(prog)
    return LMMResult(Js, l2_normGs, ss, max_it, context, u)
end

function update_progress_bar!(prog, nev, J, G, s, it)
    ProgressMeter.update!(prog, showvalues = () ->
                          [ (:ev_count, nev), 
                            (:iterations, it), 
                            (:J, J), 
                            (:G, G),
                            (:s, s)], 
                          valuecolor = :white)
    ProgressMeter.next!(prog)
end



