using Contour
using Gridap
using LinearAlgebra

function cheeger_ratio(xs, ys)
    return perimeter(xs, ys) / area(xs, ys)
end

function area(xs, ys)
    N = length(xs)
    xs = [xs..., xs[1]]
    ys = [ys..., ys[1]]
    return abs(sum(det([xs[k] xs[k+1]; ys[k] ys[k+1]]) 
                   for k = 1:N) / 2)
end

function montecarloarea(xx, yy, zz, level)
    vol = (maximum(xx) - minimum(xx))*(maximum(yy) - minimum(yy))
    return vol * count(>(level), zz) / (length(xx) * length(yy))
end

function perimeter(xs, ys; circular=true)
    if circular
        xs = [xs..., xs[1]]
        ys = [ys..., ys[1]]
    end
    dxs = xs[2:end] - xs[1:end-1]
    dys = ys[2:end] - ys[1:end-1]
    dds = vcat.(dxs, dys)
    return sum(norm.(dds))
end

function dynamic_perimeter(xs, ys, flow;circular = true)
    pts = vcat.(xs, ys)
    Tpts = flow.(pts)
    Txs = getindex.(Tpts, 1)
    Tys = getindex.(Tpts, 2)
    return perimeter(xs, ys; circular) + perimeter(Txs, Tys; circular)
end

function best_cheeger_set(u, xx, yy; nlevels = 100, 
                          flow = nothing, 
                          smoothing = false, 
                          areamethod=:shoelace)

    if areamethod ∉ [:shoelace, :montecarlo]
        throw(ArgumentError("Unknown area method. Use :shoelace or :montecarlo"))
    end
    pts = [Point(x,y) for x in xx, y in yy]

    zs = u.(pts)
    zs ./= sign(sum(zs)) * maximum(abs.(zs))
   
    ret_xs = []
    ret_ys = []
    xs = Float64[]
    ys = Float64[]
    current_best = Inf

    lvls = Float64[]
    ratios = Float64[]
    
    for cl in levels(contours(xx, yy, zs, nlevels))
        push!(lvls, level(cl))
        total_perim = 0.0
        total_area = 0.0
        for line in lines(cl)
            xs, ys = coordinates(line)
            smoothing && (smooth!(xs), smooth!(ys))
            if isnothing(flow)
                total_perim += perimeter(xs, ys)
            else
                total_perim += dynamic_perimeter(xs, ys, flow)
            end

            if areamethod == :shoelace
                total_area += area(xs, ys)
            end
        end

        if areamethod == :montecarlo
            total_area = montecarloarea(xx,yy,zs, level(cl))
        end

        current_ratio = total_perim/(2total_area) 
        if current_ratio < current_best
            ret_xs = xs
            ret_ys = ys
            current_best = current_ratio
        end
    end
    return ret_xs, ret_ys
end

function cheeger_ratios(xx, yy, zz::Array; nlevels = 100, 
                        flow = nothing, 
                        smoothing=false, 
                        areamethod=:shoelace, 
                        circular = true)
    if areamethod ∉ [:shoelace, :montecarlo]
        throw(ArgumentError("Unknown area method. Use :shoelace or :montecarlo"))
    end

    lvls = Float64[]
    ratios = Float64[]
    zz ./= sign(sum(zz)) * maximum(abs.(zz))

    for cl in levels(contours(xx, yy, zz, nlevels))
        push!(lvls, level(cl))
        total_perim = 0.0
        total_area = 0.0
        for line in lines(cl)
            local xs, ys
            xs, ys = coordinates(line)
            smoothing && (smooth!(xs), smooth!(ys))
            if isnothing(flow)
                total_perim += perimeter(xs, ys; circular)
            else
                total_perim += dynamic_perimeter(xs, ys, flow; circular)
            end

            if areamethod == :shoelace
                total_area += area(xs, ys)
            end
        end

        if areamethod == :montecarlo
            total_area = montecarloarea(xx,yy,zz, level(cl))
        end

        push!(ratios, total_perim/(2total_area))
    end
    return lvls, ratios
end


function cheeger_ratios(u, xx, yy; nlevels = 100, 
                        flow = nothing, 
                        smoothing=false, 
                        areamethod= :shoelace)
    pts = [Point(x,y) for x in xx, y in yy]
    zs = u.(pts)
    return cheeger_ratios(xx, yy, zs; nlevels, flow, smoothing, areamethod)
end


function smooth!(vals)
    # periodic continuation
    per_vals = [vals[end], vals..., vals[1]]
    vals .= (per_vals[1:end-2] .+ per_vals[2:end-1] .+ per_vals[3:end]) ./ 3
    return vals
end



