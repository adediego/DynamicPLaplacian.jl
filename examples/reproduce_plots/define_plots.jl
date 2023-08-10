using PGFPlotsX
using JLD2
using ProgressMeter
using Statistics
using DynamicPLaplacian
using Gridap
using Base.Iterators
using Printf
using LaTeXStrings
using Contour

function eig_func_plot(par)
    @info "eig_func_plot $(par.name) p=$(par.p) n=$(par.n)"
    prefix = "figures/$(par.name)_eig_func"
    fig_name = "$(prefix)_n=$(par.n)_p=$(par.p).pdf"
    loadname = "results/$(par.name)/data/n=$(par.n)/p=$(par.p).jld2"
    dict = load(loadname)
    xx = collect(dict["xx"][1:5:end])
    yy = collect(dict["yy"][1:5:end])
    zzs = dict["zzs"][par.nev][1:5:end,1:5:end]

    height = 1.0 # * min(xx[end] - xx[1], yy[end] - yy[1])
    zzs ./= sign(sum(zzs)) * maximum(abs.(zzs)) / height

    if par.name == "standard_map"
        zzs = zzs[:, 1:end-1]
        #shift to put maximum in the middle
        kk = findmax(zzs)[2]
        y_shift = -kk[2] + round(Int64, size(zzs)[2]/2)
        zzs = circshift(zzs, (0, y_shift))
        zzs = hcat(zzs, zzs[:,1])
    end

    # optics for eigenfunction plots
    ratio = par.name == "standard_map" ||
            par.name == "cylinder_flow" ?
            "1 1 3" :
            "5 5 2" 

    xtick = Dict("standard_map" => [0,π,2π], 
                 "cylinder_flow" => [0,π, 2π], 
                 "static_unit_square" => [0, 0.5, 1], 
                 "rot_gyre" => [0,0.5, 1])[par.name]

    xticklabels = Dict("standard_map" =>       ["0",L"\pi",L"2\pi"], 
                       "cylinder_flow" =>      ["0",L"\pi",L"2\pi"],
                       "static_unit_square" => ["0", L"\frac12", "1"], 
                       "rot_gyre" =>           ["0",L"\frac12", "1"])[par.name]

    ytick = Dict("standard_map" => [0,π, 2π], 
                 "cylinder_flow" => [0,π], 
                 "static_unit_square" => [0, 0.5, 1], 
                 "rot_gyre" => [0,0.5, 1])[par.name]

    yticklabels = Dict("standard_map" =>       ["0",L"\pi",L"2\pi"], 
                       "cylinder_flow" =>      ["0",L"\pi"],
                       "static_unit_square" => ["0",L"\frac12", "1"], 
                       "rot_gyre" =>           ["0",L"\frac12", "1"])[par.name]

    xlims = Dict("standard_map" =>       [0, 2π], 
                 "cylinder_flow" =>      [0, 2π], 
                 "static_unit_square" => [0, 1], 
                 "rot_gyre" =>           [0, 1])[par.name]

    ylims = Dict("standard_map" =>       [0, 2π], 
                 "cylinder_flow" =>      [0, π], 
                 "static_unit_square" => [0, 1], 
                 "rot_gyre" =>           [0, 1])[par.name]


    if par.name == "standard_map"
        height = "7cm"
    else
        height = "5cm"
    end
    
    figure = @pgf Axis(
        {
          "unit vector ratio" = ratio,
          "unit rescale keep size" = false, 
          "colormap name" = "viridis", 
          "colorbar",
          xmin = xlims[1], 
          xmax = xlims[2], 
          ymin = ylims[1], 
          ymax = ylims[2], 
          zmin = minimum(zzs), 
          zmax = maximum(zzs),
          xtick = xtick, 
          ytick = ytick, 
          xticklabels = xticklabels, 
          yticklabels = yticklabels, 
          height = height
        }, 
        Plot3(
          { surf, 
            shader = "faceted", 
            "faceted color" = "black"
          }, 
          Table(xx, yy, zzs)
         )
       )

    pgfsave(fig_name, figure) 

    plot_height = (par.name == "cylinder_flow") ? "120pt" : "240pt"

    height = "5cm"
    width  = "5cm"
    if par.name == "cylinder_flow"
        height = "3cm"
    end

    figure2 = @pgf Axis(
        {
          "axis equal", 
          width="240pt", 
          height=plot_height, 
          "colormap name" = "viridis", 
          view = (0,90), 
          xmin = xlims[1], 
          xmax = xlims[2], 
          ymin = ylims[1], 
          ymax = ylims[2], 
          xtick = xtick, 
          ytick = ytick, 
          xticklabels = xticklabels, 
          yticklabels = yticklabels, 
          height=height, 
          width=width
        }, 
        Plot3(
              { "contour lua" = {number = 17, labels = false}
              }, 
          Table(xx, yy, zzs)
         )
       )

    fig_name2 = "$(prefix)_contour_n=$(par.n)_p=$(par.p).pdf"
    pgfsave(fig_name2, figure2) 
end


function mean_std_min_plot(pars)
    prefix = "figures/$(pars[1].name)"
    fig_name = "$(prefix)_n=$(pars[1].n)_statistics.pdf"

    ps = [par.p for par in pars]
    @info "Getting isoperimetry statistics"
    stats = @showprogress map(get_isoper_statistics, pars)

    median_isop = [s[2] for s in stats]
    min_isop =    [s[3] for s in stats]
    ymin = 0.9minimum(min_isop)
    figure = @pgf Axis(
              {
               xlabel=L"$p$", 
               "legend pos" = "north west", 
               "axis y discontinuity" = "crunch",
               height="6cm",
               xmin = 1.0, 
               ymin = ymin
              },
              PlotInc({mark="*"},Coordinates(zip(ps, median_isop))), 
              LegendEntry("Median"),
              PlotInc({mark="square*"},Coordinates(zip(ps, min_isop))),
              LegendEntry("Minimum"))

    if pars[1].name == "static_unit_square"
        push!(figure, @pgf HLine({"black", "dashed"}, (4-π)/(2-sqrt(π))))
        push!(figure, raw"\addlegendimage{black, dashed, mark=none}")
        push!(figure, @pgf LegendEntry(L"h(M)"))
    end

    pgfsave(fig_name, figure)

end

function ratios_plot(pars)
    prefix = "figures/$(pars[1].name)"
    fig_name = "$(prefix)_n=$(pars[1].n)_ratios.pdf"
    data = get_ratios.(pars)
    coords = [zip(d[1], d[2]) for d in data]
    lstyles = @pgf [{"dotted", thick}, {"dashdotted", thick}, 
                    {"solid", thick}]
    plots = [ PlotInc(style, Coordinates(coord))
             for (style, coord) in zip(lstyles,coords)]
    pstr(p) = @sprintf "%.1f" p
    entries = [ LegendEntry(
                  L"$p=%$(pstr(par.p))$"
               ) for par in pars]
    pairs = flatten(zip(plots, entries))
    ymax = 10
    legpos = "south east"
    if pars[1].name == "static_unit_square"
        legpos = "north west"
    end
    if pars[1].name == "standard_map"
        legpos = "north east"
        ymax = 4
    end
    if pars[1].name == "rot_gyre"
        ymax = 13
    end
    figure = @pgf Axis(
          { 
            xlabel="Level", 
            ylabel="Ratio", 
            "no markers", 
            "legend pos" = legpos, 
            "axis y discontinuity" = "crunch", 
            height="6cm",
            ymax=ymax
          },
          pairs...)
    pgfsave(fig_name, figure)
    if pars[1].name == "static_unit_square"
        figure = @pgf Axis(
              { 
                xlabel="Level", 
                xmax = 0.2, 
                ylabel="Ratio", 
                xtick = [0,0.1,0.2],
                xticklabels = ["0","0.1","0.2"], 
                "no markers", 
                "legend pos" = legpos, 
                "axis y discontinuity" = "crunch", 
                height="6cm"
              },
              pairs...)
        fig_name = "$(prefix)_n=$(pars[1].n)_ratios_closeup.pdf"
        pgfsave(fig_name, figure)

        open("$(prefix)_ratios_plot_caption", "w") do file
            for (par, coord) in zip(pars, coords)
                cs = collect(coord)
                min_k = findmin(getindex.(cs, 2))[2]
                write(file, "p = " * pstr(par.p) * "\n")
                write(file, @sprintf "level = %.4f \n" cs[min_k][1])
                write(file, @sprintf "ratio = %.4f \n" cs[min_k][2])
            end
        end
    end
end


function get_ratios(par)
    loadname = "results/$(par.name)/data/n=$(par.n)/p=$(par.p).jld2"
    dict = load(loadname)
    xx = dict["xx"]
    yy = dict["yy"]
    zz = dict["zzs"][end]

    # for the standard map we are calculating
    # the second eigenfunction, so we dont want
    # to do the following normalizations.
    if par.name != "standard_map"
        # get rid of sign ambiguity
        zz ./= sign(sum(zz))
        # get rid of small negative overshoots
        # in the numerical approximation of the
        # first eigenfunction (they lead to 
        # spurious sets with small Cheeger ratio
        # as their superlevelsets are almost everything). 
        min_z = minimum(zz)
        zz .= max.(0, zz)
    end

    flow = x->par.flow(x,1)
    areamethod = any(par.isperiodic) ? :montecarlo : :shoelace
    circular = !any(par.isperiodic)
    lvls, ratios = cheeger_ratios(xx,yy,zz;flow,areamethod,
                                  circular, nlevels=200)
    if all(par.isperiodic)
        c_ratios = complement_ratios(lvls, ratios, xx,yy,zz)
        ratios = max.(ratios, c_ratios)
    end
    return lvls, ratios
end

include("cheeger.jl")

function plot_best_level_set(pars)
    xmin = pars[1].domain[1]
    xmax = pars[1].domain[2]
    ymin = pars[1].domain[3]
    ymax = pars[1].domain[4]

    plot_options_table = Dict( 
        "standard_map" => ( 
            xtick = [0,π,2π],
            ytick = [0,π, 2π], 
            xticklabels = ["0",L"\pi",L"2\pi"], 
            yticklabels = ["0",L"\pi",L"2\pi"], 
            plot_width = "6cm",
            legend_position = "north east", 
            inset = false, 
            inset_options = nothing
           ),
        "cylinder_flow" => (
            xtick = [0,π,2π], 
            ytick = [0,π, 2π], 
            xticklabels = ["0",L"\pi",L"2\pi"], 
            yticklabels = ["0",L"\pi",L"2\pi"], 
            legend_position = "north east", 
            plot_width = "12cm",
            inset = true, 
            inset_options = (
                outline = true, 
                center = (2.0, 0.5), 
                size   = 0.6
            )
           ), 
        "static_unit_square" => (
            xtick = [0, 0.5, 1], 
            ytick = [0, 0.5, 1], 
            xticklabels = ["0", L"\frac12", "1"], 
            yticklabels = ["0",L"\frac12", "1"], 
            legend_position = "north east", 
            plot_width = "6cm",
            inset = true, 
            inset_options = (
                outline = false, 
                center = (0.1, 0.1), 
                size   = 0.2
            )
           ),

        "rot_gyre" => (
            xtick = [0,0.5, 1],
            ytick = [0,0.5, 1], 
            xticklabels = ["0",L"\frac12", "1"], 
            yticklabels = ["0",L"\frac12", "1"], 
            legend_position = "north west", 
            plot_width = "6cm", 
            inset = true, 
            inset_options = (
                outline = true, 
                center = (0.75, 0.25), 
                size   = 0.5
            )
           )
       ) # Dict

    options = plot_options_table[pars[1].name]

    figure = @pgf Axis({
                  "axis equal", 
                  width= options.plot_width, 
                  height= "6cm", 
                  xmin  = xmin, 
                  xmax  = xmax, 
                  ymin  = ymin, 
                  ymax  = ymax, 
                  xtick = options.xtick, 
                  ytick = options.ytick, 
                  xticklabels = options.xticklabels, 
                  yticklabels = options.yticklabels, 
                  "legend pos" = options.legend_position})

    inset_square(k) = (options.inset_options.center[1] + 
                       (-1,1,1,-1)[k] * options.inset_options.size /2 , 
                       options.inset_options.center[2] + 
                       (-1,-1,1,1)[k] * options.inset_options.size /2 )

    if options.inset
        inset = @pgf Axis({
                     "axis equal", 
                      width = "6cm", 
                      height= "6cm", 
                      xtick = [ inset_square(1)[1],
                                options.inset_options.center[1],
                                inset_square(2)[1]],
                      ytick = [ inset_square(2)[2],
                                options.inset_options.center[2],
                                inset_square(3)[2]],
                      xmin = inset_square(1)[1], 
                      xmax = inset_square(2)[1], 
                      ymin = inset_square(2)[2], 
                      ymax = inset_square(3)[2]}) 
    end


    for par in pars
        lvls, ratios = get_ratios(par)
        _,k = findmin(ratios)
        lvl = lvls[k]
        xs_collection, ys_collection = get_level_set(par, lvl)
        pl_xs = []
        pl_ys = []

        for (xs,ys) in zip(xs_collection, ys_collection)
            append!(pl_xs, xs)
            append!(pl_ys, ys)
            push!(pl_xs, NaN)
            push!(pl_ys, NaN)
        end
        pl = @pgf PlotInc(
                          {mark="none"}, Coordinates(pl_xs, pl_ys))

   
        push!(figure, pl)

        pstr(p) = @sprintf "%.1f" p
        push!(figure, LegendEntry(
               L"$p=%$(pstr(par.p))$"
              )
             )

        if options.inset
            inset_plot = @pgf PlotInc(
                      {mark="none"}, Coordinates(pl_xs, pl_ys))
            push!(inset, inset_plot)
            push!(inset, LegendEntry(L"$p=%$(pstr(par.p))$"))
        end
    end

    if pars[1].name == "static_unit_square"
        pl_xs, pl_ys = get_cheeger_set_coords_rect(1,1)
        pl = @pgf PlotInc(
            {mark="none", "dashed"}, Coordinates(pl_xs, pl_ys))

        pl_close = @pgf PlotInc(
                      {mark="none", "dashed"}, Coordinates(pl_xs, pl_ys))

        push!(figure, pl)
        push!(figure, LegendEntry("Cheeger set"))
        push!(inset, pl_close)
        push!(inset, LegendEntry("Cheeger set"))
    end

    if options.inset && options.inset_options.outline
        closeup_boundary = 
          @pgf PlotInc({"dashed", mark="none"},
                       Coordinates(inset_square.([1,2,3,4,1])))
        push!(figure, closeup_boundary)
    end

    prefix = "figures/$(pars[1].name)"
    fig_name = "$(prefix)_best_set.pdf"
    save(fig_name, figure)

    if options.inset
        save("$(prefix)_best_set_closeup.pdf", inset)
    end
end

function get_cheeger_set_coords_rect(domain)
end


function get_level_set(par, lvl)
    loadname = "results/$(par.name)/data/n=$(par.n)/p=$(par.p).jld2"
    dict = load(loadname)
    xx = dict["xx"]
    yy = dict["yy"]
    zz = dict["zzs"][end]
    zz ./= sign(sum(zz)) * maximum(abs.(zz))

    if par.name == "standard_map"
        # for periodic boundary conditions 
        # there is one redundant
        # value at the boundary that would
        # lead to kinks in the plot 
        # when shifiting
        xx = xx[1:end-1]
        yy = yy[1:end-1]
        zz = zz[1:end-1,1:end-1]

        #shift to put maximum in the middle
        kk = findmax(zz)[2]
        y_shift = -kk[2] + round(Int64, size(zz)[2]/2)
        zz = circshift(zz, (0, y_shift))
    end

    cl = first(levels(contours(xx,yy,zz,  [lvl])))

    xs_collection = []
    ys_collection = []
    for line in lines(cl)
        xs, ys = coordinates(line)
        push!(xs_collection, xs)
        push!(ys_collection, ys)
    end
    return xs_collection, ys_collection
end

 
function get_isoper_statistics(par)
    lvls, ratios = get_ratios(par)
    return Statistics.mean(ratios), median(ratios), minimum(ratios)
end

function complement_ratios(lvls, ratios, xx, yy, zz)
    ret = Float64[]
    for (lvl, ratio) in zip(lvls, ratios)
        perim = ratio * DynamicPLaplacian.montecarloarea(xx,yy,zz,lvl)
        otherratio = perim / DynamicPLaplacian.montecarloarea(xx,yy,-zz,-lvl)
        push!(ret, otherratio)
    end
    return ret
end


function plot_iterations_dashed(par_arrays)
    figure = @pgf Axis(
              {
               xlabel=L"$p$", 
               ylabel="Iterations", 
               xmin=1, 
               xmax=2, 
               xtick = [1.0, range(1.3, 2.0, length=8)], 
               ytick = [10, 100, 1000, 10000], 
               yticklabels=["10", "100", raw"1\,000", raw"10\,000"], 
               ymax = 10000, 
               height = "7cm", 
               ymode = "log", 
               "scale only axis", 
               "axis y line" = "left", 
               "legend pos" = "outer north east", 
              }
          )
    legendentries = Dict( "rot_gyre" => "Transitory double gyre",
                          "static_unit_square" => "Static unit square",
                          "cylinder_flow" => "Cylinder flow", 
                          "standard_map" => "Standard map")
    for pars in par_arrays
        ps = []
        its = []
        for par in pars
            push!(ps, par.p)
    
            loadname = "results/$(par.name)/data/n=$(par.n)/p=$(par.p).jld2"
            dict = load(loadname)
            iters = dict["iterations"][end]
            push!(its, iters)
        end
    
        itercoords =  Coordinates(zip(ps, its))
        plot = @pgf PlotInc({ "only marks"}, itercoords)
        push!(figure, plot)
        push!(figure, LegendEntry(legendentries[pars[1].name]))
    end
  
#    prefix = "figures/$(pars[1].name)"
    pgfsave("figures/iterations.pdf", figure)
end

