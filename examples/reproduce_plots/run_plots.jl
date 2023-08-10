@time begin
    mkpath("figures/")
    include("define_plots.jl")
    
    # make additional colormaps available
    push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usepgfplotslibrary{colormaps}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"""
          \pgfplotsset{
                /pgfplots/colormap={modifiedthermal}{rgb255=(0,0,179) 
                                                     rgb255=(255,51,0)
                                                     rgb255=(255,255,0) 
                                                     }}""")
    
    
    include("define_experiments.jl")
    
    
    include("plot_pseudograds.jl")
    include("plot_minimiser.jl")
    
    # takes a parameter generating function
    # and plots figures for them 
    function plot_the_plots(parsmaker)
        pars = parsmaker()
        plot_best_level_set(pars[[1,5,8]])
        ratios_plot(pars[[1,5,8]])
        mean_std_min_plot(pars)
        eig_func_plot.(pars[[1,5,8]])
    end
    
    parsmakers = [
                  make_static_unit_square_pars, 
                  make_std_map_pars,
                  make_rot_gyre_pars, 
                  make_cylinder_flow_pars]
    
    plot_iterations_dashed([f() for f in parsmakers])
    plot_the_plots.(parsmakers)
    @info "Time for plotting: "
end
    


