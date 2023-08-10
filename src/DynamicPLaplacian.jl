module DynamicPLaplacian 

include("context.jl")
export GridapContext, CartesianGridapContext

include("dynamic_context.jl")
export dynamic_LMMContext

include("local_min_max.jl")
export local_min_max

include("dynamic_laplacian_egs.jl")
export get_dyn_lap_eigs, get_dyn_plap_eigs

include("save_vtk.jl")
export save_vtk

include("rot_gyre.jl")
export rot_gyre_flow

include("cheeger_ratio.jl")
export cheeger_ratios, best_cheeger_set

end # module DynamicPLaplacian
