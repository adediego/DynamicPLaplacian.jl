t = time()
try
    include("run_experiments.jl")
    include("run_plots.jl")
catch e
    @error "Exception after $(time()-t) seconds"
    rethrow(e)
end
@info "Total time: $(time() - t)s"
