include("define_experiments.jl")

parsmakers = [make_cylinder_flow_pars,
              make_rot_gyre_pars,
              make_static_unit_square_pars,
              make_std_map_pars]

for pm in parsmakers
    pars = pm()
    map(run_with_parameters, pars)
end
