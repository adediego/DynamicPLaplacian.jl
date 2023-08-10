function save_functions_to_vkt(name, us, lmm_context)
    writevtk(lmm_context.context.Ω, 
             name, 
             cellfields = ["$k" => us[k] for k in eachindex(us)])

end
