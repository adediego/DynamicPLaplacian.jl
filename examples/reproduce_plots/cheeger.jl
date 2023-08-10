
# Proposition 5.1. in "An introduction to the Cheeger Problem" Parini, 2011               
function cheeger_const_rect(lx, ly)                                                       
    R = cheeger_R(lx, ly)                                                                 
    return isoper_const(lx, ly, R)                                                        
end                                                                                       
                                                                                          
                                                                                          
function isoper_const(lx, ly, R)                                                          
    boundary = 2π*R + 2lx + 2ly - 8R                                                      
    area = π*R^2 + 2*R*(lx-2R) + 2*R*(ly-2R) + (lx-2R) * (ly - 2R)                        
                                                                                          
    return boundary / area                                                                
end                                                                                       
                                                                                          
function cheeger_R(lx, ly)                                                                
    a = 4 - π                                                                             
    b = -2*(lx + ly)                                                                      
    c = lx*ly                                                                             
    R = (-b - sqrt(b^2 - 4*a*c)) / (2a)                                                   
end

function get_cheeger_set_coords_rect(lx, ly)
    R = cheeger_R(lx, ly)
    xs = []
    ys = []
    push!(xs, R)
    push!(ys, 0)

    narc = 30 

    append!(xs, lx - R .+ cos.(range(3π/2, 2π, narc)) * R)
    append!(xs, lx - R .+ cos.(range(0, π/2, narc)) * R)
    append!(xs,      R .+ cos.(range(π/2, π, narc)) * R)
    append!(xs,      R .+ cos.(range(π, 3π/2, narc)) * R)

    append!(ys,      R .+ sin.(range(3π/2, 2π, narc)) * R)
    append!(ys, ly - R .+ sin.(range(0, π/2, narc)) * R)
    append!(ys, ly - R .+ sin.(range(π/2, π, narc)) * R)
    append!(ys,      R .+ sin.(range(π, 3π/2, narc)) * R)

    push!(xs, R)
    push!(ys, 0)

    return xs, ys
end

