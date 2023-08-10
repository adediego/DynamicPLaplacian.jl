using PGFPlotsX
include("cheeger.jl")
using DynamicPLaplacian
using LinearAlgebra
using Contour

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


function first_lap_eig(x)
    return 0.7sin(π * x[1]) * sin(π * x[2])
end

function incheeger(pt) 
    x = pt[1]
    y = pt[2]

    R = cheeger_R(1.0, 1.0)
    if R < x < 1 - R
        return  0 < y < 1
    end

    if R < y < 1 - R
        return 0 < x < 1
    end
    
    center_pts = [ [R,   R], 
                   [1-R, R], 
                   [R,   1-R], 
                   [1-R, 1-R] ]

    dists = norm.( [[x,y]] .- center_pts)
    return any(<(R), dists)
end

xx = range(0, 1, length=40)
yy = range(0, 1, length=40)
pts = [[x,y] for x in xx, y in yy]
zzs1 = first_lap_eig.(pts)

xxx = range(0, 1, length=100)
yyy = range(0, 1, length=100)
ppts = [[x,y] for x in xxx, y in yyy]
zzz = first_lap_eig.(ppts)

lvls, ratios = cheeger_ratios(xxx, yyy, zzz)
best_k = findmin(ratios)[2]
lvl = lvls[best_k]

cl = first(levels(contours(xxx,yyy,zzz,  [lvl])))

linexs = []
lineys = []
for line in lines(cl)
    xs, ys = coordinates(line)
    append!(linexs, xs)
    append!(lineys, ys)
end

figure = @pgf Axis(
    {
      "unit vector ratio" = [1,1], 
      "unit rescale keep size" = false, 
      xmin = minimum(xx), 
      xmax = maximum(xx), 
      ymin = minimum(yy), 
      ymax = maximum(yy), 
      zmin = minimum(zzs1), 
      zmax = maximum(zzs1) 
    }, 
    Plot3(
      { surf, 
        "colormap name" = "viridis",
        shader = "faceted", 
        "faceted color" = "black"
      }, 
      Table(xx, yy, zzs1)
     )
   )

fig_name = "figures/p2.pdf"
pgfsave(fig_name, figure) 

figure2 = @pgf Axis(
    {
      "axis equal", 
      width="207pt", 
      height="207pt", 
      view = (0,90), 
      xmin = minimum(xx), 
      xmax = maximum(xx), 
      ymin = minimum(yy), 
      ymax = maximum(yy) 
    }, 
    Plot3(
          { "contour lua" = {number = 17, labels = false},
            "colormap name" ="viridis" 
          }, 
      Table(xx, yy, zzs1)
     ),
    Plot(
         {
          "line width" = 3, 
          color = "red"
         }, 
         Coordinates(linexs, lineys)
        )
    )

fig_name2 = "figures/p2_contour.pdf"
pgfsave(fig_name2, figure2) 



################################################
#                                              #
#        plot characteristic function of       #
#        cheeger set                           #
#                                              #
#################################################

xx = range(0, 1, length=100)
yy = range(0, 1, length=100)
pts = [[x,y] for x in xx, y in yy]
zzs2 = 0.7incheeger.(pts)



# construct outline of cheeger set
# by hand
czs = exp.(im .* range(0, 2π, length=100))

linexs = []
lineys = []
R = cheeger_R(1.0, 1.0)
czs = R * exp.(im .* range(0, 2π, length=100))
push!(linexs, 1.0)
push!(lineys, R)

append!(linexs, 1 .- R .+ real.(czs[1:25]))
append!(lineys, 1 .- R .+ imag.(czs[1:25]))

push!(linexs, R)
push!(lineys, 1.0)

append!(linexs, R .+ real.(czs[26:50]))
append!(lineys, 1 .- R .+ imag.(czs[26:50]))

push!(linexs, 0.0)
push!(lineys, R)

append!(linexs, R .+ real.(czs[51:75]))
append!(lineys, R .+ imag.(czs[51:75]))

push!(linexs, 1-R)
push!(lineys, 0.0)

append!(linexs, 1 .- R .+ real.(czs[76:100]))
append!(lineys, R .+ imag.(czs[76:100]))

linezs = [0.35 for _ in linexs]

# shrink a tiny bit for line to fit in viewport
linexs .= 0.5 .+ 0.985* (linexs .- 0.5)
lineys .= 0.5 .+ 0.985 * (lineys .- 0.5)



figure = @pgf Axis(
    {
      "unit vector ratio" = [1,1], 
      "unit rescale keep size" = false, 
      xmin = minimum(xx), 
      xmax = maximum(xx), 
      ymin = minimum(yy), 
      ymax = maximum(yy), 
      zmin = minimum(zzs2), 
      zmax = maximum(zzs2) 
    }, 
    Plot3(
      { surf, 
        "colormap name"="viridis", 
        shader = "flat", 
#        "faceted color" = "black"
      },
      Table(xx, yy, zzs2)
     )
   )

fig_name = "figures/p1.pdf"
pgfsave(fig_name, figure) 

figure2 = @pgf Axis(
    {
      "axis equal", 
      width="207pt", 
      height="207pt", 
      view = (0,90), 
    }, 
    Plot3(
          { 
           "colormap name"="viridis", 
           shader = "interp",
           "contour filled" = {number = 10, labels = false}
          }, 
      Table(xx, yy, zzs2)
     ), 
    Plot(
         {
          "line width" = 3, 
          color = "red"
         }, 
         Table(linexs, lineys)
        )
   )

fig_name2 = "figures/p1_contour.pdf"
pgfsave(fig_name2, figure2) 

