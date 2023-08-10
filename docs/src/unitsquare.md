# First two eigenfunctions on the unit square
We calculate the first two eigenfunctions of 
``\Delta_p`` on the unit square ``[0,1]^2``. 

While this could be done with the algorithm for the 
dynamic ``p``-Laplacian by setting ``T:=id``, the
default constructor for `LMMContext` also fills in 
the right functionals for ``\Delta_p``. 

```julia
using DynamicPLaplacian, LinearAlgebra, Gridap, 
const DPL = DynamicPLaplacian

p = 1.4

# create context
gctx = CartesianGridapContext((0,1,0,1), (100,100))
lmmctx = DPL.LMMContext(p, gctx)
U = gctx.U

# initial guess
initial1 = interpolate(x -> sin( π * x[1]) * sin(π * x[2]), U).free_values
initial2 = interpolate(x -> sin(2π * x[1]) * sin(π * x[2]), U).free_values

result1 = local_min_max([], initial1, lmmctx; 
                       verbose = false, max_it=90, 
                       pgrad_stop=1e-3)

result2 = local_min_max([result1.u], initial2, lmmctx;
                         verbose = false, max_it = 45, 
                         pgrad_stop = 1e-3)
```                         
For plotting we use `GridapMakie`:
```julia
using Gridap.CellData, GridapMakie, GLMakie, FileIO

model = simplexify(gctx.Ω.model)
Ω = Triangulation(model)
U_plot = FESpace(model, ReferenceFE(lagrangian, Float64, 1))

u1 = interpolate(Interpolable(FEFunction(U, result1.u)), U_plot)
u2 = interpolate(Interpolable(FEFunction(U, result2.u)), U_plot)

# plot
fig = Figure(resolution=(800,400))
ax = Axis(fig[1, 1], aspect = 1)
ax = Axis(fig[1, 2], aspect = 1)
pl1 = plot!(fig[1,1], Ω, u1, colorrange=(-1.5,1.5), colormap=:viridis)
pl2 = plot!(fig[1,2], Ω, u2, colorrange=(-1.5,1.5), colormap=:viridis)
Colorbar(fig[1,3], pl1)
```
![input](assets/unit_square.png)

```@eval
using DynamicPLaplacian, Gridap, Gridap.CellData, LinearAlgebra, GridapMakie, GLMakie, FileIO
const DPL = DynamicPLaplacian

p = 1.4

# create context
gctx = CartesianGridapContext((0,1,0,1), (100,100))
lmmctx = DPL.LMMContext(p, gctx)
U = gctx.U

# initial guess
initial1 = interpolate(x -> sin( π * x[1]) * sin(π * x[2]), U).free_values
initial2 = interpolate(x -> sin(2π * x[1]) * sin(π * x[2]), U).free_values

result1 = local_min_max([], initial1, lmmctx; 
                       verbose = false, max_it=90, 
                       pgrad_stop=1e-3)

result2 = local_min_max([result1.u], initial2, lmmctx;
                         verbose = false, max_it = 45, 
                         pgrad_stop = 1e-3)

model = simplexify(gctx.Ω.model)
Ω = Triangulation(model)
U_plot = FESpace(model, ReferenceFE(lagrangian, Float64, 1))

u1 = interpolate(Interpolable(FEFunction(U, result1.u)), U_plot)
u2 = interpolate(Interpolable(FEFunction(U, result2.u)), U_plot)

# plot
fig = Figure(resolution=(800,400))
ax = Axis(fig[1, 1], aspect = 1)
ax = Axis(fig[1, 2], aspect = 1)
pl1 = plot!(fig[1,1], Ω, u1, colorrange=(-1.5,1.5), colormap=:viridis)
pl2 = plot!(fig[1,2], Ω, u2, colorrange=(-1.5,1.5), colormap=:viridis)
Colorbar(fig[1,3], pl1)
save("assets/unit_square.png", fig)
```

