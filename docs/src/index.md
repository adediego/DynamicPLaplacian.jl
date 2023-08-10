# DynamicPLaplacian.jl

Numerical methods for computing eigenfunctions
of the dynamic p-Laplacian and code reproducing
the results of [[1]](@ref References).

This package implements the local min-max algorithm for computing
eigenpairs described by Yao & Zhou [[2]](@ref References), with adjustments for
treating a p-Laplacian involving dynamics.

## Installation
Run
```bash
julia -e 'using Pkg; Pkg.add("git@github.com:adediego/DynamicPLaplacian.jl.git")'
``` 
or
```
]add git@github.com:adediego/DynamicPLaplacian.jl.git
```
from the julia REPL. 

## Reproducing
The results from [[1]](@ref References) were obtained using julia v1.8.3.
To reproduce them, make sure that a working
``\LaTeX`` distribution including the PGFPlots package is installed, 
as the plots are generated using the julia wrapper [PGFPlotsX](https://kristofferc.github.io/PGFPlotsX.jl/stable/).
Then run:
```bash
cd examples/reprocude_plots
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project run_all.jl
```
This will take a few hours and will save the plots in
`examples/reproduce_plots/figures/`. 

If there is a problem during plotting, the experiments do not 
have to be run again, as the results are cached on disc in `./results`.
Run `julia --project run_plots.jl` to just run the plotting code.

## References
[1] Alvaro de Diego, Gary Froyland, Oliver Junge, Peter Koltai, "A dynamic p-Laplacian" <arxiv preprint>\ 
[2] X. Yao & Zhou, Numerical Methods for Computing Nonlinear Eigenpairs: Part I. Iso-Homogeneous Cases,
SIAM J. SCI. COMPUT., Vol. 29, No. 4, pp. 1355–1374

