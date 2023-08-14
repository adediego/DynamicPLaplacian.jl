# DynamicPLaplacian.jl

A loose collection of numerical methods reproducing the figures of [[1]](#cit1).

This package implements the local min-max algorithm for computing
eigenpairs described by Yao & Zhou [[2]](#cit2), with adjustments for
treating a p-Laplacian involving dynamics.

## Reproducing figures 
The results from [[1]](#cit1) were obtained using julia v1.9.0.on Fedora 37.
To reproduce them, make sure that a working
LaTeX distribution including the PGFPlots package is installed, 
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
<a id="cit1">[1]</a> Alvaro de Diego, Gary Froyland, Oliver Junge, Peter Koltai, "A dynamic p-Laplacian", 
[arXiv:2308.05947](https://arxiv.org/abs/2308.05947)  
<a id="cit2">[2]</a>
X. Yao & Zhou, Numerical Methods for Computing Nonlinear Eigenpairs: Part I. Iso-Homogeneous Cases,
SIAM J. SCI. COMPUT., Vol. 29, No. 4, pp. 1355–1374
