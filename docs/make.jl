push!(LOAD_PATH,"../src/")
using Documenter, DynamicPLaplacian


pages = ["Home"    => "index.md", 
         "Examples"=> Any[
            "Unit square" => "unitsquare.md",
            "Unit disk"   => "unitdisk.md",
            "Transitory double gyre" => "rotdoublegyre.md"
           ]
        ]
makedocs(sitename="DynamicPLaplacian.jl"; pages)

deploydocs(
    repo = "github.com/adediego/DynamicPLaplacian.jl"
   )
