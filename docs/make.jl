push!(LOAD_PATH,"../src/")
using Documenter
using ExpandLGF

makedocs(
    sitename = "ExpandLGF.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = ["Overview" => "index.md", 
             "Examples" => "examples.md",
             "Docstrings" => "docstrings.md"]
    )
# deploydocs(
#     repo = "github.mit.edu/vanreeslab/ExpandLGF.jl"
# )