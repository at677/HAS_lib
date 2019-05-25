using Documenter
using HASlib

makedocs(
    sitename = "HASlib",
    modules = [HASlib],
    format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/feanor12/HASlib.jl.git",
)

