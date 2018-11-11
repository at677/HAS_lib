using Documenter
using HASlib

makedocs(
    sitename = "HASlib",
    format = :html,
    modules = [HASlib]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
