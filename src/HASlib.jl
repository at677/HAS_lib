module HASlib
export hexgrid, h2m, cc, fox_goodwin, CCSim

using LinearAlgebra

include("lattice.jl")
include("channels.jl")
include("close_coupling.jl")
include("drift_scan.jl")
include("fourier.jl")
end
