using Test
using HASlib

const mel = 1.6021766208000002e-22
const A = 1e-10
const m_He3 = 5.0082343984338715e-27

include("test_python.jl")
include("test_cc.jl")
include("test_driftscan.jl")

#include("fourier.jl")

@test size(HASlib.hexgrid(1)[1])==(7,2)
@test isapprox(h2m(1), 5.56e-69;atol=1e-71,rtol=0)

l = Generic2DLattice([-1,0.5].*10,[0,1])
r = calc_radius_factor(l)
ch =  get_channels(l.b1,l.b2,3,π/2,Int(ceil(3*sin(π/2)/r)))
@test 10 == length(sort_and_cut(ch,5,5).g)




#z = [1,2,3,4,5.0]
#w(z)=rand(3,3)*z
#@test length(HASlib.fox_goodwin(z,w))==length(w(0))
#show(HASlib.cc(3,2,3*1.6e-27,z,4,2,2))
