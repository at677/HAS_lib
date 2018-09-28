using Test
using HASlib

@test size(HASlib.hexgrid(1)[1])==(7,2)
@test isapprox(h2m(1), 5.56e-69;atol=1e-71,rtol=0)
println("cc:", HASlib.cc(3,2,3*1.6e-27,4))
