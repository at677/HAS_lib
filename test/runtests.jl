using Test
using HASlib

@test size(HASlib.hexgrid(1)[1])==(7,2)
@test isapprox(h2m(1), 5.56e-69;atol=1e-71,rtol=0)

z = [1,2,3,4,5.0]
w(z)=rand(3,3)*z
@test length(HASlib.fox_goodwin(z,w))==length(w(0))
show(HASlib.cc(3,2,3*1.6e-27,z,4,2,2))
