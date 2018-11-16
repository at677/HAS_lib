using Test

ϕ=60
a1 = [1,0].*2
a2 = [cos(π*ϕ/180),sin(π*ϕ/180)].*2
l = Generic2DLattice(a1,a2)
conf = CorrugatedMorse(5,0.4,0.1,l)

x = collect(0:0.1:10)
y = collect(0:0.1:10)
v = Array{Float64,2}(undef, length(x),length(y))

for (i1,xi) in enumerate(x), (i2,yi) in enumerate(y)
    z = 0
    v[i1,i2] = potential(conf,xi,yi,z)
end

using Plots
gr()
#plot(x,y,v,seriestype=:heatmap,aspect_ratio=:equal)
#gui()

start = 0
V0 = []
for s in 5:20:100
push!(V0, real(calc_VG([0,0],l,conf,0.2,s)))
end
plot(V0,seriestype=:line)
gui()

#start = 0
#z = 0.3 
#V0 = fill(0.,7,7)
#for gx in -3:3, gy in -3:3
#    V0[gx+4,gy+4] = real(calc_VG([gx,gy],l,conf,z,10))
#end
#@show V0
#plot(V0,seriestype=:heatmap)
#gui()
