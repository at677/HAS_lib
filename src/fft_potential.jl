using Revise
using Plots
using Parameters
using LinearAlgebra
using HASlib

function potential(conf::CorrugatedMorse,x,y,z)
    @unpack D,κ,h,lattice = conf
    @unpack a1,a2 = lattice
    a = norm(a1)
    R = [x,y]
    # TODO: only for 60 deg(a1,a2)
    ξ =  h*(cos(2π/a*(x-y/sqrt(3)))
            +cos(2π/a*(x+y/sqrt(3)))
            +cos(2π/a*(2y/sqrt(3))))
    D*(exp(-2*(z-ξ)/κ)-2*exp(-z/κ))
end

function coupling(a,h,xi,g,divide = 1 ;kmax = 10)
    α = 2*xi*h*a/3
    VG = 0.0 
    for k in -kmax:kmax
        VG += besseli(k,α)*(besseli(k+g[2],α)*besseli(k-g[1],α) + besseli(k-g[2],α)*besseli(k+g[1],α))
    end
    VG / divide
end

function main()
@info "start"

a = 1

x = range(0,a*(1+cosd(60)),length=100)
y = range(0,a*sind(60),length=100)
z = zeros(length(x),length(y))

hg = HASlib.HexGrid(a)
pot = HASlib.CorrugatedMorse(1,10,1,hg.l)

for (i,xi) in enumerate(x), (j,yj) in enumerate(y)
    z[i,j] = potential(pot,x[i],y[j],-1)
end

plot(y,x,z,seriestype=:heatmap,aspect_ratio=:equal,size=(300,300))

g,gi = HASlib.hexgrid(5,hg.l.a1,hg.l.a2)
scatter!(g[:,2],g[:,1],legend=:none)

data = zeros(11*11,3)
for i in 0:10, j in 0:10
    a = - i*hg.l.a1/11 - j*hg.l.a2/11
    row = (i+1)+j*11
    data[row,1] = a[1]
    data[row,2] = a[2]
    data[row,3] = potential(pot,a[1],a[2],-1)
end

scatter!(data[:,2],data[:,1],marker=2*(-minimum(data[:,3]).+data[:,3])) |> display

k = zeros(ComplexF64,size(g,1))
for i in 1:size(g,1)
    for row in data
        r = hcat(data[:,1], data[:,2])
        k[i] = sum(data[:,3].*exp.(-1im*g[i,:]'*r'))
    end
end

scatter(g[:,1],g[:,2],marker=(log.(abs.(k)).-9)*10,aspect_ratio=:equal)

x = range(-a,a,length=100)
y = range(-a,a,length=100)
z = zeros(length(x),length(y))

for (i,xi) in enumerate(x), (j,yj) in enumerate(y)
    z[i,j] =  begin
        r = [xi,yj]
        @show 1im*g*r
        real.(sum(k .* exp.(1im*g*r)))
    end
end

plot(y,x,z,seriestype=:heatmap,aspect_ratio=:equal,size=(300,300))


end
main()