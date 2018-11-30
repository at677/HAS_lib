export potential,CorrugatedMorse,calc_VG

using SpecialFunctions

struct CorrugatedMorse{T <: Abstract2DLattice}
    D
    κ
    h
    lattice::T
end

struct HexGrid <: Abstract2DLattice
    a
    l::Generic2DLattice
end

HexGrid(a) = HexGrid(a,Generic2DLattice([a,0],[cos(π/3),sin(π/3)].*a))

function get_w(z,ch::ScatteringChannels,ki,conf::CorrugatedMorse)
    get_V(z,ch,conf) + Diagonal(ch.kz²)
end

function get_V(z,ch::ScatteringChannels,conf::CorrugatedMorse)
    @unpack gi = ch
    l = length(gi)
    V = Array{Float64,2}(undef,l,l)
    for i in 1:l, j in 1:l
        V[i,j] = get_VG(z,gi[i]-gi[j],conf)
    end
    V
end

function get_VG(z,g,conf::CorrugatedMorse{HexGrid})
    @unpack κ,h,D = conf
    if g[1]==0 && g[2]==0
        return D*(exp(-2*κ*z)+2*exp(-κ*z))
    end

    V0 = calc_VG_factor([0,0],conf,1)
    VG = calc_VG_factor(g,conf,V0)

    VG *= exp(-2*κ*z)*D

    return VG
end

function calc_VG_factor(g,conf::CorrugatedMorse{HexGrid},div)
    @unpack κ,h,D = conf
    @unpack a = conf.lattice

    kmax = 10
    α = 2*κ*h*a/3
    VG = 0.0 
    for k in -kmax:kmax
        VG += besseli(k,α)*(besseli(k+g[2],α)*besseli(k+g[1],α) + besseli(k-g[2],α)*besseli(k-g[1],α))
    end
    VG / div
end

function calc_VG(g,l::Generic2DLattice,conf::CorrugatedMorse,z,steps=100)
    VG = 0+0im
    rng = range(0.0,length=steps,1)
    for i in rng , j in rng
        r1 = l.a1.*i
        r2 = l.a2.*j
        r = r1+r2
        VG += potential(conf,r[1],r[2],z)*exp(-1im*(g'*r))
    end
    VG *= 2π/(steps^2)
end


function potential(conf::CorrugatedMorse,x,y,z)
    @unpack D,κ,h,lattice = conf
    @unpack a1,a2 = lattice
    R = [x,y]

    # TODO: only for 60 deg(a1,a2)
    ξ = h*(cos(2π*(a1'*R)/norm(a1)^2)+
           cos(2π*(a2'*R)/norm(a2)^2)+
           cos(2π*((a1-a2)'*R)/norm(a1-a2)^2))
    @. D*(exp(-2*(z-ξ)/κ)-2*exp(-z/κ))

end


