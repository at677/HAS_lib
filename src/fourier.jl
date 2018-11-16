export potential,CorrugatedMorse,calc_VG

struct WCache
    M_no_z
end

function get_w(z,ch::ScatteringChannels,ki)
    norm²_g = @. norm(ch.g)^2
    get_V(z,ch) + I*norm(ki)^2 + Diagonal(norm²_g) 
end

function get_V(z,ch::ScatteringChannels)
    l = length(ch.g)
    V = Array{Float64,2}(undef,l,l)
    for i in 1:l, j in 1:l
        V[i.j] = get_VG(z,gi[i]-gi[j])
    end
end

function get_VG(gi,z)
end

function build_VG_Cache(max_gi,z)

end

struct CorrugatedMorse
    D
    κ
    h
    lattice::Generic2DLattice
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


