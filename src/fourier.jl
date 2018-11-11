export potential,CorrugatedMorse

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

function calc_VG(g,l::Generic2DLattice,conf::CorrugatedMorse)
    rr1 = 0:1.*l.a1
    rr2 = 0:1.*l.a2

    VG = 0
    for r1 in rr1, r2 in rr2
        r = r1.+r2
        VG += potential(conf,r[1],r[2],z)*exp(-1im*g'r)
    end
    VG *= 2π

end


struct CorrugatedMorse
    D
    κ
    h
end

function potential(conf::CorrugatedMorse,x,y,z)
    @unpack D,κ,h = conf

    @. D*(exp(-2z/κ)-2*exp(-z/κ))

end


