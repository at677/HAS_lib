using SpecialFunctions, LinearAlgebra

struct CCSim
    n_total
    n_open
    n_closed
    kz2
    g
    gi
end

function Base.show(result::CCSim)
    print("total channels: ",result.n_total)
    print(" (open: ",result.n_open)
    print(", closed: ",result.n_closed,")")
    print(", kz2: ",result.kz2,")")
end

function h2m(m)
    hb = 1.054571800e-34
    hb^2/(2*m)
end

function fox_goodwin_step!(w_p1,w_0,w_m1,r)
        #aa = 2*I+10*w_0
        #bb = (I-w_m1)*r
        #cc = I-w_p1
        #r .= (aa.-bb)\cc
        r .= (2*I+10*w_0 - (I-w_m1)*r )\(I-w_p1)
        w_m1 .= w_0
        w_0 .= w_p1
        r
end

function fox_goodwin(z,w!,n)
    h = z[2]-z[1]
    c1 = h^2/12

    w_m1 = Array{typeof(z[1]),2}(undef,n,n)
    w_0 = similar(w_m1)
    w_p1 = similar(w_m1)
    r  = fill!(similar(w_0),0)

    w!(w_m1,z[1])
    lmul!(c1,w_m1)
    w!(w_0,z[2])
    lmul!(c1,w_0)

    for i in 3:length(z)
        w!(w_p1,z[i])
        lmul!(c1,w_p1)
        fox_goodwin_step!(w_p1,w_0,w_m1,r)
    end
    r
end

function coupling(a,h,xi,g,divide = 1 ;kmax = 10)
    α = 2*xi*h*a/3
    VG = 0.0 
    for k in -kmax:kmax
        VG += besseli(k,α)*(besseli(k+g[2],α)*besseli(k-g[1],α) + besseli(k-g[2],α)*besseli(k+g[1],α))
    end
    VG / divide
end

function get_w_fun(channels,a,h,D,xi,m)
	n = length(channels.kz2)
	V0(z) = D*(exp(-2*xi*z) - 2*exp(-xi*z))
	V1(z) = D*exp(-2*xi*z)
        M_E = Array{Float64,2}(undef,n,n)
	c0 = coupling(a,h,xi,(0,0))
	for i in 1:size(M_E,1), j in 1:size(M_E,2)
	    g1,g2 = channels.gi[i].-channels.gi[j]
	    M_E[j,i] = coupling(a,h,xi,(g1,g2),c0)
	end
	for i in 1:size(M_E,1)
	    M_E[i,i] = 0
	end
        function w!(w,z)
            v0 = V0(z)/h2m(m)
            v1 = V1(z)/h2m(m)
            for i in 1:size(w,1), j in 1:size(w,2)
                w[i,j] = v1 * M_E[i,j]
            end
            for i in 1:size(w,1)
                w[i,i] = v0 - channels.kz2[i]
            end
        end
end

function cc(channels::ScatteringChannels,a,h,D,xi, theta,m,z;st_p_wl=100)
    #find lowest open channel - split open close channels
    last_open = findlast(x -> x > 0, channels.kz2)

    open_index = 1:last_open
    close_index = last_open+1:length(channels.kz2)

    open_kz2 = @view channels.kz2[open_index]
    open_gi = @view channels.gi[open_index,:]

    # calculate z steps
    step = 2*π / (st_p_wl*sqrt(maximum(channels.kz2)))
	zrange = collect(z[1]:step:z[2])
    
	# get W matrix function
    w = get_w_fun(channels,a,h,D,xi,m)
	
	# solve trasition matrix r
    r = fox_goodwin(zrange,w,length(channels.kz2))

    S0_m1 = zeros(size(r,1))
    S0_00 = zeros(size(r,1))
    CE_m1 = zeros(size(r,1))
    CE_00 = zeros(size(r,1))

    for (i,erg) in enumerate(channels.kz2)
        if erg > 0 # open
            kz = sqrt(erg)
            S0_m1[i] = sin(kz*zrange[end-1])/sqrt(kz)
            S0_00[i] = sin(kz*zrange[end])/sqrt(kz)
            CE_m1[i] = cos(kz*zrange[end-1])/sqrt(kz)
            CE_00[i] = cos(kz*zrange[end])/sqrt(kz)
        else # closed
            kz = sqrt(abs(erg))
            CE_m1[i] = exp(-kz*zrange[end-1])/sqrt(kz)
            CE_00[i] = exp(-kz*zrange[end])/sqrt(kz)
        end
    end

    S0_m1 = Diagonal(S0_m1)
    S0_00 = Diagonal(S0_00)
    CE_m1 = Diagonal(CE_m1)
    CE_00 = Diagonal(CE_00)

    SS = S0_m1 .- r*S0_00
    CC = r*CE_00 .- CE_m1

    KK = CC\SS

    K_open = KK[open_index,open_index]
    S = broadcast(abs2,(I+K_open*(0+1im))\(I-K_open*(0+1im)))

    spec = findfirst(x -> x==[0,0], channels.g[open_index])

    result = hcat(collect(S[spec,:]),open_gi,open_kz2)
end

