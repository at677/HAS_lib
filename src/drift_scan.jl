# change incident energy at 0 order peak

using Parameters

export DriftScan, driftscan

struct DriftScanResult
    energies
    zeroOrderSignal
end

@with_kw struct DriftScan
    energies
    a
    h
    D
    xi
    theta
    m
    z
    ch_open = -1
    ch_closed = 5
end



function driftscan(param::DriftScan)
    @unpack ch_open, ch_closed, energies, theta, a, h, D, xi, m, z = param
    signal = []
    lattice = HexGrid(a)
    b1, b2 = rotate_reciprocal(lattice.l, [1, 0])
    for energy in energies
        ki = sqrt(energy / h2m(m))
        channels = get_channels(b1, b2, ki, theta)
        sortedchannels = HASlib.sort_and_cut(channels, ch_open, ch_closed)
        res = cc(sortedchannels, a, h, D, xi, theta, m, z)

        found = false
        for i = 1:size(res, 1)
            if res[i, 2] == [0, 0]
                push!(signal, res[i, 1])
                found = true
            end
        end
        if found == false
            @warn "No channel [0,0] found"
        end

    end
    DriftScanResult(energies, signal)
end
