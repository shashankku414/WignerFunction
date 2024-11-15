module tomo

using JSON
# using Interpolations
using Trapz
using Plots
using Base.Threads
using DelimitedFiles
using KernelDensity

export tomography, dataload, dataload2

function dataload(filename :: String, kde_lims :: Tuple, num_points :: Int64)
    data_array = readdlm(filename)
    phi_data = range(0, 2*pi, length=size(data_array, 2))
    pr = zeros(size(data_array, 2), num_points)
    # for i in 1:size(data_array, 2)
    for i in axes(data_array, 2)
        y = kde(data_array[:, i], boundary=kde_lims, npoints=num_points)
        pr[i, :] = y.density
    end
    return range(kde_lims[1], kde_lims[2], length=num_points), phi_data, pr
end

function dataload2(signal_data, kde_lims :: Tuple, num_points :: Int64)
    # signal_data is 2D array of signal and motor_postion_T
    phi_data = range(0, 2*pi, length=size(signal_data, 2))
    pr = zeros(size(signal_data, 2), num_points)
    # for i in 1:size(data_array, 2)
    for i in axes(signal_data, 2)
        y = kde(signal_data[:, i], boundary=kde_lims, npoints=num_points)
        pr[i, :] = y.density
    end
    return range(kde_lims[1], kde_lims[2], length=num_points), phi_data, pr
end

function K(arg, kc)
    return ((kc^2)/2.)*(1-((kc^2)*(arg^2)/4.)+((kc^4)*(arg^4)/72.)-((kc^6)*(arg^6)/2880.)+((kc^8)*(arg^8)/201600.))
end

function Kor(arg, kc)
    return (cos(kc*arg) + kc*arg*sin(kc*arg) - 1)/(arg^2)
end

function Kcomp(q, p, angle, volt, kc)
    turn = 0.01
    arg = (q*cos(angle)) + (p*sin(angle)) .- volt

    for (i, val) in enumerate(arg)
        if abs(val*kc) < turn
            arg[i] = K(arg[i], kc)
        else
            arg[i] = Kor(arg[i], kc)
        end
    end
    return arg
end

function wigner(args)
    iq, ip, q, p, m, angles, volt, kc = args
    int_s = zeros(Float64, length(angles))
    for (i, angle) in enumerate(angles)
        int_s[i] = trapz(volt, vec(m[i, :] .* Kcomp(q, p, angle, volt, kc)))
    end
    int_phi = trapz(angles, int_s)
    return iq, ip, int_phi/(2*pi^2)
    # return angles, temp
end

function tomography(m, angles, volt; q1=-2, q2=2, p1=-2, p2=2, density=50, kc = 5)
    que = []
    Q = range(q1, q2, length = density)
    P = range(p1, p2, length = density)
    W = zeros(Float64, (density, density))

    for q in 1:density
        for p in 1:density
            push!(que, (q, p, Q[q], P[p], m, angles, volt, kc))
        end
    end
    results = [wigner(var) for var in que]

    for (iq, ip, w) in results
        W[iq, ip] = w
    end
    return collect(Q), collect(P), Matrix(W)
end

end

