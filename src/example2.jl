using DelimitedFiles
include("./tomo.jl")
using .tomo
using Printf
using KernelDensity


load_dir = "../2024_09_20/Experiment/MgO_sig_with_ref/i_Values" # Specify the location of the experimetal data
work_dir = "../work_2024_09_20"
motor_tau_range = -390e-6:30e-6:390e-6 #motor position in mm for time delay between probe and gate
motor_T_range = -200e-6:4e-6:200e-6 #motor position in mm for time delay between probe and gate

motor_tau = [@sprintf("%.6f", v) for v in motor_tau_range]
motor_T = [v == 0.0 ? "-0.000000" : @sprintf("%.6f", v) for v in motor_T_range]


# 2*motor position = Path length 
# 1nm = 1e-6 mm motor postion
# 1 nm Path length = 10/3 as
tau_range = 10/3 * 1e6 * 2 .* motor_tau_range
T_range = 10/3 * 1e6 * 2 .* motor_T_range
scan_range = 1:1:9

num_frames = 100 #This is the number of frames that is measured for each T, tau and scan

println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1
for i in eachindex(tau_range)
    Threads.@threads for j in eachindex(scan_range)
        save_dir = "$(work_dir)/tau_$(tau_range[i])as/scan_$(scan_range[j])"
        mkpath("$(save_dir)")
        data_array = zeros((num_frames, length(T_range)))
        for k in eachindex(T_range)
            filename = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range[j])).txt"
            data_array[:, k] = readdlm(filename)./10000 #To normalize the values (Can be number of pixels as well)
        end
        x_range, phi_range, pr_data = dataload2(data_array, (-2, 2), 1000) #Calculate the PDF from data
        writedlm("$(save_dir)/i-.dat", collect(x_range))
        writedlm("$(save_dir)/phi.dat", collect(phi_range))
        writedlm("$(save_dir)/pr.dat", pr_data) #Save the PDF data to plot histograms
        Q, P, W = tomography(pr_data, collect(phi_range), collect(x_range); q1=-5, q2=5, p1=-5, p2=5, density=100, kc=5)
        writedlm("$(save_dir)/Q.dat", Q)
        writedlm("$(save_dir)/P.dat", P)
        writedlm("$(save_dir)/W.dat", W) #Save the values for Wigner plots
    end
end

writedlm("$(work_dir)/timedelay.dat", collect(tau_range))