using DelimitedFiles
include("./tomo.jl")
using .tomo
using Printf
using KernelDensity


load_dir =  "/Users/shashank/Downloads/tau_-390_files" # Specify the location of the experimetal data
# motor_tau_range = -390e-6:30e-6:390e-6 #motor position in mm for time delay between probe and gate
motor_tau_range = -390e-6 # The first tau and first T is empty
motor_T_range = -200e-6:8e-6:200e-6 #motor position in mm for time delay between probe and gate

# motor_tau_range = -390e-6 #motor position in mm for time delay between probe and gate
# motor_T_range = -168e-6 #motor position in mm for time delay between probe and gate

motor_tau = [@sprintf("%.6f", v) for v in motor_tau_range]
# motor_T = [v == 0.0 ? "-0.000000" : @sprintf("%.6f", v) for v in motor_T_range]
motor_T = [@sprintf("%.6f", v) for v in motor_T_range]

# 2*motor position = Path length 
# 1nm = 1e-6 mm motor postion
# 1 nm Path length = 10/3 as

tau_range = 10/3 * 1e6 * 2 .* motor_tau_range
T_range = 10/3 * 1e6 * 2 .* motor_T_range
scan_range = 0


num_frames = 1000 #This is the number of frames that is measured for each T, tau and scan
# 679 is the freq_index at which the pixel intesity are maximum i.e. wavelength=808nm
# Each pixel have a grid spaceing of 0.064nm
# itask = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
# freq_index = 679-(itask-1)*31 # 31 index spaceing corresponds to 31*0.064 ~ 2nm

work_dir = "../test_work/"

println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1

# test_data = readdlm("$(load_dir)/T_$(motor_T[2])_tau_$(motor_tau[1]).txt")
# println(size(sum(test_data, dims=2)))

Threads.@threads for i in eachindex(tau_range)
    if motor_tau_range[i] == -120e-6 # Neglect the bad data points
        continue
    end
    for j in eachindex(scan_range)
        save_dir = "$(work_dir)/tau_$(tau_range[i])as/scan_$(scan_range[j])"
        mkpath("$(save_dir)")
        data_rf = zeros((num_frames, length(T_range)))
        for k in eachindex(T_range)
            # println(motor_T[k])
            filename = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i]).txt"
            data_rf[:, k] = sum(readdlm(filename), dims=2) #To normalize the values (Can be number of pixels as well)
        end
        x_range, phi_range, pr_data = dataload2(data_rf, (-400, 400), 100000) #Calculate the PDF from data
        writedlm("$(save_dir)/i-.dat", collect(x_range))
        writedlm("$(save_dir)/phi.dat", collect(phi_range))
        writedlm("$(save_dir)/pr-.dat", pr_data) #Save the PDF data to plot histograms
        Q, P, W = tomography(pr_data, collect(phi_range), collect(x_range); q1=-500, q2=500, p1=-500, p2=500, density=100, kc=5)
        writedlm("$(save_dir)/Q.dat", Q)
        writedlm("$(save_dir)/P.dat", P)
        writedlm("$(save_dir)/W.dat", W) #Save the values for Wigner plots
    end
end

writedlm("$(work_dir)/timedelay.dat", filter(x->x!=-800, collect(tau_range)))