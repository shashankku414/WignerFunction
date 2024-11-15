using DelimitedFiles
include("./tomo.jl")
using .tomo
using Printf
using KernelDensity


load_dir = "/Volumes/Seagate EXP/quant_optics_exp_data/2024_10_05/spectra" # Specify the location of the experimetal data
main_dir = "/Volumes/One Touch/Russ_Data/2024_10_05/Experiment/MgO_sig_with_ref_full_spectra/Spectras"
work_dir = "../work_2024_10_05_w_761nm"
motor_tau_range = -390e-6:30e-6:390e-6 #motor position in mm for time delay between probe and gate
motor_T_range = -200e-6:4e-6:200e-6 #motor position in mm for time delay between probe and gate

# motor_tau_range = -390e-6 #motor position in mm for time delay between probe and gate
# motor_T_range = -168e-6 #motor position in mm for time delay between probe and gate

motor_tau = [@sprintf("%.6f", v) for v in motor_tau_range]
motor_T = [v == 0.0 ? "-0.000000" : @sprintf("%.6f", v) for v in motor_T_range]


# 2*motor position = Path length 
# 1nm = 1e-6 mm motor postion
# 1 nm Path length = 10/3 as

tau_range = 10/3 * 1e6 * 2 .* motor_tau_range
T_range = 10/3 * 1e6 * 2 .* motor_T_range
scan_range = 1

num_frames = 1000 #This is the number of frames that is measured for each T, tau and scan
freq_index = 70 # This is calculated from finding the maxium intensity of pixel

println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1
Threads.@threads for i in eachindex(tau_range)
    save_dir = "$(work_dir)/tau_$(tau_range[i])as/scan_$(scan_range)"
    mkpath("$(save_dir)")
    data_top = zeros((num_frames, length(T_range)))
    data_bottom = zeros((num_frames, length(T_range)))
    for k in eachindex(T_range)
        # println(motor_T[k])
        filename_top = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_top_spectra.txt"
        filename_bottom = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_bottom_spectra.txt"
        data_top[:, k] = readdlm(filename_top)[:, freq_index]./10000 #To normalize the values (Can be number of pixels as well)
        data_bottom[:, k] = readdlm(filename_bottom)[:, freq_index]./10000 #To normalize the values (Can be number of pixels as well)
    end
    _, _, pr_top_data = dataload2(data_top, (-4, 4), 1000) #Calculate the PDF from data
    _, _, pr_bottom_data = dataload2(data_bottom, (-4, 4), 1000) #Calculate the PDF from data
    x_range, phi_range, pr_data = dataload2(data_top - data_bottom, (-4, 4), 1000) #Calculate the PDF from data
    writedlm("$(save_dir)/pr_top.dat", pr_top_data) #Save the PDF data to plot histograms
    writedlm("$(save_dir)/pr_bottom.dat", pr_bottom_data) #Save the PDF data to plot histograms
    writedlm("$(save_dir)/i-.dat", collect(x_range))
    writedlm("$(save_dir)/phi.dat", collect(phi_range))
    writedlm("$(save_dir)/pr-.dat", pr_data) #Save the PDF data to plot histograms
    Q, P, W = tomography(pr_data, collect(phi_range), collect(x_range); q1=-5, q2=5, p1=-5, p2=5, density=100, kc=5)
    writedlm("$(save_dir)/Q.dat", Q)
    writedlm("$(save_dir)/P.dat", P)
    writedlm("$(save_dir)/W.dat", W) #Save the values for Wigner plots
end

writedlm("$(work_dir)/timedelay.dat", collect(tau_range))


###################### copy large files ###########################


# main_dir = "/Volumes/One Touch/Russ_Data/2024_10_05/Experiment/MgO_sig_with_ref_full_spectra/Spectras"
# load_dir = "/Volumes/Seagate EXP/quant_optics_exp_data/2024_10_05/spectra"

# println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1
# for i in eachindex(tau_range)
#     Threads.@threads for j in eachindex(scan_range)
#         for k in eachindex(T_range)
#             for position in ["top", "bottom"]
#                 filename = "T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_$(position)_spectra.txt"
#                 cp("$(main_dir)/$(filename)", "$(load_dir)/$(filename)")
#             end
#             println("T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range)) done")
#         end
#     end
# end


# println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1
# Threads.@threads for i in eachindex(tau_range)
#     data_top = zeros((num_frames, length(T_range)))
#     data_bottom = zeros((num_frames, length(T_range)))
#     for k in eachindex(T_range)
#         filename_top = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_top_spectra.txt"
#         filename_bottom = "$(load_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_bottom_spectra.txt"
#         try
#             data_top[:, k] = readdlm(filename_top)[:, freq_index]./10000 #To normalize the values (Can be number of pixels as well)
#         catch e
#             println("There is something wrong with T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_top_spectra.txt")
#             cp("$(main_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_top_spectra.txt", filename_top; force=true)
#             continue
#         end
#         try
#             data_bottom[:, k] = readdlm(filename_bottom)[:, freq_index]./10000 #To normalize the values (Can be number of pixels as well)
#         catch e
#             println("There is something wrong with T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_bottom_spectra.txt")
#             cp("$(main_dir)/T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range))_bottom_spectra.txt", filename_bottom; force=true)
#             continue
#         end
#         println("T_$(motor_T[k])_tau_$(motor_tau[i])_scan_$(string(scan_range)) done")
#     end
# end