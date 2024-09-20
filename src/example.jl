using DelimitedFiles
include("./tomo.jl")
using .tomo


load_dir = "../example_data" # Specify the location of the experimetal data
time_delay = readdlm("$(load_dir)/timedelay.dat") #import the time delay array from same
work_dir = "../work" # Specify the location where do you want to save the computed data
println("Number of available threads = $(Threads.nthreads())") # For parallel computing > 1
Threads.@threads for i in eachindex(time_delay)
    save_dir = "$(work_dir)/$(i-1)"
    # run(`mkdir -p $(save_dir)`) # Make a different directory for each time delay
    mkpath("$(save_dir)")
    x_range, phi_range, pr_data = dataload("$(load_dir)/$(i-1)", (-1000, 1000), 1000) #Calculate the PDF from data
    writedlm("$(save_dir)/i-.dat", collect(x_range))
    writedlm("$(save_dir)/phi.dat", collect(phi_range))
    writedlm("$(save_dir)/pr.dat", pr_data) #Save the PDF data to plot histograms
    # Calculate the field quadratures and the Wigner function values
    # q1,q2 are the lower and higher limits of Quadrature Q
    # p1,p2 are the lower and higher limits of Quadrature P
    # density will specify the number of grid points
    # kc is the cut off limit for filter Kernal (Do not change)
    Q, P, W = tomography(pr_data, collect(phi_range), collect(x_range); q1=-5, q2=5, p1=-5, p2=5, density=100, kc=5)
    writedlm("$(save_dir)/Q.dat", Q)
    writedlm("$(save_dir)/P.dat", P)
    writedlm("$(save_dir)/W.dat", W) #Save the values for Wigner plots
end