# module JuliaPhotochem

using DifferentialEquations
using DiffEqBase
using Catalyst
using Random
using Surrogates
using Plots, GraphRecipes

export ChemSys, random_sims, JPM
ENV["PATH"] = "/opt/homebrew/bin:" * ENV["PATH"] # to get the graphs.jl access
press=1.0
tempk=298.0

# seed the random number generator
Random.seed!(42)

function newmodel(press=1.0::Float64, tempk=298.0::Float64)
    # function definition to convert from cm3 molec-1 s-1 to (ppm min)-1
	# 7.34e15 is Avogadro's number divided by R in L atm K-1 mol-1 divided by 1000 to convert to cm3 from L and 1e6 to convert to (ppm)-1
	# These units are valid for bimolecular reactions
    # rkppmmin(rk,press,tempk) = rk * 60.0 *  7.34e+15 * press / tempk 

    # make a function arrhenius(A,ER,tempk) to calculate the rate constant
    # A is the pre-exponential factor, ER is the activation energy in J/mol, and tempk is the temperature in K
    arrhenius(A,EdivR,tempk) = A * exp(-EdivR/tempk)
    molcm3frompressure(tempk,press) = 6.022e23 * press / (.0821 * tempk) / 1e3 # 6.022e23 is Avogadro's number, .0821 is R in L atm K-1 mol-1, and there are 1e3 cm3 in a L
    function termolecular(k₀298,n, k∞298, m, tempk,press)
        # calculate rate constant in the low pressure limit
        k₀ = k₀298 * (tempk/298.0)^(-n) # cm6 molec-2 s-1
        # calculate rate constant in the high pressure limit
        k∞ = k∞298 * (tempk/298.0)^(-m) # cm3 molec-1 s-1
        # calculate M, the total gas Concentration
        M = molcm3frompressure(tempk,press) # molec/cm3
        # calculate the effective second-order rate constant
        kf = ( k∞ * k₀ * M ) / (k∞ + k₀ * M ) * 0.6^(1.0 / ( 1.0 + (log10(k₀*M/k∞))^2.0 ))
        return kf 
    end

    # Reaction rate constants

    # Reaction 2 is O + O2 → O3
    # updated June 2024, JPL 19-5 page 434 Table 2-1
    # original was rkppmmin(6.0e-34 * ( tempk/300. )^(-2.3) * 7.34e+21 * press / tempk, press, tempk)
    rk2(tempk,press) = molcm3frompressure(tempk,press)*6.1e-34*(tempk/298)^(-2.4)

    # Reaction 3 is O3 + NO → NO2 + O2
    # updated June 2024, JPL 19-5 page 84 Table 1C 
    # original was rkppmmin(2.e-12 * exp(-1400/tempk), press, tempk)
    rk3(tempk) = arrhenius(3.0e-12,1500,tempk) 

    # Reaction 6 is HCHO + HO + O2 → HO2 + CO + H2O
    # Kept original rate rkppmmin(1.2e-14 * tempk^1 * exp(+287/tempk), press, tempk)
    rk6(tempk) = tempk*arrhenius(1.2e-14,-287,tempk)

    # Reaction 7 is HO2 + NO → HO + NO2
    # updated June 2024, JPL 19-5 page 83 Table 1C
    # original was rk7(press, tempk) = rkppmmin(3.7e-12 * exp(+250/tempk), press, tempk)
    rk7(tempk) = arrhenius(3.44e-12,-260,tempk)

    # Reaction 8 is HO + NO2 → HNO3
    # updated June 2024, JPL 19-5 page 434 Table 2-1
    # the original rate is below
    # function rk8(press, tempk)
    #     rk0 = 2.0e-30 * ( tempk/300 )^(-3.0)
    #     rki = 2.5e-11 * ( tempk/300 )^0
    #     rm = 7.34e+21 * press / tempk
    #     rb = 0.6^(1.0 / ( 1.0 + (log10(rk0*rm/rki))^2.0 ))
    #     rk8 = rk0 * rm * rb / (1.0 + rk0*rm/rki)
    #     rk8 = rkppmmin(rk8,press,tempk)
    # end
    rk8(tempk,press) = termolecular(1.8e-30,3.0,2.8e-11,0,tempk,press)

    # Reaction 10 is HO2H + HO → H2O + HO2
    # updated June 2024, JPL 19-5 page 73
    # "the recommendation is a temperature independent value of 1.8e-12 cm3 molec-1 s-1"
    # original was rkppmmin(3.3e-12 * exp(-200.0/tempk),press,tempk)
    rk10() = 1.8e-12

    # hv(t) = abs(ceil(0.1*sin(t/720*pi - pi/2)))
    # hv_1(t) = max(sin((t+hv_shift)/720*pi - pi/2),0.0)
    # hv_2(t) = max((sin((t+hv_shift)/720*pi - pi/2) + 1.0) / 2.0,0.0)
    # hv(t) = hv_1(t) * 0.9 + hv_2(t) * 0.1

# Set the value of hv
hv = 1/60 #0 # 60.0*1e40

# Evaluate the rate constant functions to get their values
p = [hv, rk2(press,tempk), rk3(tempk), rk6(tempk), 
    rk7(tempk), rk8(press,tempk), rk10()]


rn = @reaction_network begin
    0.5 * hv, NO2 → NO + O3 + O
    # rk2, O + O2 → O3  # keep oxygen in this rate law
    rk3, O3 + NO → NO2 + O2
    0.015  * hv / O2^2, HCHO + 2O2 → 2HO2 + CO  # rate law not proportional to O2
    # 0.015  * hv, HCHO → 2HO2 + CO 
    0.022  * hv, HCHO → H2 + CO
    rk6 / O2, HCHO + HO + O2 → HO2 + CO + H2O # rate law not proportional to O2
    rk7, HO2 + NO → HO + NO2
    rk8, HO + NO2 → HNO3
    0.0003  * hv, HO2H → 2HO
    rk10, HO2H + HO → H2O + HO2
end hv rk2 rk3 rk6 rk7 rk8 rk10

    nspecs = size(species(rn),1)

    specmap = Dict()
    for (k, v) in speciesmap(rn)
        specmap[replace(string(k), r"\(t\)$"=>"")] = v # remove '(t)' from names
    end
    natoms = 4 # C N H O
    atoms = zeros(Float64, nspecs, natoms)
    atoms[specmap["O3"], :] = [0 0 0 3]
    atoms[specmap["NO"], :] = [0 1 0 1]
    atoms[specmap["NO2"], :] = [0 1 0 2]
    atoms[specmap["HCHO"], :] = [1 0 2 1]
    atoms[specmap["HO2"], :] = [0 0 1 2]
    atoms[specmap["HO2H"], :] = [0 0 2 2]
    atoms[specmap["HO"], :] = [0 0 1 1]
    atoms[specmap["O"], :] = [0 0 0 1]
    atoms[specmap["HNO3"], :] = [0 1 1 3]
    atoms[specmap["CO"], :] = [1 0 0 1]
    atoms[specmap["H2"], :] = [0 0 2 0]
    atoms[specmap["H2O"], :] = [0 0 2 1]
    atoms[specmap["O2"], :] = [0 0 0 2]

    return rn, atoms, specmap, p
end

struct ChemSys
    rn
    atoms::Array{Float64, 2}
    specmap::Dict
    p::Array{Float64, 1}
end


function ppb_to_molec_cm3(ppb::Float64, tempk::Float64, press::Float64)
    molec_cm3 = ppb / 1e9 * press / .0821 / tempk * 6.022e23 / 1e3
    return molec_cm3
end
function molec_cm3_to_ppb(molec_cm3::Float64, tempk::Float64, press::Float64)
    ppb = molec_cm3 / 6.022e23 * 1e3 * .0821 * tempk / press * 1e9
    return ppb
end



ChemSys(press=1.0::Real, tempk=298.0::Real) = ChemSys(newmodel(press, tempk)...)
jpm = ChemSys()
graph = Graph(jpm.rn)


# initialize c0
c0 = zeros(Float64,length(jpm.specmap))
# c0[jpm.specmap["O3"], :] .= .02 + 0.001*10^(2*rand())   # 03 range: 0.021 - 0.12 ppm
# c0[jpm.specmap["NO"], :] .= 0.0015*10^(2*rand())        # NO range: 0.0015 - 0.15
# c0[jpm.specmap["NO2"], :] .= 0.0015*10^(2*rand())       # NO2 range: 0.0015 - 0.15
# c0[jpm.specmap["HCHO"], :] .= 0.02*10^(2*rand())        # HCHO range: 0.02 - 2 ppm
# c0[jpm.specmap["HO2"], :] .= 1.0e-05 * rand()           # HO2. range: 1 - 10 *ppt*
# c0[jpm.specmap["HO2H"], :] .= 0.01*rand()	            # HO2H range: 0.001 - 0.01 ppm
c0[jpm.specmap["O2"], :] .= 0.21e6	                    # O2: 0.21e6 ppm

# Initialize concentrations from uniform distribution, ranges per SI of https://doi.org/10.1029/2020JD032759
c0[jpm.specmap["O3"], :]   .= 100*rand()*1e-3            
c0[jpm.specmap["NO"], :]   .= 10*rand()*1e-3     
c0[jpm.specmap["NO2"], :]  .= 10*rand()*1e-3
c0[jpm.specmap["HCHO"], :] .= 20*rand()*1e-3    
c0[jpm.specmap["HO2"], :]  .= 0.5*rand()*1e-3
c0[jpm.specmap["HO2H"], :] .= 10*rand()*1e-3     
c0 = ppb_to_molec_cm3.(c0*1e3, tempk, press)


tspan = (0, 900)
prob = ODEProblem(jpm.rn, c0, tspan, jpm.p)

# solve the problem using Rosenbrock23
sol = solve(prob, Rosenbrock23(), saveat=1.0)



# plot O3, NO, NO2, HCHO trajectories
# Define the variables to plot and get their indices
vars_to_plot = ["O3", "NO", "NO2", "HCHO"]
indices = [jpm.specmap[var] for var in vars_to_plot]
# Extract the variables at every time step
selected_vars = [state[indices] for state in sol.u]
# Transpose the matrix to get each variable as a separate vector
selected_vars_transposed = transpose(hcat(selected_vars...))
# convert to ppb
selected_vars_transposed = molec_cm3_to_ppb.(selected_vars_transposed, tempk, press)
# Plot each variable as its own subplot
# Create a new plot with a subplot for each variable
p = plot(layout = (length(vars_to_plot), 1))

# Add each variable to the plot
for (i, var) in enumerate(vars_to_plot)
    plot!(p[i], sol.t, selected_vars_transposed[:,i],legend=false)
    ylabel!(p[i], var * " [ppb]")

end
display(p)


# run the specified number of reference model simulations for the given timelength,
# and return results spaced at the given time step length,
# both in minutes.
# function random_sims(m::ChemSys, nruns::Int, minutes, saveat, idx)
#     c0 = zeros(Float64,length(m.specmap))
#         c0[m.specmap["O3"], :] .= 0.001*10^(2*rand())         # 03 range: 0.001 - 0.1 ppm
#         c0[m.specmap["NO"], :] .= 0.0015*10^(2*rand())        # NO range: 0.0015 - 0.15
#         c0[m.specmap["NO2"], :] .= 0.0015*10^(2*rand())        # NO2 range: 0.0015 - 0.15
#         c0[m.specmap["HCHO"], :] .= 0.02*10^(2*rand())          # HCHO range: 0.02 - 2 ppm
#         c0[m.specmap["HO2"], :] .= 1.0e-05 * rand()            # HO2. range: 1 - 10 *ppt*
#         c0[m.specmap["HO2H"], :] .= 0.01*rand()	               # HO2H range: 0.001 - 0.01 ppm
#         #c[7:11] = zeros(Float64,1,5)
#     c0[m.specmap["O2"], :] .= 0.21e6	               # O2: 0.21e6 ppm

#     tspan = (0, minutes) # minutes
#     prob = ODEProblem(m.rn, c0, tspan, m.p)
#     # 
    
#     # press, tempk, emission1, emission2, emission3, hv_shift
#     p_lower = [0.9, 288.0, 0.5, 0.5, 0.5, 0.0 ,0.0, 0.0, 0.0]
#     p_upper = [1.1, 308.0, 1.5, 1.5, 1.5, 2*pi , 2*pi, 2*pi, 1440.0]
#     # select set of Sobol-sampled parameter vectors
#     #Random.seed!(42)
#     para = Surrogates.sample(3750, p_lower, p_upper, SobolSample())[nruns * idx+1 : nruns * idx+nruns]
#     #println(para)

#     function prob_func(prob,i,repeat)
#         # prob.u0[:] = c0[:, i]
#         #press = 0.95 + rand() * 0.1                 # 0.9  -  1.1
#         #tempk = 298.0 + rand() * 20 - 10           # 288  -  308
#         #emission = 0.2 .+ rand(3) * 1.8              # 0.2  -  2
#         #hv_shift = rand()*1440                        # 0 - 1440
#         # println("emission\n")
#         m = ChemSys(para[i][1], para[i][2], collect(para[i][3:8]), para[i][9])
#         prob = ODEProblem(m.rn, c0, tspan, m.p)
#     end

#     ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

#     res = solve(ensemble_prob, Rosenbrock23(), trajectories=nruns, saveat=saveat,maxiters=Int(1e10), progress=true)
#     return res, para
# end


