# module JuliaPhotochem
using Measures
using DifferentialEquations
using Catalyst
using Random
using Plots, GraphRecipes
using CSV
using DataFrames


export JPMStruct, jpm
ENV["PATH"] = "/opt/homebrew/bin:" * ENV["PATH"] # to get the graphs.jl access
press=1.0
tempk=298.0
# Set the value of hv
hv = 1

# seed the random number generator
Random.seed!(42)


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

function bimolecular_uncertainty(f298,g,tempk)
    # calculate the uncertainty in the bimolecular rate constant
    # f298 is the rate constant at 298 K
    # g is additional uncertainty arising from temperature extrapolation
    # tempk is the temperature in K
    f = f298 * exp( abs( g * (1/tempk - 1/298.0) ) )
    return f
end


function julia_photochem(press=1.0::Float64, tempk=298.0::Float64)
    # function definition to convert from cm3 molec-1 s-1 to (ppm min)-1
	# 7.34e15 is Avogadro's number divided by R in L atm K-1 mol-1 divided by 1000 to convert to cm3 from L and 1e6 to convert to (ppm)-1
	# These units are valid for bimolecular reactions
    # rkppmmin(rk,press,tempk) = rk * 60.0 *  7.34e+15 * press / tempk 

    # make a function arrhenius(A,(E/R),tempk) to calculate the rate constant

    # Reaction rate constants

    # Old reaction 2 was O + O₂ → O₃
    # removed July 2024, JPL 19-5 page 434 Table 2-1
    # original was rkppmmin(6.0e-34 * ( tempk/300. )^(-2.3) * 7.34e+21 * press / tempk, press, tempk)
    # rk2(tempk,press) = molcm3frompressure(tempk,press)*6.1e-34*(tempk/298)^(-2.4)

    # Reaction 1 is NO₂ → NO + O₃
    # photolysis reaction, rate constant from original mechanism below
    # rk1(hv) = 0.5 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk1(hv) = 9.258186912960969E-003 * hv
    f1 = 0.0

    # Reaction 2 is O₃ + NO → NO₂ + O₂
    # updated July 2024, JPL 19-5 page 84 Table 1C 
    # original was rkppmmin(2.e-12 * exp(-1400/tempk), press, tempk)
    rk2(tempk) = arrhenius(3.0e-12,1500,tempk) 
    f2 = bimolecular_uncertainty(1.07,130.0,tempk)

    # Reaction 3 is HCHO + O₂ → HO₂ + CO
    # photolysis reaction, rate constant taken from original mechanism
    # rk3(hv) = 0.015 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk3(hv) = 3.442262702424394E-005 * hv
    f3 = 0.0

    # Reaction 4 is HCHO → H₂ + CO
    # photolysis reaction, rate constant from original mechanism
    # rk4(hv) = 0.022 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk4(hv) = 5.132233465281091E-005 * hv
    f4 = 0.0


    # Reaction 5 is HCHO + OH + O₂ → HO₂ + CO + H₂O
    # updated July 2024, JPL 19-5 page 106 Table 1D
    # OH + H₂CO (not HCHO or CH₂O, lol) → H₂O + HCO 
    # (formyl radical in oxygen makes HO₂ and CO, GEOS-Chem also parameterizes this as instantaneous)
    # original rate was rkppmmin(1.2e-14 * tempk^1 * exp(+287/tempk), press, tempk)
    rk5(tempk) = arrhenius(5.5e-12,-125.0,tempk)
    f5 = bimolecular_uncertainty(1.15,50.0,tempk)


    # Reaction 6 is HO₂ + NO → OH + NO₂
    # updated July 2024, JPL 19-5 page 83 Table 1C
    # original was rk7(press, tempk) = rkppmmin(3.7e-12 * exp(+250/tempk), press, tempk)
    rk6(tempk) = arrhenius(3.44e-12,-260,tempk)
    f6 = bimolecular_uncertainty(1.10,0.0,tempk) # yes that value for g of 0.0 was doublechecked with JPL 19-5 page 83 Table 1C

    # Reaction 7 is OH + NO₂ → HNO₃
    # updated July 2024, JPL 19-5 page 434 Table 2-1
    # the original rate is below
    # function rk8(press, tempk)
    #     rk0 = 2.0e-30 * ( tempk/300 )^(-3.0)
    #     rki = 2.5e-11 * ( tempk/300 )^0
    #     rm = 7.34e+21 * press / tempk
    #     rb = 0.6^(1.0 / ( 1.0 + (log10(rk0*rm/rki))^2.0 ))
    #     rk8 = rk0 * rm * rb / (1.0 + rk0*rm/rki)
    #     rk8 = rkppmmin(rk8,press,tempk)
    # end
    rk7(tempk,press) = termolecular(1.8e-30,3.0,2.8e-11,0,tempk,press)
    f7 = bimolecular_uncertainty(1.3,100.0,tempk) 

    # Reaction 8 is H₂O₂ → 2HO
    # this is a photolysis reaction, rate constant from original mechanism below
    # rk8(hv) = 0.0003 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk8(hv) = 7.373314014460711E-006 * hv
    f8 = 0.0

    # Reaction 9 is H₂O₂ + OH → H₂O + HO₂
    # updated July 2024, JPL 19-5 page 73
    # "the recommendation is a temperature independent value of 1.8e-12 cm3 molec-1 s-1"
    # original was rkppmmin(3.3e-12 * exp(-200.0/tempk),press,tempk)
    rk9() = 1.8e-12
    f9 = bimolecular_uncertainty(1.15,45.0,tempk) 

    # Reaction 10 is ALD2 + OH → MCO₃ + H₂O
    # updated July 2024, JPL 19-5 page 107 Table 1D
    # follows an Arrhenius equation with A = 4.63e-12, E/R = -350
    # note that GCARR_ac has the negative in it already so 
    # this is equivalent to GCARR_ac(4.63d-12, 350.0d0);
    rk10(tempk) = arrhenius(4.63e-12,-350.0,tempk)
    f10 = bimolecular_uncertainty(1.05,20.0,tempk) 

    # Reaction 11 is MGLY + hv → MCO₃ + CO + HO₂ 
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk11(hv) = 1.863860764162095E-004 * hv 
    f11 = 0.0


    # Reaction 12 is MCO₃ + NO₂ {+M} → PAN
    # updated July 2024, JPL 19-5 page 435 Table 2-1
    # MCO3 is peroxyacetyl radical CH3C(O)O2 in GEOS-Chem secret code
    rk12(tempk,press) = termolecular(7.3e-29,4.1,9.5e-12,1.6,tempk,press)
    f12 = bimolecular_uncertainty(1.05,0.0,tempk) # yes that value for g of 0.0 was doublechecked with JPL 19-5 page 435 Table 2-1 

    # Reaction 13 is thermal decomposition of PAN → MCO₃ + NO₂ 
    # updated July 2024, JPL 19-5 page 521 Table 3-1
    # using the equilibrium constant divided by the same termolecular function as reaction 12
    # to get the backward rate constant
    rk13(tempk,press) = termolecular(7.3e-29,4.1,9.5e-12,1.6,tempk,press) / ( 9.0e-29*exp(14000/tempk) )
    f13 = f12 * bimolecular_uncertainty(1.2,2000.0,tempk) 

    # Reaction 14 is PAN + OH → HCHO + CO + NO₂ + H₂O
    # updated July 2026, from GEOS-Chem v14.7.0
    # 3.00d-14 is less than the JPL 19-5 recommendation <4.0d-14 in reaction D46
    rk14() = 3.0e-14
    f14 = 0.0

    # Reaction 15 is O3 + OH → HO2 + O2
    # updated July 2026, JPL 19-5 page 67 table 1B reaction B6
    rk15(tempk) = arrhenius(1.7e-12,940.0,tempk)
    f15 = bimolecular_uncertainty(1.15,50.0,tempk)



    # Evaluate the rate constant functions to get their values
    p = [rk1(hv), rk2(tempk), rk3(hv), rk4(hv), rk5(tempk),
        rk6(tempk), rk7(press,tempk), rk8(hv), rk9(), rk10(tempk), 
        rk11(hv), rk12(tempk,press), rk13(tempk,press), rk14(), rk15(tempk)]
    f_uncertainty = [
        f1, f2, f3, f4, f5,
        f6, f7, f8, f9, f10,
        f11, f12, f13, f14, f15
    ]


    rn = @reaction_network begin
        rk1/O₂, NO₂ + O₂ → NO + O₃
        rk2, O₃ + NO → NO₂ + O₂
        rk3/(O₂^2), HCHO + 2O₂ → 2HO₂ + CO  # rate law not proportional to O₂
        rk4, HCHO → H₂ + CO
        rk5/O₂, HCHO + OH + O₂ → HO₂ + CO + H₂O # rate law not proportional to O₂
        rk6, HO₂ + NO → OH + NO₂
        rk7, OH + NO₂ → HNO₃
        rk8, H₂O₂ → 2OH
        rk9, H₂O₂ + OH → H₂O + HO₂
        rk10/O₂, ALD2 + OH + O₂ → MCO₃ + H₂O  # rate law not proportional to O₂
        rk11/(O₂^2), MGLY + 2O₂ → MCO₃ + CO + HO₂  # rate law not proportional to O₂
        rk12, MCO₃ + NO₂ → PAN
        rk13, PAN → MCO₃ + NO₂
        rk14, PAN + OH → HCHO + CO + NO₂ + H₂ + O₂
        rk15, O₃ + OH → HO₂ + O₂
        # All species go to null for testing
        # rk2*1e12, NO₂ → ∅;
        # rk5*1e12, NO → ∅;
        # rk6*1e12, O₃ → ∅;
        # rk7*1e12, OH → ∅;
        # rk9*1e12, HO₂ → ∅;
        # rk2*1e12, HCHO → ∅;
        # rk6*1e12, H₂O₂ → ∅;
        # rk7*1e12, HNO₃ → ∅;
        # rk9*1e12, CO → ∅;
        # rk12*1e12, H₂ → ∅;
        # rk13*1e12, H₂O → ∅;
        # rk13*1e12, O₂ → ∅;
        # rk12*1e12, MCO₃ → ∅;
        # rk13*1e12, PAN → ∅;
        # rk13*1e12, ALD2 → ∅;
        # rk13*1e12, MGLY → ∅;
    end rk1 rk2 rk3 rk4 rk5 rk6 rk7 rk8 rk9 rk10 rk11 rk12 rk13 rk14 rk15

    nspecs = size(species(rn),1)

    specmap = Dict()
    for (k, v) in speciesmap(rn)
        specmap[replace(string(k), r"\(t\)$"=>"")] = v # remove '(t)' from names
    end
    natoms = 4 # C N H O
    atoms = zeros(Float64, nspecs, natoms)
    atoms[specmap["O₃"],   :]   = [0 0 0 3]
    atoms[specmap["NO"],   :]   = [0 1 0 1]
    atoms[specmap["NO₂"],  :]   = [0 1 0 2]
    atoms[specmap["HCHO"], :]   = [1 0 2 1]
    atoms[specmap["HO₂"],  :]   = [0 0 1 2]
    atoms[specmap["H₂O₂"], :]   = [0 0 2 2]
    atoms[specmap["OH"],   :]   = [0 0 1 1]
    atoms[specmap["HNO₃"], :]   = [0 1 1 3]
    atoms[specmap["CO"],   :]   = [1 0 0 1]
    atoms[specmap["H₂"],   :]   = [0 0 2 0]
    atoms[specmap["ALD2"], :]   = [2 0 4 1]
    atoms[specmap["MGLY"], :]   = [3 0 4 2]
    atoms[specmap["MCO₃"], :]   = [2 0 3 3]
    atoms[specmap["PAN"],  :]   = [2 1 3 5]
    atoms[specmap["H₂O"],  :]   = [0 0 2 1]
    atoms[specmap["O₂"],   :]   = [0 0 0 2]

    return rn, atoms, specmap, p, f_uncertainty
end

struct JPMStruct
    rn
    atoms::Array{Float64, 2}
    specmap::Dict
    p::Array{Float64, 1}
    f_uncertainty::Array{Float64, 1}
end


function ppb_to_molec_cm3(ppb::Float64, tempk::Float64, press::Float64)
    molec_cm3 = ppb / 1e9 * press / .0821 / tempk * 6.022e23 / 1e3
    return molec_cm3
end
function molec_cm3_to_ppb(molec_cm3::Float64, tempk::Float64, press::Float64)
    ppb = molec_cm3 / 6.022e23 * 1e3 * .0821 * tempk / press * 1e9
    return ppb
end


JPMStruct(press=1.0::Real, tempk=298.0::Real) = JPMStruct(julia_photochem(press, tempk)...)
jpm = JPMStruct()
graph = Graph(jpm.rn)
# save the graph to a file
# savegraph(graph,"jpm-graph.pdf","pdf")
savegraph(graph,"jpm-graph.png","png")

function plot_trajectories(jpm, tempk, press, vars_to_plot=["O₃", "NO", "NO₂", "HCHO", "HO₂", "H₂O₂", "OH", "HNO₃", "CO", "H₂", "ALD2", "MGLY", "MCO₃", "PAN", "H₂O", "O₂"])
    tempk = 283.0 + 20.0*rand() # random temperature between 283 and 303 K
    jpm = JPMStruct(press, tempk)

    # Initialize c0 with zeros
    c0 = zeros(Float64, length(jpm.specmap))
    # Set specific initial conditions
    c0[jpm.specmap["O₂"], :] .= 0.21e9  # O₂: 0.21e9 ppb

    # Initialize concentrations from uniform distribution
    c0[jpm.specmap["O₃"], :]   .= 50*rand()         
    c0[jpm.specmap["NO"], :]   .= 20*rand()   
    c0[jpm.specmap["NO₂"], :]  .= 20*rand()
    c0[jpm.specmap["HCHO"], :] .= 50*rand()  
    c0[jpm.specmap["HO₂"], :]  .= 0.01*rand()
    c0[jpm.specmap["H₂O₂"], :] .= 10*rand()
    c0[jpm.specmap["OH"], :]   .= 0.001*rand() 
    c0[jpm.specmap["MCO₃"], :] .= 0.01*rand()
    c0[jpm.specmap["PAN"], :]  .= 25*rand()
    c0[jpm.specmap["ALD2"], :] .= 10*rand()
    c0[jpm.specmap["MGLY"], :] .= 50*rand()
    # Convert concentrations to molecules/cm^3
    c0 = ppb_to_molec_cm3.(c0, tempk, press)

    # Define the time span for the simulation
    tspan = (0, 900*2)
    # Create the ODE problem
    prob = ODEProblem(jpm.rn, c0, tspan, jpm.p,combinatoric_ratelaws = false)

    # Solve the ODE problem
    sol = solve(prob, Rosenbrock23(), saveat=1.0, progress=true)

    # Get indices of variables to plot
    global indices = [jpm.specmap[var] for var in vars_to_plot]
    # Extract the variables at every time step
    selected_vars = [state[indices] for state in sol.u]
    # Transpose the matrix to get each variable as a separate vector
    selected_vars_transposed = transpose(hcat(selected_vars...))
    # Convert to ppb
    selected_vars_transposed = molec_cm3_to_ppb.(selected_vars_transposed, tempk, press)

    # Plot each variable as its own subplot
    num_vars = length(vars_to_plot)
    fig_height = 200 * num_vars  # Adjust height per variable
    fig_width = 800  # Standard width
    p = plot(layout = (num_vars, 1), size = (fig_width, fig_height))
    # Add each variable to the plot
    for (i, var) in enumerate(vars_to_plot)
        plot!(p[i], sol.t, selected_vars_transposed[:,i], legend=false, left_margin=6Measures.mm)
        ylims!(p[i], (0, maximum(selected_vars_transposed[:,i])))
        xlabel!(p[i], "Time [seconds]")
        ylabel!(p[i], var * " [ppb]")
        title!(p[i], var)
    end
    display(p)
end

function plot_selected_vars_transposed(selected_vars_transposed, vars_to_plot, sol_t, tempk, press)
    # Plot each variable as its own subplot
    num_vars = length(vars_to_plot)
    fig_height = 200 * num_vars  # Adjust height per variable
    fig_width = 800  # Standard width
    p = plot(layout = (num_vars, 1), size = (fig_width, fig_height))
    # Add each variable to the plot
    for (i, var) in enumerate(vars_to_plot)
        plot!(p[i], sol_t, selected_vars_transposed[:,i], legend=false, left_margin=6Measures.mm)
        # ylims!(p[i], (0, maximum(selected_vars_transposed[:,i])))
        xlabel!(p[i], "Time [seconds]")
        ylabel!(p[i], var * " [ppb]")
        title!(p[i], var)
    end
    display(p)
end

plot_trajectories(jpm, tempk, press)

z = rand(Float64)*2000 # up to 2 km
function barometric(z::Float64)
    # This function calculates the pressure and temperature at a given altitude
    # assuming an vertically averaged temperature of 288.2 
    # and a dry adiabatic lapse rate of 9.8e-3 ˚C/m
    tempk = Float64(298.0) - Float64(9.8e-3)*z   # assume dry adiabatic lapse rate of 9.8e-3 ˚C/m
    press = Float64(1.0)*exp(-z*9.81/8.31*29e-3/288.2) #  hypsometric equation, assume average temp of 298 + (298-9.8*2)
    return tempk, press
end


# function to run N simulations of 1 hour, with random initial conditions, sampling every 5 minutes
function run_simulations(num_tests::Int, duration::Float64, sampling_interval::Float64)
    # num_tests = 1
    # duration = 3600.0
    # sampling_interval = 300.0
    # print number of experiments

    # create JPM based on random temperature conditions

    
    indices = [jpm.specmap[var] for var in ["O₃", "NO", "NO₂", "HCHO", "HO₂", "H₂O₂", "OH", "HNO₃", "CO", "H₂", "ALD2", "MGLY", "MCO₃", "PAN", "H₂O", "O₂"]]
    label = ["Experiment number", "Temperature (K)", "Time (min)", "O3 [ppb]", "NO [ppb]", "NO2 [ppb]", "HCHO [ppb]", "HO2 [ppb]", "H2O2 [ppb]", "OH [ppb]", "HNO3 [ppb]", "CO [ppb]", "H2 [ppb]", "ALD2 [ppb]", "MGLY [ppb]", "MCO3 [ppb]", "PAN [ppb]", "H2O [ppb]", "O2 [ppb]"]
    num_time_points = Int(duration / sampling_interval)+1
    experiments_array = Array{Any}(undef, (num_tests * num_time_points + 1), length(label))
    experiments_array[1,:] = label

    row_index = 2
    for i in 1:num_tests
        tempk = 283.0 + 20.0*rand() # random temperature between 283 and 303 K
        jpm = JPMStruct(press, tempk) # initialize the mechanism with a new temperature
        println("Experiment number: ", i)
        # Initialize c0 with zeros
        c0 = zeros(Float64, length(jpm.specmap))
        # Set specific initial conditions
        c0[jpm.specmap["O₂"], :] .= 0.21e9  # O₂: 0.21e6 ppm
        # Initialize concentrations from uniform distribution
        c0[jpm.specmap["O₃"], :]   .= 50*rand()         
        c0[jpm.specmap["NO"], :]   .= 20*rand()   
        c0[jpm.specmap["NO₂"], :]  .= 20*rand()
        c0[jpm.specmap["HCHO"], :] .= 50*rand()  
        c0[jpm.specmap["HO₂"], :]  .= 0.01*rand()
        c0[jpm.specmap["H₂O₂"], :] .= 10*rand()
        c0[jpm.specmap["OH"], :]   .= 0.001*rand() 
        c0[jpm.specmap["MCO₃"], :] .= 0.01*rand()
        c0[jpm.specmap["PAN"], :]  .= 25*rand()
        c0[jpm.specmap["ALD2"], :] .= 10*rand()
        c0[jpm.specmap["MGLY"], :] .= 50*rand()
        # Convert concentrations to molecules/cm^3
        c0 = ppb_to_molec_cm3.(c0, tempk, press)
        # Define the time span for the simulation
        tspan = (0, duration)
        # Create the ODE problem
        # prob = ODEProblem(jpm.rn, c0, tspan, jpm.p)
        prob = ODEProblem(jpm.rn, c0, tspan, jpm.p, combinatoric_ratelaws = false)
        # Solve the ODE problem
        global sol = solve(prob, Rosenbrock23(), saveat=sampling_interval, progress=false)
        # Extract the variables at every time step
        selected_vars = [state[indices] for state in sol.u]
        # create array
        global selected_vars_transposed = transpose(hcat(selected_vars...))
        # Convert to ppb
        selected_vars_transposed = molec_cm3_to_ppb.(selected_vars_transposed, tempk, press)
        for (t, time_point) in enumerate(0:sampling_interval:duration)
            # print the time point
            # println("Time point: ", time_point)
            experiments_array[row_index, :] = [i, tempk, time_point, selected_vars_transposed[t, :]...]
            row_index += 1
        end
    end
    return experiments_array
end

# number_of_experiments = Integer(1e6)
number_of_experiments = Integer(1000)
# number_of_experiments = 3
experiments_array = run_simulations(number_of_experiments, 3600.0, 300.0)

# mass balance checker
M = zeros(Int, 16, 4)
M[1,  :] = [0, 0, 0, 3]  # O₃
M[2,  :] = [0, 1, 0, 1]  # NO
M[3,  :] = [0, 1, 0, 2]  # NO₂
M[4,  :] = [1, 0, 2, 1]  # HCHO
M[5,  :] = [0, 0, 1, 2]  # HO₂
M[6,  :] = [0, 0, 2, 2]  # H₂O₂
M[7,  :] = [0, 0, 1, 1]  # OH
M[8,  :] = [0, 1, 1, 3]  # HNO₃
M[9,  :] = [1, 0, 0, 1]  # CO
M[10, :] = [0, 0, 2, 0]  # H₂
M[11, :] = [2, 0, 4, 1]  # ALD2
M[12, :] = [3, 0, 4, 2]  # MGLY
M[13, :] = [2, 0, 3, 3]  # MCO₃
M[14, :] = [2, 1, 3, 5]  # PAN
M[15, :] = [0, 0, 2, 1]  # H₂O
M[16, :] = [0, 0, 0, 2]  # O₂

# get the tendencies for all experiments
tendencies_diff = diff(experiments_array[2:end, 4:end], dims=1)
# skip every 12 rows to get the tendencies for each time point
tendencies = tendencies_diff[1:13:end, :]
maximum(tendencies*M, dims=1)


# save the experiments array to a csv file
CSV.write("1000experiments_JPMv1.2_varyingtemperature.csv", 
          DataFrame(experiments_array[2:end, :], 
          Symbol.(experiments_array[1, :])))

# convert csv array to selected_vars_transposed
# selected_vars_transposed = experiments_array[2:end, 4:end]
# time_points = experiments_array[2:end, 3]
# plot_selected_vars_transposed(selected_vars_transposed, ["O₃", "NO", "NO₂", "HCHO", "HO₂", "H₂O₂", "OH", "HNO₃", "CO", "H₂", "ALD2", "MGLY", "MCO₃", "PAN", "H₂O", "O₂"], time_points, tempk, press)


# end

# B = netstoichmat(jpm.rn)
# B = B[indices,:]
# use B to get an edgelist. Every row indexes a species. Every column index is a reaction
# The value in the matrix is the stoichiometric coefficient of the species in the reaction
# The edgelist is a list of tuples (species, reaction, stoichiometric coefficient)
# edgelist = []
# species_names = species(jpm.rn)
# for j in 1:size(B, 2)
#     for i in 1:size(B, 1)
#         species_name = replace(string(species_names[i]), "(t)" => "")
#         # replace ₂ with 2 and same for 3
#         species_name = replace(species_name, "₂" => "2")
#         species_name = replace(species_name, "₃" => "3")
    
#         if B[i, j] < 0
#             push!(edgelist, (i, species_name, "R"*string(j), abs(B[i, j])))
#         elseif B[i, j] > 0
#             push!(edgelist, (i, "R"*string(j), species_name, abs(B[i, j])))
#         end
#     end
# end

# println(edgelist)

# save the edgelist to a csv file
# with the header:
# header = ["species index (starts from 1)","from","to","stoichiometric value "]
# Create the DataFrame without specifying the header
# df = DataFrame(edgelist)

# Rename the columns to the desired string headers
# rename!(df, header)
# CSV.write("jpm_edgelist.csv", df)

# rewrite rk1,..., as REPL-scope values from jpm.p

