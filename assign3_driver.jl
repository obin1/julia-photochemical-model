module mechanism
include("modlspc.jl")
include("assign2_rk.jl")
include("euler.jl")
include("assign3a.jl")
include("chemODE.jl")
include("stiff.jl")
#using DifferentialEquations
using LinearAlgebra
import Random
export gas

function gas(t1,t2,tm,days::Int64)
#=**************************************************************************
This is the main function in an ozone NOx formaldehyde mechanism, adapted from Fortran
by Obin Sturm (UC Davis) in October 2019.  It has been edited in order to generate training
data for machine learning and restricted inverse algorithms.

Notable changes:
Oct 2019 -- Random seeding of concentrations, looping over days, incorporating
            fluctuating daylight constant (HV) values.
         -- Option to use stiff equation solver.  To use this
            integration method, set method to "stiff!" instead of "euler!"

Mar 2020 -- Implementation of a NaN in the first column of S output matrix
            when using "euler" method.  This is used as a flag for when
            concentrations are negative and a mass balance violation occurs
            (a limitation of the original Fortran mechanism).


The original comments and file description are below.
     created by: Mike Kleeman (April 1999)
             University of California, Davis
             ECI 241 Air Quality Modeling
             Assignment #3
            The purpose of this program is to act as a driver for the evaluation and
            integration of a model ozone mechanism involving NOx and formaldehyde.
            Note: the "active" array contains active species in the first "maxact"
            locations and buildup species in the remaining "maxbo" locations.
**************************************************************************=#
# names for active species.
#    name = ["O3", "NO", "NO2","HCHO","HO2.","HO2H","HO.","O","HNO3","CO","H2"]

    Random.seed!(43) # Seed, for replicability. 43 was arbitrary 😁

    constn[1] = Float32(2.09e+05) # O2 concentration in ppm


# get the reaction rate constants at 1 atm and 25C--
    tempk = Float64(298.0)
    press = Float64(1.0)
    #      tempk = 283.
    #      press = 0.7
    global rk = getrk!(tempk,press,rk)


# call the integration subroutine to advance the solution
#    t1 = Float64(0.0)
#    t2 = Float64(60.0*24.0*1.0) # 60 mins/hr * 24 hrs/day *1 day
    dt = Float64(1.0e-3)
#    tm = Float64(6.0)
    C = zeros(1,11)
    S = zeros(1,10)
    J = zeros(1,1)

    Random.seed!(42)
    for d = 1:days
        println("Running ozone/NOx/formaldehyde mechanism: day ", d , "  of ", days)
        global c = zeros(Float64,1,maxact+maxbo)
        c[1] = 0.001*10^(2*rand())         # 03 range: 0.001 - 0.1 ppm
        c[2] = 0.0015*10^(2*rand())        # NO range: 0.0015 - 0.15
        c[3] = 0.0015*10^(2*rand())        # NO2 range: 0.0015 - 0.15
        c[4] = 0.02*10^(2*rand())          # HCHO range: 0.02 - 2 ppm
        c[5] = 1.0e-05 * rand()            # HO2. range: 1 - 10 *ppt*
        c[6] = 0.01*rand()	               # HO2H range: 0.001 - 0.01 ppm
        c[7:11] = zeros(Float64,1,5)

        # c[1] = 0.0  # For testing with original Fortran mechanism
        # c[2] = 0.0421
        # c[3] = 0.0679
        # c[4] = 0.1
        # c[5] = 0.0
        # c[6] = 0.01
        # c[7:11] = zeros(Float64,1,5)

        daymax = rand(Float64)

        # set daylight HV value
        constn[2] = 0.0


        global C_day,S_day,J_day = euler!(t1,t2,dt,tm,c,daymax,maxact+maxbo,rk,diff!)

        C = vcat(C,C_day)
        S = vcat(S,S_day)
        J = vcat(J,J_day)


    end
    C = C[2:end,:]
    S = S[2:end,:]
    J = J[2:end,:]
    return C,S,J
end
end
