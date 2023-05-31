module mechanism
include("modlspc.jl")
include("assign2_rk.jl")
include("euler.jl")
include("assign3a.jl")
using LinearAlgebra
import Random
export gas

function gas(t1,t2,tm,days::Int64)
#=**************************************************************************
This is the top level function in an ozone NOx formaldehyde mechanism, rewritten from Fortran
by Obin Sturm (UC Davis) in October 2019.  It has been adapted in order to generate training
data for machine learning and restricted inverse algorithms.
Written in Julia Version 1.2.0, tested up to 1.4.2

For questions please contact Obin Sturm at posturm@ucdavis.edu

Notable changes:
Oct 2019 -- Random seeding of concentrations, looping over days, incorporating
            fluctuating daylight constant (HV) values.
         -- Option to use stiff equation solver.  To use this
            integration method, set method to "stiff!" instead of "euler!"
            Note: deprecated March 2020

Mar 2020 -- Implementation of a NaN in the first column of S output matrix
            when using "euler" method.  This is used as a flag for when
            concentrations are negative and a mass balance violation occurs
            (a limitation of the original Fortran mechanism).
            Note: removed January 2021

July 2021 -- Fix of a printing bug inherited from original Fortran code,
            that printed pseudo-steady state species as zero. This is fixed
            in euler.jl

            Inputs:
             t1   - start time in minutes
             t2   - stop time in minutes
             tm   - concentration recording timestep
             days - number of independent trials from t1 to t2
            Outputs:
             C    -  species concentrations, each column a new species
             S    -  integrated rate of reaction, each column a new reaction
             J    -  daylight constant influencing photolytic reactions

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
    T = zeros(1,1)
    P = zeros(1,1)
    K = zeros(1,10)

    Random.seed!(42)
    for d = 1:days
        println("Running ozone/NOx/formaldehyde mechanism: day ", d , "  of ", days)

        global c = zeros(Float64,1,maxact+maxbo)
        # c[1] = 0.001*10^(2*rand())         # 03 range: 0.001 - 0.1 ppm
        # c[2] = 0.0015*10^(2*rand())        # NO range: 0.0015 - 0.15
        # c[3] = 0.0015*10^(2*rand())        # NO2 range: 0.0015 - 0.15
        # c[4] = 0.02*10^(2*rand())          # HCHO range: 0.02 - 2 ppm
        # c[5] = 1.0e-05 * rand()            # HO2. range: 1 - 10 *ppt*
        # c[6] = 0.01*rand()	               # HO2H range: 0.001 - 0.01 ppm

        # Initialize concentrations from uniform distribution, ranges per SI of https://doi.org/10.1029/2020JD032759
        c[1] = 100*rand()*1e-3            
        c[2] = 10*rand()*1e-3              
        c[3] = 10*rand()*1e-3             
        c[4] = 20*rand()*1e-3              
        c[5] = 0.5*rand()*1e-3             
        c[6] = 10*rand()*1e-3   
        # Randomly initialize z, which changes P,T which change rk
        z = rand(Float64)*2000 # up to 2 km
        tempk = Float64(298.0) - Float64(9.8e-3)*z   # assume dry adiabatic lapse rate of 9.8e-3 ˚C/m
        press = Float64(1.0)*exp(-z*9.81/8.31*29e-3/288.2) #  hypsometric equation, assume average temp of 298 + (298-9.8*2)

        
        # varyTP = rand(Float64)
        # tempk = Float64(273.0)+25*varyTP
        # press = Float64(0.7)+0.3*varyTP
        #      tempk = 283.
        #      press = 0.7
        global rk = getrk!(tempk,press,rk)
        # println("rk: ", rk)
        # println("P = ",press,", T = ",tempk)

        # make sure ss species are also printed at the beginning
        global s, r, fr, rlr = difun!(constn,c,rpssa,rk,r,fr,rlr)
        # print("steady state species at beginning: ",s)
        c[7] = s[1]
        c[8] = s[2]
        c[9:11] = zeros(Float64,1,3)       # buildup species start at zero
        # c[11] = 0.0011637709108759861 # check H2 invariance for Ziming

         # c[1] = 0.0  # For testing with original Fortran mechanism
         # c[2] = 0.0421
         # c[3] = 0.0679
         # c[4] = 0.1
         # c[5] = 0.0
         # c[6] = 0.01
         # c[7:11] = zeros(Float64,1,5)

        daymax = rand(Float64)

        # set daylight HV value
        constn[2] = 1.0 # constn[2] = 0.0

        global C_day,S_day,J_day = euler!(t1,t2,dt,tm,c,daymax,maxact+maxbo,rk,diff!)

        C = vcat(C,C_day)
        S = vcat(S,S_day)
        J = vcat(J,J_day)
        T = vcat(T,tempk)
        P = vcat(P,press)
        K = vcat(K,rk[:,1:10])


    end
    C = C[2:end,:]
    S = S[2:end,:]
    J = J[2:end,:]
    T = T[2:end,:]
    P = P[2:end,:]
    K = K[2:end,:]
    return C,S,J,T,P,K
end
end
