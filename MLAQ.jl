#=**************************************************************************
This is a wrapper for a photochemical mechanism used to generate training
data for machine learning and restricted inverse algorithms.  It loads the
module "FormO3NOx" as mechanism.gas.  An adjustable parameter is number of days
to run the mechanism.

last edited June 30 2020
Obin Sturm posturm@ucdavis.edu
**************************************************************************=#

#mechanism = include("/Users/obin/Dropbox/ML Share/Photochemical Mechanisms/FormO3NOx/assign3_driver.jl
# For when in the same path
mechanism = include("assign3_driver.jl")
using DelimitedFiles

# The system ∆C = AS
A_reduced = [ 0  1 -1  0  0  0  0  0  0  0;  # O3
              1  0 -1  0  0  0 -1  0  0  0;  # NO
             -1  0  1  0  0  0  1 -1  0  0;  # NO2
              0  0  0 -1 -1 -1  0  0  0  0;  # HCHO
              0  0  0  2  0  1 -1  0  0  1;  # HO2.
              0  0  0  0  0  0  0  0 -1 -1];  # HO2H
#              0  0  0  0  0 -1  1 -1  2 -1;  # OH.
#              1 -1  0  0  0  0  0  0  0  0;  # O
#              0  0  0  0  0  0  0  1  0  0;  # HNO3
#              0  0  0  1  1  1  0  0  0  0;  # CO
#              0  0  0  0  1  0  0  0  0  0]; # H2


# Change this to adjust the number of full days the mechanism should run for
numdays = 365

# P
C,S,J = mechanism.gas(0.0,24.0*60,6.0,numdays)

# Create a data structure for day, time, daylight constant, and concentrations
Cday = repeat(1:numdays, inner = 241)
Ctime = repeat(collect(0.0:6.0:1440.0), outer = numdays)
Clabel = permutedims(["Day","Time [min]","Daylight Constant [0-1]","O3","NO","NO2","HCHO","HO2.","HO2H","OH.","O","HNO3","CO","H2"])
C_and_J = vcat(Clabel, hcat(Cday,Ctime,J,C))

# Create a data structure for S, related to adjacent C's as A*S = (C(t+∆t)-C(∆t))
Sday = repeat(1:numdays, inner = 240)
Stime = repeat(collect(0.0:6.0:1434.0), outer = numdays)
Slabel = permutedims(["Day","Time [min]","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10"])
S = vcat(Slabel, hcat(Sday,Stime,S))

open("C_and_J.txt", "w") do io
                        writedlm(io, C,',')
                  end
open("S.txt", "w") do io
                        writedlm(io, S,',')
                  end
