#=**************************************************************************

This script imports and runs an ozone NOx formaldehyde mechanism, producing 3
text files to be used in training physics constrained neural networks.

Written in Julia Version 1.2.0, tested up to 1.4.2

Written by Obin Sturm (posturm@ucdavis.edu)


The inverse problem tries to solve for S according to

  A*S = ∆C

 where S is the integrated reaction rates over a given timestep, and ∆C is the
 change in concentration.

The output files are "C.txt", "S.txt" and "delC.txt" for use in machine
learning research and inverse problems.  This was first used in a paper
introducing a framework for conserving mass in machine learning applications
for air quality and climate modeling: https://doi.org/10.5194/gmd-13-4435-2020

 C and delC correspond to the following species:
         1) O3
         2) NO
         3) NO2
         4) HCHO
         5) HO2
         6) HO2H
         7) HO.
         8) O
         9) HNO3
        10) CO
        11) H2

 S corresponds to the following reactions:
         1) NO2 + HV = NO + O
         2) O + O2 = O3
         3) O3 + NO = NO2 + O2
         4) HCHO + HV = 2 HO2. + CO
         5) HCHO + HV = H2 + CO
         6) HCHO + HO. = HO2. + CO + H2O
         7) HO2. + NO = HO. + NO2
         8) HO. + NO2 = HNO3
         9) HO2H + HV = 2 HO.
        10) HO2H + HO. = H2O + HO2
The stoichiometry matrix A, which can also be interpreted as the directed,
weighted incidence matrix of the species-reaction network, is

A =           [0  1 -1  0  0  0  0  0  0  0;
               1  0 -1  0  0  0 -1  0  0  0;
              -1  0  1  0  0  0  1 -1  0  0;
               0  0  0 -1 -1 -1  0  0  0  0;
               0  0  0  2  0  1 -1  0  0  1;
               0  0  0  0  0  0  0  0 -1 -1;
               0  0  0  0  0 -1  1 -1  2 -1;
               1 -1  0  0  0  0  0  0  0  0;
               0  0  0  0  0  0  0  1  0  0;
               0  0  0  1  1  1  0  0  0  0;
               0  0  0  0  1  0  0  0  0  0]
****************************************************************************=#

mechanism = include("assign3_driver.jl") # Include photochemical mechanism module
using DelimitedFiles

mechanism = include("assign3_driver.jl")
C,S,J = mechanism.gas(0.0,24.0*60.0,6.0,5000) # Run 5000 full days, sampling every 6 minutes
# delC = diff([C ; zeros(1,11)],dims = 1)

  open("S.txt", "w") do io
                         writedlm(io, S,',')
                   end
  open("C.txt", "w") do io
                          writedlm(io, C,',')
                   end
  open("J.txt", "w") do io
                         writedlm(io, J,',')
                   end
