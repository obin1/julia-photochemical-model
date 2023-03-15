#=**************************************************************************

This script imports an abbreviated photochemical mechanism involving
ozone, formaldehyde, and NOx.

The original version was written in Julia Version 1.2.0, 
and subsequent versions have been tested up to 1.8.1

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
        12) H2O

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
               0  0  0  0  1  0  0  0  0  0;
               0  0  0  0  0  1  0  0  0  1]
****************************************************************************=#

mechanism = include("assign3_driver.jl") # Include photochemical mechanism module
include("visualizeCQ.jl")
using DelimitedFiles
using Plots
using StatsBase
using Measures # for cleaner figure margins


# Note: currently configured for Ziming, see photolytic constant set on:
        # euler.jl line 53
        # assign3_driver.jl line 106
num_exp = 1000
C,S,J,T,P,K = mechanism.gas(0.0,1*20.0,Float64(1.0e-2),num_exp) # mechanism.gas(0.0,24.0*60.0,60.0,5000) # Run 5000 full days, sampling every 6 minutes

Cday = repeat(1:num_exp, inner = 2001)
Ctime = repeat(collect(0.0:1.0e-2:20.0), outer = num_exp)
Clabel = permutedims(["Day","Time [min]","T [K]","P [atm]","O3","NO","NO2","HCHO","HO2.","HO2H","OH.","O","HNO3","CO","H2","H2O"])
X = vcat(Clabel, hcat(Cday,Ctime,repeat(T[:,1],inner=2001),repeat(P[:,1],inner=2001),C))

# To write model output to CSV
#   open("X_Ziming.txt", "w") do io
#                           writedlm(io, X,',')
#                    end

# delC = diff([C ; zeros(1,11)],dims = 1)

#   open("S_Ziming.txt", "w") do io
#                          writedlm(io, S,',')
#                    end
#   open("C_Ziming.csv", "w") do io
#                           writedlm(io, C,',')
#                    end
#   open("J_Ziming.txt", "w") do io
#                          writedlm(io, J,',')
#                    end
# V_old = readdlm("V.txt",' ')     


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
               0  0  0  0  1  0  0  0  0  0;
               0  0  0  0  0  1  0  0  0  1]


# V = vcat([ 3.49497182e-11, -5.44313032e-11, -1.94818953e-11, -7.07106781e-01,
#          5.24143203e-11,  3.49495370e-11,  1.74603652e-11,  3.49498516e-11,
#         -2.01927364e-12, -7.07106781e-01, -1.26962632e-11]',
#         [ 4.32105740e-13,  5.77350269e-01,  5.77350269e-01, -3.07830081e-11,
#          6.47877585e-13,  4.32966163e-13,  2.16420631e-13,  4.32237579e-13,
#          5.77350269e-01, -3.12153220e-11, -1.56357566e-13]',
#         [-0.37014454,  0.30845379, -0.06169076, -0.18507227, -0.55521681,
#         -0.37014454, -0.18507227, -0.37014454, -0.24676303,  0.18507227,
#          0.13459802]')

# CQ3_snapped = [-6, 5, -1, -3, -9, -6, -3, -6, -4, 3, 2.213115]
# CQ3_snapped = [-6, 5, -1, -3, -9, -6, -3, -6, -4, 3, 2]


# experiment = 34 #21 #3 might be good 30 and 32 is interesting
# p1 = visualizeCQ(C[2001*(experiment-1)+1:2001*experiment,:],V,
#                 title="P="*string(round(100*P[experiment])/100)*"atm, T="*string(round(T[experiment]-273.15))*"˚C")


# experiment = 60 #40, 35 might be good
# p2 = visualizeCQ(C[2001*(experiment-1)+1:2001*experiment,:],V,
#                  title="P="*string(round(100*P[experiment])/100)*"atm, T="*string(round(T[experiment]-273.15))*"˚C",
#                  legend=false)


plot(p1, p2, layout=(1,2))
savefig("CQ3_varyingTP.pdf")

HC = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0]
HN = [0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0]
HH = [0, 0, 0, 2, 1, 2, 1, 0, 1, 0, 2, 2]
HH[7] = 0 # excluding OH, hydrogen is much better conserved

# CV's are coefficients of variation (standard deviation divided by the mean)
cvC = zeros(num_exp)
cvN = zeros(num_exp)
cvH = zeros(num_exp)
cv1 = zeros(num_exp)
cv2 = zeros(num_exp)
cv3 = zeros(num_exp)
cv3_snapped = zeros(num_exp)

for i = 1:num_exp
        c_exp = X[X[:,1].==i,5:end] # C[2001*(i-1)+1:2001*i,:] #
        cvC[i] = std(c_exp*HC) / abs(mean(c_exp*HC))
        cvN[i] = std(c_exp*HN) / abs(mean(c_exp*HN))
        cvH[i] = std(c_exp*HH) / abs(mean(c_exp*HH))
        # cv1[i] = std(c_exp*V[1,:]) / abs(mean(c_exp*V[1,:]))
        # cv2[i] = std(c_exp*V[2,:]) / abs(mean(c_exp*V[2,:]))
        # cv3[i] = std(c_exp*V[3,:]) / abs(mean(c_exp*V[3,:]))
        # cv3_snapped[i] = std(c_exp*CQ3_snapped) / abs(mean(c_exp*CQ3_snapped))
end

# Maximum cv for each CQ, as well as 95 percentile, in percent
# cv3_max = maximum(cv3)*100
# cv3_95 = percentile(cv3,95)*100
# cv3_50 = percentile(cv3,50)*100

h1 = histogram(cv1*1e2,xlabel = "coefficient of variation [%]", ylabel = "number of cases",title="CQ1")
h2 = histogram(cv2*1e2,xlabel = "coefficient of variation [%]", ylabel = "number of cases",title="CQ2")
h3 = histogram(cv3*1e2,xlabel = "coefficient of variation [%]", ylabel = "number of cases",title="CQ3")
hall = plot(h1,h2,h3,size=(12e2,7e2))
plot!(margin=5mm)
savefig("histogram_relativevariation.pdf")

