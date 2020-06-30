function getrk!(tempk,press,rk)

#=**************************************************************************

 This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
 by Obin Sturm (UC Davis) in October 2019.
 The original comments and file description are below.

 created by: Mike Kleeman
             ECI 298 Air Quality Modeling
             April 1999

 The purpose of this subroutine is to calculate the reaction rate constants
 for the ozone NOx formaldehyde system.  The order of equations is assumed
 to be:
         1) NO2 + HV = NO + O
         2) O + O2 = O3
         3) O3 + NO = NO2 + O2
         4) HCHO + HV = 2 HO2. + CO
         5) HCHO + HV = H2 + CO
         6) HCHO + HO. = HO2. + CO + H2O
         7) HO2. + NO = HO. + NO2
         8) HO. + NO2 = HNO3
         9) HO2H + HV = 2 HO.
        10) HO2H + HO. = HO2. + H2O

 It is assumed that the photolytic constants are set at a constant value
 appropriate for mid-day in Southern California.

 Inputs:
  tempk - atmospheric temperature in degrees Kelvin
  press - atmospheric pressure in atmospheres

 Outputs:
  rk    - reaction rate constants in ppm min units

****************************************************************************=#


# function definition to convert from cm3 molec-1 s-1 to ppm min-1
    rkppmmin(rk1,press,tempk) = rk1 * 60 *  7.34e+15 * press / tempk

# begin by setting the rate constants for the photolytic reactions
    rk[1] = 0.5
    rk[4] = 0.015
    rk[5] = 0.022
    rk[9] = 0.0003

# reaction 2
    rk[2] = 6.0e-34 * ( tempk/300. )^(-2.3)
    rm = 7.34e+21 * press / tempk
    rk[2] = rk[2] * rm
    rk[2] = rkppmmin(rk[2],press,tempk)

# reaction 3
    rk[3] = 2.e-12 * exp(-1400/tempk)
    rk[3] = rkppmmin(rk[3],press,tempk)

# reaction 6
    rk[6] = 1.2e-14 * tempk^1 * exp(+287/tempk)
    rk[6] = rkppmmin(rk[6],press,tempk)

# reaction 7
    rk[7] = 3.7e-12 * exp(+250/tempk)
    rk[7] = rkppmmin(rk[7],press,tempk)

# reaction 8
# note: the following 2 lines are adapted from comments from the o.g. Fortran module
#      rk0 = 2.6e-30 * ( tempk/300 )^-3.2
#      rki = 2.4e-11 * ( tempk/300 )^-1.3
    rk0 = 2.0e-30 * ( tempk/300 )^(-3.0)
    rki = 2.5e-11 * ( tempk/300 )^0
    rm = 7.34e+21 * press / tempk
    rb = 0.6^(1.0 / ( 1.0 + (log10(rk0*rm/rki))^2.0 ))
    rk[8] = rk0 * rm * rb / (1.0 + rk0*rm/rki)
    rk[8] = rkppmmin(rk[8],press,tempk)

# reaction 10
    rk[10] = 3.3e-12 * exp(-200/tempk)
    rk[10] = rkppmmin(rk[10],press,tempk)

# return to the calling subroutine
    return rk
end
