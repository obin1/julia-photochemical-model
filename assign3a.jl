# Contained functions in file: bldup!, difun!

function bldup!(r,dy)
#=**************************************************************************
     This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
     by Obin Sturm (UC Davis) in October 2019.

     This module was recieved from Mike Kleeman (UC Davis) in October 2018.
     The original comments and file description are below.

     MODEL DEPENDENT SUBROUTINE TO CALCULATE FORMATION
     RATES OF NON-REACTING SPECIES

     Inputs:
     r   - reaction rates for the mechanism (get from difun)
         1) NO2 + HV = NO + O
         2) O + O2 = O3
         3) O3 + NO = NO2 + O2
         4) HCHO + HV = 2 HO2. + CO
         5) HCHO + HV = H2 + CO
         6) HCHO + HO. = HO2. + CO + H2O
         7) HO2. + NO = HO. + NO2
         8) HO. + NO2 = HNO3
         9) HO2H + HV = 2 HO.
        10) HO2H + HO. = H2O + HO2.

      Outputs:
         dy  - formation rates of the buildup species [ppm min-1]
**************************************************************************=#
    # HNO3
    dy[1] = r[8]

    # CO
    dy[2] = r[4]+r[5]+r[6]

    # HO
    dy[3] = r[5]

    return dy
end

function  difun!(c,a,s,rk,r,fr,lr)
    #=**************************************************************************
    This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
    by Obin Sturm (UC Davis) in October 2019.

    This module was recieved from Mike Kleeman (UC Davis) in October 2018.
    The original comments and file description are below.
    MODEL DEPENDENT SUBROUTINE TO CALCULATE RATES OF
    FORMATION AND LOSS FOR ACTIVE SPECIES

        dA/dt = fr - lr*a

    Inputs
     c   - constant species
            c(1) = O2 (set to atmospheric value = 2.090E+05ppm)
            c(2) = HV (leave set to 1)
     a   - active species concentrations
            a(1) = O3
            a(2) = NO
            a(3) = NO2
            a(4) = HCHO
            a(5) = HO2.
            a(6) = HO2H
     rk  - reaction constants for the mechanism
            1) NO2 + HV = NO + O
            2) O + O2 = O3
            3) O3 + NO = NO2 + O2
            4) HCHO + HV = 2 HO2. + CO
            5) HCHO + HV = H2 + CO
            6) HCHO + HO. = HO2. + CO + H2O
            7) HO2. + NO = HO. + NO2
            8) HO. + NO2 = HNO3
            9) HO2H + HV = 2 HO.
           10) HO2H + HO. = H2O + HO2.

    Outputs:
     s   - steady state species concentrations
            s(1) = HO.
            s(2) = O


     r   - reaction rates for mechanism listed above
     fr  - formation rates for each of the active species (ppm min-1)
     lr  - loss rates for each of the active species (ppm min-1)
    **************************************************************************=#

    # DEFINE RATES OF REACTIONS, PARTIAL RATES
    # FOR REACTIONS WITH STEAdy STATE REACTANTS


     r[1] = rk[1]*a[3]*c[2]
     r[2] = rk[2]*c[1]
     r[3] = rk[3]*a[1]*a[2]
     r[4] = rk[4]*a[4]*c[2]
     r[5] = rk[5]*a[4]*c[2]
     r[6] = rk[6]*a[4]
     r[7] = rk[7]*a[5]*a[2]
     r[8] = rk[8]*a[3]
     r[9] = rk[9]*a[6]*c[2]
     r[10] = rk[10]*a[6]

     # CALCULATE STEADY STATE CONCENTRATIONS AND
     # RATES OF REACTION WITH SS SPECIES

     # HO.
     s[1] = (r[7]+r[9]*2.0)/(1.0e-20+r[6]+r[8]+r[10])
     r[6] = s[1]*r[6]
     r[8] = s[1]*r[8]
     r[10] = s[1]*r[10]

     # O
     s[2] = r[1]/(1e-20+r[2])
     r[2] = s[2]*r[2]

     # DEFINE FORMATION RATES FOR ACTIVE SPECIES

     # O3
     fr[1] = r[2]

     # NO
     fr[2] = r[1]

     # NO2
     fr[3] = r[3] + r[7]

     # HCHO
     fr[4] = 0.0

     # HO2.
     fr[5] = 2*r[4] + r[6] + r[10]

     # HO2H
     fr[6] = 0.0

     # DEFINE LOSS RATES OF ACTIVE species

     # O3
     lr[1] = rk[3]*a[2]

     # NO
     lr[2] = rk[3]*a[1] + rk[7]*a[5]

     # NO2
     lr[3] = rk[1]*c[2] + rk[8]*s[1]

     # HCHO
     lr[4] =  rk[4]*c[2]+rk[5]*c[2]+rk[6]*s[1]

     # HO2.
     lr[5] = rk[7]*a[2]

     # HO2H
     lr[6] =  rk[9]*c[2]+rk[10]*s[1]
     return s,r,fr,lr
 end
