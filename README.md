# Julia photochemical model
The Julia photochemical model (formerly FormO3NOx) is small photochemical mechanism with reactions related to ozone formation: NOx chemistry, VOC chemistry, formation of a peroxy radical from VOC chemistry, and further reaction of the peroxy radical with NO to form NO2 and OH. It is based on a Fortran mechanism from Dr. Michael Kleeman's graduate course ECI241 on air quality modeling at University of California, Davis.

# Code Structure and Output
The top level module for the mechanism is contained in the file assign3_driver.jl.  The wrapper program main.jl runs this module for a user-specified number of cases of user-specified duration, and generating concentration data at specified intervals. Further documentation is contained in the .jl files.

The "jpm-h2o" branch adds net production of H2O as a buildup species.

# Reactions and Species

## Reactions
1)  NO2 + hv -> NO + O

2)  O + O2 -> O3

3)  O3 + NO -> NO2 + O2

4)  HCHO + hv -> 2 HO2. + CO

5)  HCHO + hv -> H2 + CO

6)  HCHO + HO. -> HO2. + CO + H2O

7)  HO2. + NO -> OH. + NO2

8)  OH. + NO2 -> HNO3

9)  HO2H + hv -> 2 HO.

10) HO2H + HO. -> H2O + HO2.

## Species
1) O3

2) NO

3) NO2

4) HCHO

5) HO2.

6) HO2H

7) OH.

8) O

9) HNO3

10) CO

11) H2

12) H2O

# Other Information
This mechanism was written in Julia Version 1.2.0 on a macOS (x86_64-apple-darwin18.6.0) and last tested with Julia v1.8.1. Full documentation is contained within the .jl files.

An archived version of this code can be found on Zenodo:
https://doi.org/10.5281/zenodo.3733503

Feel free to direct all questions to Obin Sturm, posturm@ucdavis.edu
