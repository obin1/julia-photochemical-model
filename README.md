# FormO3NOx
A small photochemical mechanism with reactions related to ozone formation: NOx chemistry, VOC chemistry, formation of a peroxy radical from VOC chemistry, and further reaction of the peroxy radical with NO to form NO2 and OH. This Julia module is based on a Fortran mechanism from Dr. Michael Kleeman's graduate course at University of California, Davis on atmospheric chemistry modeling.

# Code Structure and Output
The top level module for the mechanism is contained in the file assign3_driver.jl.  The wrapper program MLAQ.jl runs this module for 365 simulated days, generating concentration data every 6 minutes, saved as the text file C_and_J.txt. The rows of C.txt correspond to the time steps, one column to a daylight constant "J" between zero and one which affects the rate of the photochemical reactions, and the rest of the columns to the species. This mechanism also generates the integrated rate of reaction over every 6 minutes and saves this information as S.txt, rows corresponding to new time-steps and columns corresponding to the reactions.

The output files C_and_J.txt and S.txt are available for download.

Further documentation is contained in the .jl files.
