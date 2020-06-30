#=*****************************************************************************
 This is a module in an ozone NOx formaldehyde mechanism, adapted from Fortran
 by Obin Sturm (UC Davis) in October 2019.

 This module was recieved from Mike Kleeman (UC Davis) in October 2018.
 ****************************************************************************=#


#  parameter definitions common to several subroutines

const maxcon, maxact, maxrxn, maxbo, maxsts = 10, 8, 20, 3, 2
const constn = zeros(Float64,1,maxcon) # note: renamed from Fortran const due to conflict with Julia function const()
const rpssa  = zeros(Float64,1,maxsts)
const r  = zeros(Float64,1,maxrxn)
const rk = zeros(Float64,1,maxrxn)
const fr  = zeros(Float64,1,maxact)
const rlr  = zeros(Float64,1,maxact)
