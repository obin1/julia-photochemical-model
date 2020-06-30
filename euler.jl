function euler!(t1,t2,dt,tm,c,daymax,n::Int64,rk,dc::Function)

#=**************************************************************************
This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
by Obin Sturm (UC Davis) in October 2019.
The original comments and file description are below.

    created by: Mike Kleeman (April 1999)
                ECI 298 Air Quality Modeling
                University of California, Davis

    The purpose of this subroutine is to use an explicit Euler method to
    integrate a set of coupled ordinary differential equations between t1
    and t2.

    Inputs:
     t1   - start time
     t2   - stop time
     dt   - fixed integration time step
     tm   - concentration recording timestep
     c    - initial active species concentrations
     n    - number of active species
     diff - name of the subroutine to evaluate the derivative dc/dt = diff

    Outputs:
     C    -  species concentrations reported at specific time intervals
**************************************************************************=#
    #include("modlspc.jl")
    f = zeros(1,n)
    C = zeros(ceil(Int,t2/tm)+1,n)
    C[1,:] = copy(c[1,:])
    S = zeros(ceil(Int,t2/tm),length(r))
    s_aggregate = zeros(Float64,1,maxrxn)
    J = zeros(ceil(Int,t2/tm)+1,1)
    J[1,1] = copy(constn[2])

    #return C
    #c = zeros(1,n)
# preliminary steps
    if n > maxrxn
        error("too many reactions")
    end
    Time = t1
# loop over the time step until we come to t2
    icount=Int64(1)


    while (Time < t2)
        # oscillating HV value during the day
        constn[2] = daymax*max(sin(Time/720*pi - pi/2),0.0)

        f = dc(c,f,n)
        #       print*,'report from euler:c,f'
        #       do i=1,n
        #        print*,i,c(i),f(i)
        #       enddo
        #       stop 'euler marker1'

        c .= c .+ f*dt
        if any(c.<0.0)
            c[c.<0] .= 0
            #s_aggregate[1] = NaN  # Added March 2020: marker to show possible mass balance violation
            #println("Mass balance prob")
        end


        Time = Time + dt
        s_aggregate = s_aggregate .+ r.*dt

        if mod(icount,floor(Int,1000*tm)) == 0
            C[icount ÷ floor(Int,1000*tm)+1,:] = copy(c[1,:])
            S[icount ÷ floor(Int,1000*tm),:]   = copy(s_aggregate)
            J[icount ÷ floor(Int,1000*tm)+1,1] = copy(constn[2])
            s_aggregate = zeros(Float64,1,maxrxn)
        end
        icount=icount+1

    end

    return C,S[:,1:10],J
end


function diff!(c,f,n)

#=**************************************************************************
    This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
    by Obin Sturm (UC Davis) in October 2019.
    The original comments and file description are below.
    created by: Mike Kleeman (April 1999)
             ECI 298 Air Quality Modeling
             University of California, Davis

 The purpose of this subroutine is to evaluate the right hand side of
 the equation
                    dc/dt = f
 Inputs:
  c    - active species concentrations
  n    - number of active species

 Outputs:
  f    - net rate of change of c
**************************************************************************=#
    #include("modlspc.jl")

# call the subroutine to evaluate the rates for the chemical reaction sys
    global s, r, fr, rlr = difun!(constn,c,rpssa,rk,r,fr,rlr)
    f[maxact+1:(maxbo+maxact)] = bldup!(r,f[maxact+1:(maxact+maxbo)])

# calculate the net rate of reaction for active species
    f[1:maxact] .= fr[1:maxact] .- rlr[1:maxact].*c[1:maxact]
    return f
end

function jcalc(c,n,dc,rj)
#=**************************************************************************
This is a function in an ozone NOx formaldehyde mechanism, adapted from Fortran
by Obin Sturm (UC Davis) in October 2019.

Note that this subroutine was never called in the original mechanism.

The original comments and file description are below.

    created by: Mike Kleeman (April 1999)
                ECI 298 Air Quality Modeling
                University of California, Davis

   The purpose of this subroutine is to numerically evaluate the Jacobian
   of the right hand side of the equation
                      dc/dt = f

   Inputs:
    c   - active species concentration
    n   - number of active species
    dc  - name for the function f

   Outputs:
    rj  - jacobian matrix for function f
**************************************************************************=#
end
