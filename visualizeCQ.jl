function visualizeCQ(C,V;title=" ",legend=true)
    C = C*1000 # ppm to ppb 
    # L1 normalization, so units aren't distorted (so we can still plot in ppb)  
    CQ3 = V[3,:]/sum(V[3,:])

    time = 0:1e-2:20


    p = plot(time,C[:,1],label=s"$O_3$",xlabel="Time [minutes]",title=title,xguidefontsize=15)
    p = plot!(time,C[:,2],label=s"$NO$",ylabel="Mixing ratio [ppb]",yguidefontsize=15)
    p = plot!(time,C[:,3],label=s"$NO_2$")

    if legend==true 
        plot!(time,C*CQ3,label=s"$CQ_3$",legend=:topright,linewidth=2)
    else
        plot!(time,C*CQ3,label=s"$CQ_3$",legend=legend,linewidth=2)
    end
end