#---------------------------------------------------
# prop is a simple propagator to final time or
# to next pericentre
#---------------------------------------------------

function prop(x0::Vector,et0::FloatingPoint,etf::FloatingPoint;
              refbod="earth",stopco=falses(1,2))

    mu=planets[refbod]["mu"]
    tra_t,tra_x=ode45(der,x0,[et0;etf])
    mjd2k=tra_t/86400+0.5
    et=tra_t
    tra=(tra_t,tra_x)

    if stopco[1] # pericentre stop condition
        n=length(tra_t)
        if (refbod=="earth")
            delta=[tra_x[j][:] for j=1:n]
        else
            xm=jpleph(mjd2k,refbod,"earth")
            delta=[tra_x[j][:]-xm[:,j] for j=1:n]
        end
        
        deltar=zeros(3,n)
        for j=1:n
            deltar[:,j]=delta[j][1:3]
        end
        r=sqrt(diag(deltar'*deltar))
        dr=(r[2:end]-r[1:end-1])./(et[2:end]-et[1:end-1])
        # detect all sign changes from negative to positive in dr
        sgnchg=find((sign(dr[2:end]).*sign(dr[1:end-1]) .< 0.0) 
                    & (sign(dr[1:end-1]).<0.0))
        if (~isempty(sgnchg))
            #use the first change of sign
            idx=sgnchg[1]
            x0=delta[idx]
            t0=tra_t[idx]
            el0=elements(x0,mu)
            f0=rem(el0[6],2pi)-2pi
            M0=ecctomean(truetoecc(f0,el0[2]),el0[2])

            f=round(f0/pi)*pi
            el=[el0[1:5];f]

            M=ecctomean(truetoecc(f,el[2]),el[2])
            n=el[1] > 0.0 ? sqrt(mu/el[1]^3) : sqrt(mu/(-el[1]^3))
            dt=(M-M0)/n

            t=t0+dt

            if (refbod=="earth")
                x=cartesian(el,mu)
            else
                mjd2k=t/86400+0.5
                xm=jpleph(mjd2k,refbod,"earth")
                x=cartesian(el,mu) + xm
            end

            # the trajectory is truncated at the pericentre
            newt=[tra_t[1:idx];t]
            newx=Array(Any,idx+1)
            newx[1:idx]=tra_x[1:idx]
            newx[idx+1]=x

            tra=(newt,newx)
        end
    elseif stopco[2] # xz crossing
        n=length(tra_t)
        rottra=[eci2emrotpulse(tra_x[j],mjd2k[j]) for j=1:n]
        # detect all sign changes of yr
        yr=[rottra[j][2] for j=1:n]
        syrright=sign(yr[2:end])
        syrleft=sign(yr[1:end-1])
        sgnchg=find((syrright.*syrleft .< 0.0))
        if ~isempty(sgnchg)
            #use the first change of sign
            idx=sgnchg[1]
            ## bisection to determine x-z crossing
            # initialise bisection
            et0=tra_t[idx]
            et1=tra_t[idx+1]
            x0=tra_x[idx]
            x1=tra_x[idx+1]
            Dt=et1-et0
            const bisect_eps=1.0 #s
            while abs(Dt) > bisect_eps
                # new point: in the middle of the interval
                et=.5(et0+et1)
                #propagate to new point - start bisection loop here
                arc=ode45(der,x0,[et0;et])
                x=arc[2][end]
                mjd2k=et/86400+0.5
                xrot=eci2emrotpulse(arc[2][end],mjd2k);
                # re-set interval limits based on sign change or not of yrot
                if sign(xrot[2])*yr[idx]<0.0
                    # sign change took place
                    Dt=et1-et
                    et1=et
                    x1=x
                else
                    # sign change did not take place
                    Dt=et-et0
                    et0=et
                    x0=x
                end
            end
            # concatenate new trajectory with x-z crossing as new final point
            newt=[tra_t[1:idx];et]
            newx=Array(Any,idx+1)
            newx[1:idx]=tra_x[1:idx]
            newx[idx+1]=x
            tra=(newt,newx)
        end
    end
    return(tra)
end
