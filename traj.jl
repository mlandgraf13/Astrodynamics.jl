function frtraj(dv1::Vector=[0;0;0],dv2::Vector=[0;0;0],
               dv3::Vector=[0;0;0])

    const x0=[4179.47;
              -4455.25;
              -2399.48;  
              6.63923;
              8.03688;
              -3.35817;
              ]
    const et0=5.68022588928017e8
    const et1=5.68474771104017e8+2*86400
    const etm1=et0+86400
    const dtm3=86400

    # first arc from perigee to manoeuvre 1
    tra1=prop(x0,et0,etm1)
    t0=tra1[1][end]
    xf=tra1[2][end,:]
    # add delta-v1
    x0=[xf[1:3];xf[4:6]+dv1]

    #second arc from manoeuvre 1 to periselenium manoeuvre
    tra2=prop(x0,t0,et1,"moon",true);
    t0=tra2[1][end]
    xf=tra2[2][end,:]
    x0=[xf[1:3];xf[4:6]+dv2]
    
    #third arc from periselenium manoeuvre to manoeuvre 3
    tra3=prop(x0,t0,t0+dtm3)
    t0=tra1[1][end]
    xf=tra1[2][end,:]
    # add delta-v3
    x0=[xf[1:3];xf[4:6]+dv3]

    #fourth arc from manoeuvre 3 to perigee
    tra4=prop(x0,t0,et1,"earth",true)

    return([tra1[1];tra2[1];tra3[1]],
           [tra1[2];tra2[2];tra3[2]])
end

function frxf(dv1::Vector=[0;0;0], dv2::Vector=[0;0;0],
              dv3::Vector=[0;0;0])

   tra=frtraj(dv1,dv2,dv3)

   return(tra[2][end,:])
end

function frdist(body::ASCIIString,dv1::Vector=[0;0;0], dv2::Vector=[0;0;0],
               dv3::Vector=[0;0;0])

    tra=frtraj(dv1,dv2,dv3)
    mjd2k_TDB=tra[1]/86400+0.5
    
    xm=jpleph(mjd2k_TDB,body,"Earth");
    delta=traj[2][:,1:3]-xm[1:3,:]';
    r=sqrt(diag(delta*delta'))
end

function frinc(dv1=[0;0;0],dv2=[0;0;0],
               dv3=[0;0;0])
    el=elements(frxf(dv1,dv2,dv2),planets["earth"]["mu"]);
    el[3]
end

#---------------------------------------------------
# prop is a simple propagator to final time or
# to next peri/apocentre
#---------------------------------------------------
function prop(x0,et0,etf,refbod="earth",stopco=falses(1))

    mu=planets[refbod]["mu"]
    tra=ode45(der,[et0;etf],x0)
    mjd2k_TDB=tra[1]/86400+0.5
    et=tra[1]
    
    if(stopco[1]) # peri/apocentre stop condition
        if (refbod=="earth")
            delta=tra[2]
        else
            xm=jpleph(mjd2k_TDB,refbod,"earth")
            delta=tra[2]-xm'
        end

        deltar=delta[:,1:3]
        r=sqrt(diag(deltar*deltar'))
        dr=(r[2:end]-r[1:end-1])./(et[2:end]-et[1:end-1])
        sgnchg=find(sign(dr[2:end]).*sign(dr[1:end-1]).<0.0)
        idx=sgnchg[1]
        #debug
        display(r[idx:idx+3])
        display(dr[idx:idx+3])
        #use the first change of sign
        x0=delta[idx,:]
        t0=tra[1][idx]
        el0=elements(x0,mu)
        f0=rem(el0[6],2pi)
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
            mjd2k_TDB=t/86400+0.5
            xm=jpleph(mjd2k_TDB,refbod,"earth")
            x=cartesian(el,mu) + xm
        end
        newt=[tra[1][1:idx];t]
        newx=[tra[2][1:idx,:];x']

#debug
        display(f0)
        display(f)
        display(M0)
        display(M)
        display(dt)
        
        tra=(newt,newx)
    end
    return(tra)
end


