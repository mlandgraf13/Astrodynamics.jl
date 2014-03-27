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
    const et1=5.68474771104017e8
    const etm1=et0+86400

    tra1=prop(x0,et0,etm1)
    t0=tra1[1][end]
    xf=tra1[2][end,:]
    x0=[xf[1:3];xf[4:6]+dv1]
    tra2=prop(x0,t0,et1)

    return([tra1[1];tra2[1]],[tra1[2];tra2[2]])
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
        if (refbod=="Earth")
            delta=traj[2]
        else
            xm=jpleph(mjd2k_TDB,refbod,"earth")
            delta=tra[2]-xm'
        end
        deltar=delta[:,1:3]
        r=sqrt(diag(deltar*deltar'))
        dr=(r[2:end]-r[1:end-1])./(et[2:end]-et[1:end-1])
        sgnchg=find(sign(dr[2:end]).*sign(dr[1:end-1]).<0.0)
        #use the first change of sign
        x0=delta[sgnchg[1],:]
        t0=tra[1][sgnchg[1]]
        el0=elements(x0,mu)
        f0=el0[6]
        f=round(f0/pi)*pi
        el=[el0[1:5];f]
        x=cartesian(el,mu)
        M0=ecctomean(truetoecc(f0,el[2]),el[2])
        M=ecctomean(truetoecc(f,el[2]),el[2])
        n=sqrt(mu/el(1)^3)
        dt=(M-M0)/n
        t=t0+dt
        tra=([tra1[1][1:sgnchg];t],[tra[2][1:sgnchg,:];x])
    end
    return(tra)
end


