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

    tra1=ode45(der,[et0;etm1],x0)
    t0=tra1[1][end]
    xf=tra1[2][end,:]
    x0=[xf[1:3];xf[4:6]+dv1]
    tra2=ode45(der,[t0;et1],x0)

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

function frrrate(body::ASCIIString,dv1::Vector=[0;0;0], dv2::Vector=[0;0;0],
               dv3::Vector=[0;0;0])

    tra=frtraj(dv1,dv2,dv3)
    mjd2k_TDB=tra[1]/86400+0.5
    
    #obtain moon states for range of epochs 
    xm=jpleph(mjd2k_TDB,body,"Earth");
    #relative position and velocity vectors are needed
    deltaR=traj[2][:,1:3]-xm[1:3,:]';
    deltaV=traj[2][:,4:6]-xm[4:6,:]';

    # rr_dot=(vec(r),vec(r_dot))
    r=sqrt(diag(deltaR*deltaR'))
    # the "diag" term is the inner product
    rrate=diag(deltaR*deltaV')./r

    return(rrate)
end

function frinc(dv1=[0;0;0],dv2=[0;0;0],
               dv3=[0;0;0])
    el=elements(frxf(dv1,dv2,dv2),planets["earth"]["mu"]);
    el[3]
end
