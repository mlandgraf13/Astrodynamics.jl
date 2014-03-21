using Astrodynamics
using DerIMDEX
using ODE
using jspice

function traj()
    const x0=[4179.47;
                          -4455.25;
                          -2399.48;  
                          6.63923;
                          8.03688;
                          -3.35817;
                          ]
    const et0=5.68022588928017e8
    const et1=5.68474771104017e8

    tra=ode45(der,[et0;et1],x0)
    tra[2][end,:]
end

function frinc()
    el=elements(traj(),planets["earth"]["mu"]);
    el[3]
end