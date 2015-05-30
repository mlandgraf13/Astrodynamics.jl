# transform state vector from Earth-centric inertial MEE2000 to rotating 
# pulsating frame
function eci2rotpulse(rv::Vector, t::FloatingPoint)
    sM=jpleph(t,"Moon","Earth");
    rM=norm(sM[1:3]);
    vM=norm(sM[4:6]);
# unit vector in Earth Moon direction is x-base vector
    ex=sM[1:3]/rM;
    ez=cross(ex, sM[4:6]/vM);
    ey=cross(ez, ex);
# rotation matrix equals matrix with base vectors as columns
    R=[ex ey ez];
    posrot=R'*rv[1:3]
    velrot=R'*rv[4:6]

# velocity compoment of rotating frame
    rotvec=cross(sM[1:3],sM[4:6]);
    omega=norm(rotvec)/rM^2;
    velrot[1]=velrot[1]-(-omega*posrot[2]);
    velrot[2]=velrot[2]-omega*posrot[1];
    
# position and velocity in natural units for scaling 
#  (i.e. Moon always at [1;0;0])
    npos=posrot/rM;
    nvel=velrot/(omega*rM);
    return [npos;nvel];
end

function elements(rv::Vector, mu::FloatingPoint)
    r, v = rv[1:3], rv[4:6]
    rm = norm(r)
    vm = norm(v)
    h = cross(r, v)
    hm = norm(h)
    k = [0, 0, 1]
    n = cross(k, h)
    nm = norm(n)
    xi = vm^2/2 - mu/rm
    ec = ((vm^2 - mu/rm)*r - v*dot(r, v))/mu
    ecc = norm(ec)
    if ecc != 1
        sma = -mu/(2*xi)
        p = sma*(1-ecc^2)
    else
        p = hm^2/mu
        sma = p
    end
    inc = acos(h[3]/hm)
    node = acos(n[1]/nm)
    peri = acos(dot(n, ec)/(ecc*nm))
    # this is a workaround to avoid domainerror - dot and/or
    #  norm must be fixed
    ano = acos(max(min(dot(ec, r)/(ecc*rm),1),-1))
    if n[2] < 0
        node = 2*pi - node
    end
    if ec[3] < 0
        peri = 2*pi - peri
    end
    if dot(r, v) < 0
        ano = 2*pi - ano
    end
    return [sma, ecc, inc, node, peri, ano]
end

function elements(rv::Matrix, mu::FloatingPoint)
    m, n = size(rv)
    if m != 6 && n != 6
        error("'rv' must be a 6xN or Nx6 matrix.")
    elseif m == 6
        rv = rv'
        m, n = n, m
    end
    ele = zeros(m,n)
    for i = 1:m
        ele[i,:] = elements(vec(rv[i,:]), mu)
    end
    return ele
end

function elements(s::State)
    return elements(s.rv, planets[s.body]["mu"])
end

function elements(s::State, deg::Bool)
    ele = elements(s)
    if deg
        ele[3:end] = ele[3:end]*180/pi
        return ele
    else
        return ele
    end
end

function cartesian(ele::Vector, mu::FloatingPoint)
    sma, ecc, inc, lan, per, ano = ele
    u = per + ano
    if ecc == 1
        p = sma
    else
        p = sma*(1 - ecc^2)
    end
    r = p/(1 + ecc*cos(ano))
    x = r*(cos(lan)*cos(u) - sin(lan)*cos(inc)*sin(u))
    y = r*(sin(lan)*cos(u) + cos(lan)*cos(inc)*sin(u))
    z = r*sin(inc)*sin(u)
    vr = sqrt(mu/p)*ecc*sin(ano)
    vf = sqrt(mu*p)/r
    vx = ((vr*(cos(lan)*cos(u) - sin(lan)*cos(inc)*sin(u))
        - vf*(cos(lan)*sin(u) + sin(lan)*cos(u)*cos(inc))))
    vy = ((vr*(sin(lan)*cos(u) + cos(lan)*cos(inc)*sin(u))
        - vf*(sin(lan)*sin(u) - cos(lan)*cos(u)*cos(inc))))
    vz = vr*sin(inc)*sin(u) + vf*cos(u)*sin(inc)
    return [x, y, z, vx, vy, vz]
end

function cartesian(ele::Matrix, mu::FloatingPoint)
    m, n = size(ele)
    if m != 6 && n != 6
        error("'ele' must be a 6xN or Nx6 matrix.")
    elseif m == 6
        ele = ele'
        m, n = n, m
    end
    rv = zeros(m,n)
    for i = 1:m
        rv[i,:] = cartesian(vec(ele[i,:]), mu)
    end
    return rv
end
