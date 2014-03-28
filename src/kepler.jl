function meantoecc(M::FloatingPoint, ecc::FloatingPoint)
    kepler(E) = ecc < 1.0 ? E - ecc*sin(E) - M : ecc*sinh(E) - E
    kepler_der(E) = ecc < 1.0 ? 1 - ecc*cos(E) : ecc*cosh(E) - 1
    return newton(M, kepler, kepler_der)
end

function ecctomean(E::FloatingPoint, ecc::FloatingPoint)
    return ecc<1.0 ? E - ecc*sin(E) : ecc*sinh(E)-E
end

function ecctotrue(E::FloatingPoint, ecc::FloatingPoint)
    return ecc < 1.0 ? 
           2*atan2(sqrt(1 + ecc)*sin(E/2), sqrt(1 - ecc)*cos(E/2)) :
           2*atan2(sqrt(ecc + 1)*sin(E/2), sqrt(ecc - 1)*sin(E/2))
end

function truetoecc(T::FloatingPoint, ecc::FloatingPoint)
    return ecc<1.0 ?  
        2*atan2(sqrt(1 - ecc)*sin(T/2), sqrt(1 + ecc)*cos(T/2)) :
        2*atan2(sqrt(ecc - 1)*sin(T/2), sqrt(ecc + 1)*cos(T/2))
end

function period(a::FloatingPoint, mu::FloatingPoint)
    return sqrt(4*a^3*pi^2/mu)
end

function period(s::State)
    mu = planets[s.body]["mu"]
    ele = elements(s)
    return period(ele[1], mu)
end

function kepler(ele::Vector, dt::FloatingPoint, mu::FloatingPoint)
    E0 = truetoecc(ele[6], ele[2])
    M0 = ecctomean(E0, ele[2])
    n = 2*pi/period(ele[1], mu)
    M = M0 + n*dt
    E = meantoecc(M, ele[2])
    T = ecctotrue(E, ele[2])
    return [ele[1:5], T]
end

function kepler(ele::Vector, dt::Vector, mu::FloatingPoint)
    out = Array(FloatingPoint, length(dt), 6)
    for i = 1:length(dt)
        out[i,:] = kepler(ele, dt[i], mu)
    end
    return out
end
