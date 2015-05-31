module Astrodynamics

    using Datetime
    using jspice

    importall Base

    export State
    export julian, gregorian, JD2000, JD1950, MJD, ephtim
    export ecctomean, meantoecc, ecctotrue, truetoecc
    export propagate
    export planets
    export cartesian, elements
    export newton,findiff,fingrd
    export eci2emrotpulse,emrotpulse2eci
    export getemrpstate

    immutable State
        rv::Vector
        t::DateTime
        frame::String
        body::String

        function State(rv0::Vector,t0::FloatingPoint,f::String,b::String)
            calt0=gregorian(t0)
            dt=datetime(calt0[1],calt0[2],calt0[3],calt0[4],calt0[5],calt0[6])
            if f=="emrp"
                rv=emrotpulse2eci(rv0,t0);
            elseif f=="eci"
                rv=rv0;
            else
                println("invalid frame designator");
            end
            new(rv,dt,"eci","earth");
        end
        function State(rv0::Vector,t0::DateTime,f::String,b::String)
            dt=t0
            if f=="emrp"
                rv=emrotpulse2eci(rv0,julian(t0,"2000"))
            elseif f=="eci"
                rv=rv0
            else
                println("invalid frame designator")
            end
            new(rv,dt,"eci","earth")
        end
    end
    
    include("constants.jl")
    include("time.jl")
    include("frames.jl")
    include("math.jl")
    include("kepler.jl")

    iss = State([8.59072560e+02; -4.13720368e+03; 5.29556871e+03; 7.37289205e+00; 2.08223573e+00; 4.39999794e-01],
        julian(2013,3,18,12,0,0,"2000"), "eci", "earth")

    function propagate(s::State; method::String="kepler")
       mu = planets[s.body]["mu"] 
       if method == "kepler"
           ele = elements(s.rv, mu)
           tend = period(ele[1], mu)
           dt = [0:tend]
           return cartesian(kepler(ele, dt, mu), mu)
       end
    end
end
