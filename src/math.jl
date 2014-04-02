function newton(x0::Float64, func::Function, derivative::Function, 
                maxiter::Int=50, tol::Float64=sqrt(eps()))
    p0 = x0
    for i = 1:maxiter
        i += 1
        p = p0 - func(p0)/derivative(p0)
        if abs(p - p0) < tol
            return p
        end
        p0 = p
    end
    error("Not converged.")
end

#---------------------------------------------------
# newton method with no derivative specified using finite differences
# Written by Markus Landgraf. Last modified 20Mar14
#---------------------------------------------------
function newton(x0::Float64, func::Function, maxiter::Int=50, 
                tol::Float64=sqrt(eps()))

    deriv(x)=findiff(func,x);
    newton(x0,func,deriv,maxiter,tol)

end

#---------------------------------------------------
# finite difference for numerical differentiation
# written by Markus Landgraf. Last modified 20Mar14
#---------------------------------------------------
function findiff(F::Function,x0::Float64,tol::Float64=(eps())^(1/3))
    dx = (x0==0 ? tol:x0*tol)
    return( (F(x0+dx) - F(x0-dx)) / (2*dx) )
end

#---------------------------------------------------
# finite difference differentiation calculation of gradient
# written by Markus Landgraf. Last modified 02Apr14
# CHANGE RECORD
# 30Mar14: ML - initial coding
# 02Apr14: ML - adoption for vector functions
#---------------------------------------------------
function fingrd(F::Function,x0::Vector,tol::Float64=(eps())^(1/3))
    n=length(x0)
    F0=F(x0)
    dF=zeros(length(F0),n)

    for j=1:n
        F1(x)=F([x0[1:j-1];x;x0[j+1:end]])
        dF[:,j]=findiff(F1,x0[j],tol)
    end

    return(dF)
end

function finvgrd(F::Array{Function,1},x0::Vector,tol::Float64=(eps())^(1/3))
    m=length(F)
    n=length(x0)
    
    vgrd=zeros(m,n)
    for j=1:m
        vgrd[j,:]=fingrd(F[j],x0,tol)
    end
    
    return(vgrd)
end