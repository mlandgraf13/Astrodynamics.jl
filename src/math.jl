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