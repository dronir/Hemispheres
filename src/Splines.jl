module Splines

export Spline, SplineFunction

immutable Spline
    knots::Vector{Float64}
    yval::Vector{Float64}
    deriv::Vector{Float64}
	function Spline(xval::Vector{Float64}, yval::Vector{Float64}, deriv::Vector{Float64})
	    N = size(xval, 1)
	    
	    A = zeros(Float64, N)
	    B = ones(Float64, N)
	    C = zeros(Float64, N)
	    R = zeros(Float64, N)
	    gamma = zeros(Float64, N)
	    U = zeros(Float64, N)
	
	    R[1] = deriv[1]
	    for i = 2:N-1
	        A[i] = 1 / (xval[i] - xval[i-1])
	        B[i] = 2/(xval[i] - xval[i-1]) + 2/(xval[i+1]-xval[i])
	        C[i] = 1/(xval[i+1] - xval[i])
	        R[i] += 3 * (yval[i]-yval[i-1])/(xval[i]-xval[i-1])^2
	        R[i] += 3 * (yval[i+1]-yval[i])/(xval[i+1]-xval[i])^2
	    end
	    R[N] = deriv[2]
	    
	    # tri-diagonal solver
	    beta = B[1]
	    U[1] = R[1]/beta
	    for i = 2:N
	        gamma[i] = C[i-1] / beta
	        beta = B[i] - A[i] * gamma[i]
	        U[i] = (R[i] - A[i]*U[i-1]) / beta
	    end
	
	    for i = N-1 : -1 : 1
	        U[i] = U[i] - gamma[i+1]*U[i+1]
	    end
	    
	    S = new(copy(xval), copy(yval), copy(U))
	    return S
	end
end # type

function SplineFunction(spline::Spline, x::Real)
    n_knots = size(spline.knots,1)
    i = 0
    for j = 1:n_knots-1
        if spline.knots[j] <= x <= spline.knots[j+1]
            i = j
            break
        end
    end
    if i == 0
        println("Warning: x value outside spline!")
        return NaN
    end
    x1 = spline.knots[i]
    x2 = spline.knots[i+1]
    y1 = spline.yval[i]
    y2 = spline.yval[i+1]
    d1 = spline.deriv[i]
    d2 = spline.deriv[i+1]
    xd = x2-x1
    yd = y2-y1
    a = d1*xd - yd
    b = -d2*xd + yd
    t = (x-x1)/xd
    y = (1-t)*y1 + t*y2 + t*(1-t)*(a*(1-t)+b*t)
    return y
end

end # module
