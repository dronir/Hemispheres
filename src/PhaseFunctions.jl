module PhaseFunctions

using Abstracts
using ScatteringGeometry
using Splines

export LSintegral, LambertSphere, Numerical, Isotropic

# ---- Generic value function from Geometry ----
value(S::PhaseFunction, G::Geometry) = value(S, phase_angle(G))


# ---- Isotropic phase function ----
immutable Isotropic <: PhaseFunction
end
value(S::Isotropic, G::Geometry) = 1.0
value(S::Isotropic, alpha::Real) = 1.0


# ---- Lommel-Seeliger sphere integrated brightness ----

immutable LSintegral <: PhaseFunction
end

LSi(alpha::Real) = alpha<eps() ? 1.0 : 1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4)) 

value(S::LSintegral, alpha::Real) = LSi(alpha) / 32
value(::Type{LSintegral}, alpha::Real) = value(LSintegral(), alpha)


# ---- Lambert sphere integrated brightness ----

immutable LambertSphere <: PhaseFunction
end

value(S::LambertSphere, alpha::Real) = (sin(alpha) + (pi - alpha) * cos(alpha)) / (6 * pi)


# ---- Numerically interpolated phase function ----

immutable Numerical <: PhaseFunction
	data::Array
end

Numerical(filename::String) = Numerical(readcsv(filename))

value(S::Numerical, alpha::Real) = interpolate_phasecurve(S.data, alpha)

function interpolate_phasecurve(Curve::Array, alpha::Real)
    for i = 1:size(Curve)[1]-1
        if Curve[i+1, 1] > alpha
            delta = (alpha - Curve[i, 1]) / (Curve[i+1, 1] - Curve[i, 1])
            P = Curve[i, 2] + (Curve[i+1, 2] - Curve[i, 2]) * delta
            return P
		end
	end
end


# --- Double Henyey-Greenstein phase function ----

immutable HG2 <: PhaseFunction
    w::Float64
    g1::Float64
    g2::Float64
end

value(S::HG2, alpha::Real) = S.w * HG(alpha, S.g1) + (1-S.w) * HG(alpha, S.g2)

HG(alpha, g) = (1-g^2) / (1 + g^2 + 2*g*cos(alpha))^1.5


# ---- H, G1, G2 phase function ----

immutable HG1G2 <: PhaseFunction
	H::Float64
	G1::Float64
	G2::Float64
end
value(S::HG1G2, alpha::Real) = P_HG(alpha, S.H, S.G1, S.G2)



# ---- HG divided by another phase function ----

immutable HGreduced <: PhaseFunction
	p::Float64
	G1::Float64
	G2::Float64
	P::PhaseFunction
end

function HGreduced(p::Real, G12::Real, P::PhaseFunction)
	HGreduced(p, getG12(G12)..., P)
end

value(S::HGreduced, alpha::Real) = S.p * P_HG(alpha, 0.0, S.G1, S.G2) / (4*value(S.P, alpha))



# --- H, G1, G2 function, not a PhaseFunction but used in them ----

function getG12(G12::Real)
	if G12 < 0.2
		return (0.7527*G12 + 0.06164, -0.9612*G12 + 0.6270)
	else
		return (0.9529*G12 + 0.02162, -0.6125*G12 + 0.5572)
	end
end

function P_HG(alpha::Real, H::Real, G1::Real, G2::Real)
	if alpha >= 2.6179938779914944
		return 0.0
	end
	if alpha < 0.1308996938995747
		a = 1 - 6*alpha / pi
		b = 1 - 9*alpha / (5*pi)
	else
		a = zeta1(alpha)
		b = zeta2(alpha)
	end
	c = alpha < 0.5235987755982988 ? zeta3(alpha) : 0.0
	return 10^(-0.4*H) * (G1*a + G2*b + (1-G1-G2)*c)	
end

const S1 = Spline([7.5, 30, 60, 90, 120, 150] * pi/180, 
		    [0.75, 0.33486, 0.134106, 0.0511048, 0.0214657, 0.0036397],
			[-1.90986, -0.0913286])
const S2 = Spline([7.5, 30, 60, 90, 120, 150] * pi/180, 
			[0.925, 0.628842, 0.317555, 0.127164, 0.0223739, 0.000165057],
			[-0.572958, -8.6573138e-8])
const S3 = Spline([0.0, 0.3, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0, 30.0]*pi/180,
			[1.,0.833812,0.577354,0.421448,0.231742,0.103482,0.0617335,0.016107,0.0],
			[-0.106301, 0.0])

zeta1(alpha::Real) = SplineFunction(S1, alpha)
zeta2(alpha::Real) = SplineFunction(S2, alpha)
zeta3(alpha::Real) = SplineFunction(S3, alpha)


	
end # module
