module PhaseFunctions

using Abstracts
using ScatteringGeometry

export LSintegral, LambertSphere, NumericalPhaseFunction, Isotropic

# ---- Isotropic phase function ----
immutable Isotropic <: PhaseFunction
end
value(S::Isotropic, G) = 1.0

# ---- Lommel-Seeliger sphere integrated brightness ----

immutable LSintegral <: PhaseFunction
	omega::Float64
	P::PhaseFunction
end

value(S::LSintegral, alpha::Real) = (1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4))) / 32
value(S::LSintegral, G::Geometry) = value(S, phase_angle(G))
value(::Type{LSintegral}, G::Geometry) = value(LSintegral(), G)
value(::Type{LSintegral}, alpha::Real) = value(LSintegral(), alpha)


# ---- Lambert sphere integrated brightness ----


immutable LambertSphere <: PhaseFunction
end

value(::Type{LambertSphere}, G::Geometry) = value(LambertSphere(), G)
value(S::LambertSphere, alpha::Real) = (sin(alpha) + (pi - alpha) * cos(alpha)) / (6 * pi)
value(S::LambertSphere, G::Geometry) = value(S, phase_angle(G))


# ---- Numerically interpolated phase function ----

immutable NumericalPhaseFunction <: ScatteringLaw
	data::Array
end

NumericalPhaseFunction(filename::String) = NumericalPhaseFunction(readcsv(filename))

function value(S::NumericalPhaseFunction, G::Geometry)
	alpha = phase_angle(G)
	return interpolate_phasecurve(S.data, alpha)
end

function interpolate_phasecurve(Curve::Array, alpha::Real)
    for i = 1:size(Curve)[1]-1
        if Curve[i+1, 1] > alpha
            delta = (alpha - Curve[i, 1]) / (Curve[i+1, 1] - Curve[i, 1])
            P = Curve[i, 2] + (Curve[i+1, 2] - Curve[i, 2]) * delta
            return P
		end
	end
end


	
end # module
