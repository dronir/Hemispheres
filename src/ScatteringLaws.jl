

module ScatteringLaws

#include("Abstracts.jl")
#include("Abstracts.jl")
#include("Abstracts.jl")
using Abstracts
using ScatteringGeometry
using PhaseFunctions
using Hemispheres

using Cubature

import PhaseFunctions.value

export value, Lambert, LommelSeeliger, ParticulateMedium, AntiShadow, AntiR
export spherealbedo, geometricalbedo, integrated


# ---- Lambert ----

immutable Lambert <: AnalyticalScatteringLaw
end

value(::Type{Lambert}, G::Geometry) = value(Lambert(), G)
value(S::Lambert, G::Geometry) = bool(G) ? 1.0 : 0.0


planaralbedo(S::Lambert, theta::Real) = 1.0
geometricalbedo(S::Lambert) = 2/3
spherealbedo(S::Lambert) = 2/3
integrated(S::Lambert, alpha::Real) = (sin(alpha) + (pi-alpha)*cos(alpha))/6


# ---- Lommel-Seeliger ----

immutable LommelSeeliger <: AnalyticalScatteringLaw
	P::PhaseFunction
	omega::Float64
end
LommelSeeliger() = LommelSeeliger(Isotropic(), 1.0)

value(::Type{LommelSeeliger}, G::Geometry) = value(LommelSeeliger(), G)
function value(S::LommelSeeliger, G::Geometry)
	return bool(G) ? S.omega * value(S.P, G) / (cos(G.theta_e)+cos(G.theta_i)) / (4pi) : 0.0
end

function planaralbedo(S::LommelSeeliger, theta::Real) 
	mu = cos(theta)
	mu<eps(theta) ? 0.5*S.omega : S.omega*(mu*log(mu) - mu*log(mu+1) + 1) / 2
end

geometricalbedo(S::LommelSeeliger) = S.omega / 8 * value(S.P, 0.0)
spherealbedo(S::LommelSeeliger) = 2/3 * (1 - log(2))
integrated(S::LommelSeeliger, alpha::Real) = value(S.P, alpha)/32 * (alpha < eps() ? 1.0 : 1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4)))



# ---- Particulate Medium ----

immutable ParticulateMedium <: ScatteringLaw
	hemi::Hemisphere
	P::PhaseFunction
end
value(S::ParticulateMedium, G::Geometry) = G ? value(S.hemi,G) * value(S.P, G) : 0.0


# ---- Anti-Shadow ----

# Dividing this out of a hemiScatter output hemisphere
# gives the shadowing correction S.

immutable AntiShadow <: AnalyticalScatteringLaw
    omega0::Float64
end
AntiShadow() = AntiShadow(0.1)
function value(S::AntiShadow, G::Geometry)
	alpha = phase_angle(G)
    mu0 = cos(G.theta_i)
    mu = cos(G.theta_e)
    omegaV = 2/3 * (1 - log(2)) * S.omega0
    PV = P_LS(alpha)
	return 0.25/pi * mu0 / (mu + mu0) * omegaV * PV
end

P_LS(alpha) = 0.75*(1 - sin(alpha/2) .* tan(alpha / 2) .* log(cot(alpha / 4))) / (1 - log(2))



# ---- Anti-R ----

# Dividing this out of a hemiScatter output hemisphere
# gives a reflection coefficient, sans phase function.

immutable AntiR <: AnalyticalScatteringLaw end
value(S::AntiR, G::Geometry) = G ? 0.1 * P_LS(phase_angle(G)) * 4 * cos(G.theta_i) : 0.0




# ---- Albedo computations ----

function spherealbedo(S::ScatteringLaw)
	s = 0.0
	for mu in linspace(0.0, 0.985, 10000)
		s += mu*Hemispheres.planaralbedo(S, acos(mu))
	end
	return 2s/10000
end

function geometricalbedo(S::ScatteringLaw)
	theta = linspace(0.0, 0.95*pi/2, 100)
	s = 0.0
	for th in theta
		s += cos(th) * value(S, Geometry(th,th,0.0))
	end
	return s / 100
end


# ---- Sphere integrated brightness, generic ----

function integrated(S::ScatteringLaw, alpha::Real)
	if alpha >= pi
		return 0.0
	end
	function integrand(x)
		G = from_latlon(x[1],x[2],alpha)
		value(S, G)*cos(G.theta_e)*cos(G.theta_i)*cos(x[1]) / 4
	end
	(val,err) = hcubature(integrand, (-pi/2, alpha-pi/2), (pi/2, pi/2),
						  reltol=1e-3)
	return val
end

	
end # module
