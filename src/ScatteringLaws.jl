

module ScatteringLaws

#include("Abstracts.jl")
#include("Abstracts.jl")
#include("Abstracts.jl")
using Abstracts
using ScatteringGeometry
using PhaseFunctions
using Hemispheres

import PhaseFunctions.value

export value, Lambert, LommelSeeliger, ParticulateMedium, AntiShadow, AntiR
export spherealbedo, geometricalbedo


# ---- Lambert ----

immutable Lambert <: AnalyticalScatteringLaw
end

value(::Type{Lambert}, G::Geometry) = value(Lambert(), G)
value(S::Lambert, G::Geometry) = 1.0

planaralbedo(S::Lambert, theta::Real) = 1.0
geometricalbedo(S::Lambert) = 2/3
spherealbedo(S::Lambert) = 2/3


# ---- Lommel-Seeliger ----

immutable LommelSeeliger <: AnalyticalScatteringLaw
	P::PhaseFunction
	omega::Float64
end
LommelSeeliger() = LommelSeeliger(Isotropic(), 1.0)

value(::Type{LommelSeeliger}, G::Geometry) = value(LommelSeeliger(), G)
function value(S::LommelSeeliger, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
	return S.omega * value(S.P, G) / (mu+mu0) / (4pi)
end

function planaralbedo(S::LommelSeeliger, theta::Real) 
	mu = cos(theta)
	mu<eps(theta) ? 0.5*S.omega : S.omega*(mu*log(mu) - mu*log(mu+1) + 1) / 2
end

geometricalbedo(S::LommelSeeliger) = S.omega / 8 * value(S.P, 0.0)
spherealbedo(S::LommelSeeliger) = 2/3 * (1 - log(2))


# ---- Particulate Medium ----

immutable ParticulateMedium <: ScatteringLaw
	hemi::Hemisphere
	P::PhaseFunction
end
value(S::ParticulateMedium, G::Geometry) = value(S.hemi,G) * value(S.P, G)


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
function value(S::AntiR, G::Geometry)
	mu0 = cos(G.theta_i)
	alpha = phase_angle(G)
	return 0.1 * P_LS(alpha) * 4 * mu0
end



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


# ---- Sphere integrated brightness ----

const Niter = 100000

function integrated(S::ScatteringLaw, alpha::Real)
	s = 0.0
	for i = 1:Niter
		G = random_geometry(alpha)
		s += value(S,G) * cos(G.theta_i) * cos(G.theta_e)
	end
	pi*s/Niter
end

	
end # module
