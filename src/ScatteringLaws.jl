

module ScatteringLaws

using Abstracts
using ScatteringGeometry
using PhaseFunctions
using Hemispheres

using Cubature

import PhaseFunctions.value

export ScatteringLaw, PhaseFunctions
export value, Lambert, LommelSeeliger, ParticulateMedium, AntiShadow, AntiR
export sphere_albedo, geometric_albedo, integrated
export Hemisphere
export generate_hemisphere, save_hemisphere, load_hemisphere, plot_hemisphere
export Geometry


# ---- Lambert ----

immutable Lambert <: AnalyticalScatteringLaw
end

value(::Type{Lambert}, G::Geometry) = value(Lambert(), G)
value(S::Lambert, G::Geometry) = bool(G) ? 1.0 : 0.0


planaralbedo(S::Lambert, theta::Real) = 1.0
geometric_albedo(S::Lambert) = 2.0/3.0
sphere_albedo(S::Lambert) = 1.0
integrated(S::Lambert, alpha::Real) = (sin(alpha) + (pi-alpha)*cos(alpha)) / (6pi)


# ---- Lommel-Seeliger ----

immutable LommelSeeliger <: AnalyticalScatteringLaw
	P::PhaseFunction
	omega::Float64
end
LommelSeeliger() = LommelSeeliger(Isotropic(), 1.0)

value(::Type{LommelSeeliger}, G::Geometry) = value(LommelSeeliger(), G)
function value(S::LommelSeeliger, G::Geometry)
	return bool(G) ? S.omega * value(S.P, G) / (cos(G.theta_e)+cos(G.theta_i)) / 4 : 0.0
end

function planaralbedo(S::LommelSeeliger, theta::Real) 
	mu = cos(theta)
	mu<eps(theta) ? 0.5*S.omega : S.omega*(mu*log(mu) - mu*log(mu+1) + 1) / 2
end

geometric_albedo(S::LommelSeeliger) = S.omega / 8 * value(S.P, 0.0)
sphere_albedo(S::LommelSeeliger) = 2/3 * (1 - log(2))
integrated(S::LommelSeeliger, alpha::Real) = value(S.P, alpha)/32 * (alpha < eps() ? 1.0 : 1 - sin(alpha/2) * tan(alpha/2) * log(cot(alpha/4)))



# ---- Particulate Medium ----

immutable ParticulateMedium <: ScatteringLaw
	hemi::Hemisphere
	P::PhaseFunction
end
value(S::ParticulateMedium, G::Geometry) = bool(G) ? value(S.hemi,G) * value(S.P, G) : 0.0


# ---- Anti-Shadow ----

# Dividing this out of a hemiScatter output hemisphere
# gives the shadowing correction S.

immutable AntiShadow <: AnalyticalScatteringLaw
    omega0::Float64
end
AntiShadow() = AntiShadow(0.1)
function value(S::AntiShadow, G::Geometry)
	if !bool(G)
		return 0.0 
	end
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

immutable AntiR <: AnalyticalScatteringLaw 
	C::Float64
end
AntiR() = AntiR(1/(2pi))
function value(S::AntiR, G::Geometry) 
	if !bool(G) 
		return 0.0 
	end
	alpha = phase_angle(G)
	phiLS = (1 - sin(alpha/2) .* tan(alpha / 2) .* log(cot(alpha / 4)))
	return phiLS * S.C * cos(G.theta_i)
end




# ---- Albedo computations ----

function sphere_albedo(S::ScatteringLaw)
	s = 0.0
	for mu in linspace(0.0, 0.985, 10000)
		s += mu*Hemispheres.planaralbedo(S, acos(mu))
	end
	return 2s/10000
end

geometric_albedo(S::ScatteringLaw) = 4 * integrated(S, 0.0)


# ---- Sphere integrated brightness, generic ----

function integrated(S::ScatteringLaw, alpha::Real)
	if alpha >= pi
		return 0.0
	end
	function integrand(x)
		G = from_latlon(x[1],x[2],alpha)
		value(S, G)*cos(G.theta_e)*cos(G.theta_i)*cos(x[1])
	end
	(val,err) = hcubature(integrand, (-pi/2, alpha-pi/2), (pi/2, pi/2),
						  reltol=1e-3)
	return val / (4pi)
end

	
end # module
