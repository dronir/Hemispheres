

module ScatteringLaws

using Abstracts
using ScatteringGeometry
using PhaseFunctions
using Hemispheres

using Distributions
using Cubature

import PhaseFunctions.value

export ScatteringLaw, PhaseFunctions, AnalyticalScatteringLaw
export value, Lambert, LommelSeeliger, ParticulateMedium, Shadow, AntiR
export sphere_albedo, geometric_albedo, integrated
export Hemisphere
export Peltoniemi
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



# ---- Peltoniemi surface ----

const N1 = Normal(0,1)

immutable Peltoniemi <: AnalyticalScatteringLaw
    R::ScatteringLaw
    rho::Float64
end

Peltoniemi() = Peltoniemi(LommelSeeliger(), 0.5)

const PX =
[-5.38748,-4.60368,-3.94476,-3.34785,-2.78881,-2.25497,-1.73854,-1.23408,
-0.737474,-0.245341,0.245341,0.737474,1.23408,1.73854,2.25497,2.78881,3.34785,
3.94476,4.60368,5.38748]

const Pw =
[2.22939e-13,4.39934e-10,1.08607e-7,7.80256e-6,0.000228339,0.00324377,0.0248105,
0.109017,0.286676,0.462244,0.462244,0.286676,0.109017,0.0248105,0.00324377,
0.000228339,7.80256e-6,1.08607e-7,4.39934e-10,2.22939e-13]




function value(S::Peltoniemi, G::Geometry)
	W = PeltoniemiW(S,G)
	integral = 0.0
	rndphi = 2*pi*rand()
	vecI = [sin(G.theta_i)*cos(rndphi), sin(G.theta_i)*sin(rndphi), cos(G.theta_i)]
	vecE = [sin(G.theta_e)*cos(G.phi+rndphi), sin(G.theta_e)*sin(G.phi+rndphi), cos(G.theta_e)]
	f = 0.0
	T = [0.0, 0.0, 1.0]
	for i = 1:20
		T[1] = PX[i] * sqrt(2)*S.rho
		w1 = Pw[i]
		for j = 1:20
			T[2] = PX[j] * sqrt(2)*S.rho
			f = PeltoniemiIntegrand(S,vecI,vecE,T)
			integral += f * w1 * Pw[j]
		end
	end
	return W * integral / sqrt(pi)
end

function PeltoniemiIntegrand(S::Peltoniemi, vecI::Vector, vecE::Vector, T::Vector)
	t = norm(T)
	mux0 = dot(vecI, T) / t
	mux = dot(vecE, T) / t
	if mux0 > 0 && mux > 0
		p_i = vecI - T*mux0
		p_e = vecE - T*mux
		if norm(p_i) > 0 && norm(p_e) > 0
			phi = acos(dot(p_i, p_e) / norm(p_i) / norm(p_e))
		else
			phi = 0.0
		end
		Gx = Geometry(acos(mux0), acos(mux), phi)
		return value(S.R, Gx) * 1 * mux * t / (vecI[3]*vecE[3])
	else
		return 0.0
	end
end

function PeltoniemiG(invxi::Real)
	if invxi > 0
		xi = 1/invxi
		n = pdf(N1, xi)
		N = cdf(N1, -xi)
		return n/xi - N
	else
		return 1.0
	end
end

function PeltoniemiW(S::Peltoniemi, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
	invxi = (S.rho * sqrt(1-mu^2)) / mu
	invxi0 = (S.rho * sqrt(1-mu0^2)) / mu0
	V = 1 / (1 + PeltoniemiG(invxi))
	V0 = 1 / (1 + PeltoniemiG(invxi0))
	s = sqrt(sin(G.phi/2))
	return min(V,V0) * (1-s) + V*V0*s
end




# ---- Particulate Medium ----

immutable ParticulateMedium <: ScatteringLaw
    hemi::Hemisphere
    P::PhaseFunction
end
value(S::ParticulateMedium, G::Geometry) = bool(G) ? value(S.hemi,G) * value(S.P, G) : 0.0
integrated(S::ParticulateMedium, a::Real) = integrated(S.hemi, a) * value(S.P, a)


# ---- Shadow ----

# Gives the value of the shadowing correction
# when H is a raytracing output hemisphere.

immutable Shadow <: ScatteringLaw
    omega0::Float64
    H::Hemisphere
end

function value(S::Shadow, G::Geometry)
    if !bool(G)
        return 0.0 
    end
    alpha = phase_angle(G)
    mu0 = cos(G.theta_i)
    mu = cos(G.theta_e)
    return 8 * value(S.H, G) * (mu+mu0) / (S.omega0 * mu0 * phi_LS(alpha))
end

phi_LS(alpha) = (1 - sin(alpha/2) .* tan(alpha / 2) .* log(cot(alpha / 4)))



# ---- Anti-R ----

# Dividing this out of a hemiScatter output hemisphere
# gives a reflection coefficient, sans phase function.

immutable AntiR <: AnalyticalScatteringLaw 
    C::Float64
end
AntiR() = AntiR(1.0)
function value(S::AntiR, G::Geometry) 
    if !bool(G) 
        return 0.0 
    end
    alpha = phase_angle(G)
    return S.C * cos(G.theta_i) * phi_LS(alpha) / 0.5
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
    function integrand{T<:Real}(x::Array{T})
        G = from_latlon(x[1],x[2],alpha)
        value(S, G)*cos(G.theta_e)*cos(G.theta_i)*cos(x[1])
    end
    (val,err) = hcubature(integrand, (-pi/2, alpha-pi/2), (pi/2, pi/2),
                          reltol=1e-3)
    return val / (4pi)
end

    
end # module
