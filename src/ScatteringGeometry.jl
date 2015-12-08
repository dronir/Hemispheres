module ScatteringGeometry

export Geometry, phase_angle, beta_angle, gamma_angle, photometric_coordinates
export random_geometry, from_latlon

import Base.bool

# Illumination geometry type
immutable Geometry
	theta_i::Float64
	theta_e::Float64
	phi::Float64
end

function Geometry(i::Vector, e::Vector, n::Vector)
	mu0 = dot(i, n)
	mu = dot(e, n)
	cosa = dot(i, e)
	theta_i = acos(mu0)
	theta_e = acos(mu)
	if theta_e==0 && theta_i==0
		return Geometry(0.0, 0.0, 0.0)
	end
	cosphi = (cosa - mu0*mu) / (sin(theta_i) * sin(theta_e))
	cosphi = clamp(cosphi, -1.0, 1.0)
	phi = acos(cosphi)
	return Geometry(theta_i, theta_e, phi)
end

bool(G::Geometry) = cos(G.theta_i) > 0.0 && cos(G.theta_e) > 0.0


function phase_angle(G::Geometry)
	i = G.theta_i
	e = G.theta_e
	p = G.phi
	cosa = cos(i)*cos(e) + sin(i)*sin(e)*cos(p)
	cosa = clamp(cosa, -1+eps(), 1-eps())
	return acos(cosa)
end

function beta_angle(G::Geometry)
    i = G.theta_i
    e = G.theta_e
    p = G.phi
    top = sin(i+e)^2 - cos(p/2)^2 * sin(2*e) * sin(2*i)
    bottom = top + sin(e)^2 * sin(i)^2 * sin(p)^2
    return acos(sqrt(top/bottom))
end

gamma_angle(G::Geometry) = acos(cos(G.theta_e) / cos(beta_angle(G)))

function photometric_coordinates(G::Geometry)
    ix = G.theta_i
    ex = G.theta_e
    px = G.phi
    cosalpha = cos(ix)*cos(ex) + sin(ix)*sin(ex)*cos(px)
	alpha = acos(cosalpha)
    tanbeta = (cos(ix)/cos(ex) - cosalpha) / sin(alpha)
	beta = atan(tanbeta)
    cosgamma = cos(ex) / cos(beta)
    return (alpha, beta, acos(cosgamma))
end

random_geometry() = random_geometry(0.0)
function random_geometry(alpha::Real)
	lat = asin(2*(rand()-0.5))
	lon = (alpha - pi/2) + rand() * (pi - alpha)
	mu = cos(lon)*cos(lat)
	mu0 = cos(lon-alpha)*cos(lat)
	theta_e = acos(mu)
	theta_i = acos(mu0)
	cosphi = abs(alpha)<eps() ? 1.0 : (cos(alpha) - mu0*mu) / (sin(theta_i) * sin(theta_e))
	phi = acos(cosphi)
	Geometry(theta_i, theta_e, phi)
end

function from_latlon(lat::Real, lon::Real, alpha::Real)
	mu0 = cos(lon)*cos(lat)
	mu = cos(lon-alpha)*cos(lat)
	theta_e = acos(mu)
	theta_i = acos(mu0)
	if theta_e==0 && theta_i==0
		return Geometry(0.0, 0.0, 0.0)
	end
	cosphi = (cos(alpha) - mu0*mu) / (sin(theta_i) * sin(theta_e))
	cosphi = clamp(cosphi, -1.0, 1.0)
	phi = acos(cosphi)
	Geometry(theta_i, theta_e, phi)
end


end # module
