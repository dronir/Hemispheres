module ScatteringGeometry

export Geometry, phase_angle, beta_angle, gamma_angle, photometric_coordinates

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
	phi = acos((cosa - mu0*mu) / (sin(theta_i) * sin(theta_e)))
	return Geometry(theta_i, theta_e, phi)
end


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


end # module
