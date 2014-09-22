module ScatteringGeometry

export Geometry, phase_angle, beta_angle, gamma_angle, photometric_coordinates

# Illumination geometry type
immutable Geometry
	theta_i::Float64
	theta_e::Float64
	phi::Float64
end

function phase_angle(G::Geometry)
	i = G.theta_i
	e = G.theta_e
	p = G.phi
	return acos(cos(i)*cos(e) + sin(i)*sin(e)*cos(p))
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
    i = G.theta_i
    e = G.theta_e
    p = G.phi
    top = sin(i+e)^2 - cos(p/2)^2 * sin(2*e) * sin(2*i)
    bottom = top + sin(e)^2 * sin(i)^2 * sin(p)^2
    cosalpha = cos(i)*cos(e) + sin(i)*sin(e)*cos(p)
    cosbeta = sqrt(top/bottom)
    cosgamma = cos(e) / cosbeta
    return (acos(cosalpha), acos(cosbeta), acos(cosgamma))
end


end # module
