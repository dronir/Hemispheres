module ScatteringGeometry

export Geometry, phase_angle

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

end # module
