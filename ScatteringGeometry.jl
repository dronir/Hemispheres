module ScatteringGeometry

export Geometry

# Illumination geometry type
immutable Geometry
	theta_i::Float64
	theta_e::Float64
	phi::Float64
end

end # module
