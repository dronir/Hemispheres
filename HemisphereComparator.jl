module HemisphereComparator
	
using NetCDF

abstract ScatteringLaw
abstract AnalyticalScatteringLaw

immutable Hemisphere <: ScatteringLaw
	nData::Int64
	NTheta::Int64
	dTheta::Float64
	NPhi::Array{Int64}
	dPhi::Array{Float64}
	cIdx::Array{Int64}
	dA::Array{Float64}
	data::Array{Float64}
end

# Construct directly from file.
Hemisphere(filename::String) = load_hemisphere(filename)

function ratio(A::Hemisphere, B::Hemisphere) 
	if A.NTheta==B.NTheta
		return Hemisphere(A.nData, A.NTheta, A.dTheta, A.NPhi, A.dPhi, A.cIdx, A.dA, A.data./B.data)
	else
		error("Hemispheres don't match.")
	end
end

# This function generates a hemisphere with a given analytical scattering law.
function generate_hemisphere(S::AnalyticalScatteringLaw, Ntheta::Integer)
	nothing
end

# This function loads a Hemisphere from a hemiScat NetCDF file.
function load_hemisphere(filename::String)
	foo = ncinfo(filename)
	nTheta = ncgetatt(filename, "Global", "nThetaI")
	dTheta = ncgetatt(filename, "Global", "dTheta")
	NPhi = ncgetatt(filename, "Global", "nPhi")
	dPhi = ncgetatt(filename, "Global", "dPhi")
	cIdx = ncgetatt(filename, "Global", "cIdx")
	dA = ncgetatt(filename, "Global", "dA")
	data = ncread(filename, "Hemisphere")
	data = squeeze(data,3)
	nData = sum(NPhi)
	return Hemisphere(nData,nTheta,dTheta,NPhi,dPhi,cIdx,dA,data)
end

# This function saves a hemisphere to a (minimal) file that can be
# read with load_hemisphere().
function save_hemisphere(H::Hemisphere, filename::String)
	nothing
end

# This function makes a plot of a given hemisphere.
function plot_hemisphere(H::Hemisphere, thetaI::Real)
	nothing
end
	
end # module
