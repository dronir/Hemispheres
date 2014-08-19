module HemisphereComparator
	
using NetCDF

immutable ScatLaw
end

immutable Hemisphere
	nData::Int64
	NTheta::Int64
	dTheta::Float64
	NPhi::Array{Int64}
	dPhi::Array{Float64}
	cIdx::Array{Int64}
	dA::Array{Float64}
	data::Array{Float64}
end

Hemisphere(filename::String) = load_hemisphere(filename)

function generate_hemisphere(S::ScatLaw, Ntheta::Integer)
	nothing
end

function load_hemisphere(filename::String)
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

function plot_hemisphere(H::Hemisphere)
	nothing
end
	
end # module
