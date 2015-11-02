module Hemispheres

using Abstracts
using ScatteringGeometry

using NetCDF
using Distributions

export Hemisphere
export value, cell_index, point_in_cell
export generate_hemisphere, save_hemisphere, load_hemisphere, plot_hemisphere
export plot_primary_plane

# Hemisphere type
type Hemisphere <: ScatteringLaw
	nData::Int64
	nTheta::Int64
	dTheta::Float64
	nPhi::Array{Int64}
	dPhi::Array{Float64}
	cIdx::Array{Int64}
	dA::Array{Float64}
	data::Array{Float64}
	planar::Array{Float64}
end

# Construct Hemisphere directly from file.
Hemisphere(filename::String) = load_hemisphere(filename)


import PhaseFunctions.value


# ---- Product and division of Hemispheres ----

import Base.*, Base./

function *(A::Hemisphere, B::Hemisphere)
	if A.nTheta==B.nTheta
		H = Hemisphere(A.nData, A.nTheta, A.dTheta, A.nPhi, A.dPhi, A.cIdx, A.dA, A.data.*B.data, A.planar)
		compute_planar!(H)
		return H
	else
		error("Hemispheres don't match.")
	end
end

function /(A::Hemisphere, B::Hemisphere)
	if A.nTheta==B.nTheta
		H = Hemisphere(A.nData, A.nTheta, A.dTheta, A.nPhi, A.dPhi, A.cIdx, A.dA, A.data./B.data, A.planar)
		compute_planar!(H)
		return H
	else
		error("Hemispheres don't match.")
	end
end



# ---- Retrieving the value of a Hemisphere in a given geometry ----

function value(H::Hemisphere, G::Geometry)
	if !bool(G)
		return 0.0
	end
	idx_in, idx_out = cell_index(H,G)
	return H.data[idx_in, idx_out]
end

function cell_index(H::Hemisphere, G::Geometry)
	idx_in = int(fld(G.theta_i, H.dTheta)) + 1
	t = int(fld(G.theta_e, H.dTheta)) + 1
	p = int(fld(G.phi, H.dPhi[t]))
	idx_out = H.cIdx[t] + p
	return idx_in, idx_out
end



# ---- Generating a Hemisphere from a ScatteringLaw ----

function generate_hemisphere(S::ScatteringLaw, nTheta::Integer, nSamples::Integer)
	cIdx = ones(Int64, nTheta)
	nPhi = ones(Int64, nTheta)
	dPhi = zeros(Float64, nTheta)
	dA = zeros(Float64, nTheta)

	dTheta = pi/(2*nTheta)
	dA0 = pi*(1 - cos(dTheta))
	dA[1] = dA0
	dPhi[1] = pi
	for i = 2:nTheta
		dp = dA0 / (cos((i-1)*dTheta) - cos(i*dTheta))
		nPhi[i] = int(pi / dp)
		dPhi[i] = pi / nPhi[i]
		dA[i] = dPhi[i] * (cos((i-1)*dTheta) - cos(i*dTheta))
		cIdx[i] = cIdx[i-1] + nPhi[i-1]
	end
	nBins = sum(nPhi)
	data = zeros(nTheta, nBins)

	# go through each bin and compute scattering law values
	for i = 1:nTheta
		N = int((dA[i] / dA[1])*nSamples)
		for j = 1:nPhi[i]
			for k = 1:nTheta
				for n = 1:N
					theta_i = (k-rand())*dTheta
					ca = cos((i-1)*dTheta)
					cb = cos(i*dTheta)
					theta_e = acos(ca - rand()*(ca-cb))
					phi = (j-rand())*dPhi[i]
					G = Geometry(theta_i, theta_e, phi)
					data[k, cIdx[i]+j-1] += value(S, G)
				end
			end
		end
		data[1:end, cIdx[i]:cIdx[i]+nPhi[i]-1] /= N
	end
	H = Hemisphere(nBins, nTheta, dTheta, nPhi, dPhi, cIdx, dA, data, zeros(nTheta))
	compute_planar!(H)
	return H
end

generate_hemisphere(S::ScatteringLaw, nTheta::Integer) = generate_hemisphere(S, nTheta, 1000)




# ---- Planar albedos of Hemispheres ----

function compute_planar!(H::Hemisphere)
	for i = 1:H.nTheta
		theta = (i-0.5) * H.dTheta
		H.planar[i] = raw_planaralbedo(H, theta)
	end
end

function planaralbedo(H::Hemisphere, theta_i::Real)
	theta_i = abs(theta_i)
	if theta_i > pi/2
		return 0.0
	end
	idx = int(fld(theta_i, H.dTheta)) + 1 
	return H.planar[idx]
end


function raw_planaralbedo(H::Hemisphere, theta_i::Real)
	Ap = 0.0
	theta_i = abs(theta_i)
	if theta_i > pi/2
		return 0.0
	end
	idx = int(fld(theta_i, H.dTheta)) + 1 
	for i = 1:H.nTheta
		theta_e = (i-1)*H.dTheta
		# dA = sin(theta) dtheta dphi, basically
		t = cos(theta_e) * H.dA[i]
		for j = 0:H.nPhi[i]-1
			Ap += t * H.data[idx, H.cIdx[i]+j]
		end
	end
	return 2 * Ap / pi
end




# ---- I/O of Hemispheres ----

function load_hemisphere(filename::String)
	foo = ncinfo(filename)
	dTheta = ncgetatt(filename, "Global", "dTheta")
	nPhi = ncgetatt(filename, "Global", "nPhi")
	dPhi = ncgetatt(filename, "Global", "dPhi")
	cIdx = ncgetatt(filename, "Global", "cIdx")
	dA = ncgetatt(filename, "Global", "dA")
	data = ncread(filename, "Hemisphere")
	if ndims(data) == 3
		data = squeeze(data,3)
	end
	nTheta = size(data)[1]
	nData = sum(nPhi)
	ncclose()
	planar = zeros(nTheta)
	H = Hemisphere(nData,nTheta,dTheta,nPhi,dPhi,cIdx,dA,data,planar)
	compute_planar!(H)
	return H
end


function save_hemisphere(H::Hemisphere, filename::String)
	level = NcDim("level", 1)
	nData = NcDim("data", H.nData)
	thetaI = NcDim("thetaI", H.nTheta)
	data = NcVar("Hemisphere", [thetaI, nData, level], t=Float64)
	globals = {
		"nThetaI"=>H.nTheta,
		"dTheta"=>H.dTheta,
		"nPhi"=>H.nPhi,
		"dPhi"=>H.dPhi,
		"cIdx"=>H.cIdx,
		"dA"=>H.dA
	}
	D = ones(Float64, H.nTheta, sum(H.nPhi), 1)
	D[:,:,1] = H.data
	nc = NetCDF.create(filename, [data], mode=NC_CLASSIC_MODEL)
	NetCDF.putatt(nc, "", globals)
	NetCDF.putvar(nc, "Hemisphere", D)
	NetCDF.close(nc)
end

function get_primary_plane(H::Hemisphere, thetaI::Real, N::Integer)
	X = linspace(-85.0, 85.0, N)
	Y = zeros(N)
	for i = 1:N
		x = X[i]
		theta_e = abs(x)*pi/180
		t = int(fld(theta_e, H.dTheta))+1
		phi = x<=0 ? 10*eps() : pi-0.5*H.dPhi[t]
		G = Geometry(thetaI, theta_e, phi)
		Y[i] = value(H,G)
	end
	return X,Y
end



end # module
