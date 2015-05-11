module Hemispheres

include("ScatteringGeometry.jl")
using .ScatteringGeometry
println("Loading NetCDF...")
using NetCDF
using Distributions

include("GaussQuadrature.jl/src/GaussQuadrature.jl")
using .GaussQuadrature

abstract ScatteringLaw
abstract AnalyticalScatteringLaw <: ScatteringLaw

export Hemisphere
export value, ratio, cell_index, point_in_cell
export generate_hemisphere, save_hemisphere, load_hemisphere, plot_hemisphere
export plot_primary_plane
export ScatteringLaw
export LommelSeeliger, ModifiedLommelSeeliger, BennuNominal, Akimov, Hapke
export Peltoniemi

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
end

# Construct Hemisphere directly from file.
Hemisphere(filename::String) = load_hemisphere(filename)


immutable Lambert <: AnalyticalScatteringLaw
end

immutable LommelSeeliger <: AnalyticalScatteringLaw
	omega::Float64
end
LommelSeeliger() = LommelSeeliger(1.0)

immutable Akimov <: AnalyticalScatteringLaw
end

immutable LSintegral <: AnalyticalScatteringLaw
end

immutable Lintegral <: AnalyticalScatteringLaw
end

immutable ModifiedLommelSeeliger <: AnalyticalScatteringLaw
	a::Float64
	b::Float64
	c::Float64
end

BennuNominal() = ModifiedLommelSeeliger(-4.36e-2, 2.69e-4, -9.90e-7)

immutable Hapke <: AnalyticalScatteringLaw
	w::Float64
end


immutable Minnaert <: AnalyticalScatteringLaw
	nu::Real
end


value(::Type{Lambert}, G::Geometry) = value(Lambert(), G)
value(S::Lambert, G::Geometry) = cos(G.theta_i)

value(::Type{LommelSeeliger}, G::Geometry) = value(LommelSeeliger(), G)
function value(S::LommelSeeliger, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
	return S.omega * mu0 / (mu+mu0) / (4*pi)
end

function value(S::ModifiedLommelSeeliger, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
	alpha = phase_angle(G)
	return mu0 / (mu+mu0) / (4*pi) * exp(S.a*alpha + S.b*alpha^2 + S.c*alpha^3)
end

function value(H::Hemisphere, G::Geometry)
	i,j = cell_index(H,G)
	return H.data[j,i]
end

value(::Type{Akimov}, G::Geometry) = value(Akimov(), G)
function value(S::Akimov, G::Geometry)
    alpha,beta,gamma = photometric_coordinates(G)
    return cos(alpha/2) * cos(pi/(pi-alpha)*(gamma - alpha/2)) / cos(gamma) * cos(beta)^(alpha/(pi-alpha))
end

function value(S::Hapke, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
    return 0.25 * mu0 / (mu + mu0) * HapkeH(S.w, mu) * HapkeH(S.w, mu0)
end

HapkeH(w,mu) = (1 + 2mu) / (1 + 2*mu*sqrt(1 - w))


function value(S::Minnaert, G::Geometry)
	mu0 = cos(G.theta_i)
	mu = cos(G.theta_e)
	return mu0^S.nu * mu^(S.nu-1)
end


value(::Type{LSintegral}, G::Geometry) = value(LSintegral(), G)
value(::Type{Lintegral}, G::Geometry) = value(Lintegral(), G)

function value(S::LSintegral, G::Geometry)
	alpha = phase_angle(G)
	return (1 - sin(alpha/2) * tan(alpha / 2) * log(cot(alpha / 4))) / 32
end

function value(S::Lintegral, G::Geometry)
	alpha = phase_angle(G)
	return (sin(alpha) + (pi - alpha) * cos(alpha)) / (6 * pi)
end


immutable AnalyticalLS <: AnalyticalScatteringLaw
	omega0::Float64
	P0::Float64
end
function value(S::AnalyticalLS, G::Geometry)
	alpha = phase_angle(G)
    omegaV = 0.5 * S.omega0 * S.P0
    PV = (1 - sin(alpha/2) * tan(alpha / 2) * log(cot(alpha / 4)))
	return 0.25 * omegaV * PV / (cos(G.theta_i) + cos(G.theta_e))
end

immutable AntiShadow <: AnalyticalScatteringLaw
    omega0::Float64
end
AntiShadow() = AntiShadow(0.1)
function value(S::AntiShadow, G::Geometry)
	alpha = phase_angle(G)
    mu0 = cos(G.theta_i)
    mu = cos(G.theta_e)
    omegaV = 2/3 * (1 - log(2)) * S.omega0
    PV = (1 - sin(alpha/2) * tan(alpha / 2) * log(cot(alpha / 4))) / (4/3 * (1 - log(2)))
	return 0.25/pi * mu0 / (mu + mu0) * omegaV * PV
end



# Compute the ratio between two hemispheres.
function ratio(A::Hemisphere, B::Hemisphere)
	if A.nTheta==B.nTheta
		return Hemisphere(A.nData, A.nTheta, A.dTheta, A.nPhi, A.dPhi, A.cIdx, A.dA, A.data./B.data)
	else
		error("Hemispheres don't match.")
	end
end



# Get cell index for given scattering direction.
function cell_index(H::Hemisphere, G::Geometry)
	t = int(fld(G.theta_e, H.dTheta)) + 1
	p = int(fld(G.phi, H.dPhi[t]))
	i = H.cIdx[t] + p
	j = int(fld(G.theta_i, H.dTheta)) + 1
	return (i,j)
end



# Get random point in given cell.
function point_in_cell(H::Hemisphere, idx::Integer)
	t = 0
	for t = 1:H.nTheta-1
		if H.cIdx[t+1] > idx
			break
		end
	end
	p = idx-H.cIdx[t]+1

    ca = cos((t-1) * H.dTheta)
    cb = cos(t * H.dTheta)

    theta = acos(ca - rand()*(ca - cb))
    phi   = (p - rand()) * H.dPhi[t]

	d = zeros(3)
    d[1] = sin(theta) * cos(phi)
    d[2] = sin(theta) * sin(phi)
    d[3] = cos(theta)
	return d
end



# This function generates a hemisphere with a given analytical scattering law.
function generate_hemisphere(S::AnalyticalScatteringLaw, nTheta::Integer, nSamples::Integer)
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
					#theta_i = min(theta_i, 1.5607966601082313)
					ca = cos((i-1)*dTheta)
					cb = cos(i*dTheta)
					theta_e = acos(ca - rand()*(ca-cb))
					#theta_e = min(theta_e, 1.5607966601082313)
					phi = (j-rand())*dPhi[i]
					G = Geometry(theta_i, theta_e, phi)
					data[k, cIdx[i]+j-1] += value(S, G)
				end
			end
		end
		data[1:end, cIdx[i]:cIdx[i]+nPhi[i]-1] /= N
	end
#	data /= nSamples
	Hemisphere(nBins, nTheta, dTheta, nPhi, dPhi, cIdx, dA, data)
end

generate_hemisphere(S::ScatteringLaw, nTheta::Integer) = generate_hemisphere(S, nTheta, 1000)



# This function loads a Hemisphere from a hemiScat NetCDF file.
function load_hemisphere(filename::String)
	foo = ncinfo(filename)
	nTheta = ncgetatt(filename, "Global", "nThetaI")
	dTheta = ncgetatt(filename, "Global", "dTheta")
	nPhi = ncgetatt(filename, "Global", "nPhi")
	dPhi = ncgetatt(filename, "Global", "dPhi")
	cIdx = ncgetatt(filename, "Global", "cIdx")
	dA = ncgetatt(filename, "Global", "dA")
	data = ncread(filename, "Hemisphere")
	data = squeeze(data,3)
	nData = sum(nPhi)
	ncclose()
	return Hemisphere(nData,nTheta,dTheta,nPhi,dPhi,cIdx,dA,data)
end



# This function saves a hemisphere to a (minimal) file that can be
# read with load_hemisphere().
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



# This function makes a plot of a given hemisphere.
plot_hemisphere(H::Hemisphere, thetaI::Real) = plot_hemisphere(H,thetaI,"",[])
plot_hemisphere(H::Hemisphere, thetaI::Real, saveplot::String) = plot_hemisphere(H,thetaI,saveplot,[])
function plot_hemisphere(H::Hemisphere, thetaI::Real, saveplot::String, vrange::Array)
    patches = Any[]
	const DEG = 57.295791433
    N = H.nTheta
    width = H.dTheta * DEG
    n = 0
    for i = 1:N
        r = i*H.dTheta * DEG
        for j = 1:H.nPhi[i]
            a = (j-1)*H.dPhi[i] * DEG - 0.5
            b = j*H.dPhi[i] * DEG + 0.5
            P1 = patch.Wedge((0,0), r, a-90, b-90, width=width*1.05)
            P2 = patch.Wedge((0,0), r, 360-b-90, 360-a-90, width=width*1.05)
			push!(patches, P1, P2)
            n += 2
		end
	end

	idx = int(fld(thetaI, H.dTheta))+1

    newdata = zeros(2*H.nData)
    for i = 1:H.nData
        newdata[2*i-1] = H.data[idx,i]
        newdata[2*i] = H.data[idx,i]
	end
	if vrange != []
		p = collections.PatchCollection(patches, cmap=colormap.gray, clim=vrange, linewidths=zeros(2*N))
	else
		p = collections.PatchCollection(patches, cmap=colormap.gray, linewidths=zeros(2*N))
	end
    p[:set_array](newdata[1:n])
    fig, ax = plot.subplots()
    ax[:add_collection](p)
    plot.xlim(-90, 90)
    plot.ylim(-90, 90)
    #plot.title("foo")
	#plot.colorbar(p)
	if saveplot != ""
		println("Saving to $saveplot...")
		plot.savefig(saveplot)
	else
		plot.show()
	end
end

# Make a plot of the primary plane
# (theta_e varies from -85 to +85 with theta_i constant)
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
function plot_primary_plane(H::Hemisphere, thetaI::Real, N::Integer)
	X,Y = get_primary_plane(H,thetaI,N)
	plot.plot(X,Y)
	plot.xlim(-90,90)
	plot.axvline(0,color="black")
	plot.show()
end

# Make a plot of the emergent plane
# (theta_i varies from -85 to +85 with theta_e constant)
function plot_emergent_plane(H::Hemisphere, theta_e::Real, N::Integer)
	X = linspace(-85.0, 85.0, N)
	Y = zeros(N)
	for i = 1:N
		x = X[i]
		theta_i = abs(x)*pi/180
		phi = x<=0 ? 0.0 : pi
		G = Geometry(theta_i, theta_e, phi)
		Y[i] = value(H,G)
	end
	plot.plot(X,Y)
	plot.xlim(-90,90)
	plot.axvline(0,color="black")
	plot.show()
	return X,Y
end

# Make a plot of the disk
# (theta_i and theta_e vary from -85 to +85 with alpha constant)
function plot_disk(H::Hemisphere, alpha::Real, N::Integer)
    alpha = alpha * 180/pi
	X = linspace(-85.0+alpha, 85.0, N)
	Y = zeros(N)
	for i = 1:N
		x = X[i]
		theta_i = abs(x)*pi/180
        theta_e = abs(x-alpha)*pi/180
		phi = theta_e<=0 && theta_i > 0 ? pi : 0.0
		G = Geometry(theta_i, theta_e, phi)
		Y[i] = value(H,G)
	end
	plot.plot(X,Y)
	plot.xlim(-90,90)
	plot.axvline(0,color="black")
	plot.show()
	return X,Y
end

plot_primary_plane(H::Hemisphere, thetaI::Real) = plot_primary_plane(H,thetaI,1000)


end # module
