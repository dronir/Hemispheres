module HemispherePlots

using ScatteringLaws
using PyCall

@pyimport matplotlib.pyplot as plot


# This function makes a plot of a given hemisphere.
plot_hemisphere(H::ScatteringLaw, thetaI::Real) = plot_hemisphere(H,thetaI,"",[])
plot_hemisphere(H::ScatteringLaw, thetaI::Real, saveplot::String) = plot_hemisphere(H,thetaI,saveplot,[])
function plot_hemisphere(H::ScatteringLaw, thetaI::Real, saveplot::String, vrange::Array)
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
function get_primary_plane(H::ScatteringLaw, thetaI::Real, N::Integer)
	X = linspace(-85.0, 85.0, N)
	Y = zeros(N)
	for i = 1:N
		x = X[i]
		theta_e = abs(x)*pi/180
		if typeof(H) == Hemisphere
			t = int(fld(theta_e, H.dTheta))+1
			phi = x<=0 ? 10*eps() : pi-0.5*H.dPhi[t]
		else
			phi = x<=0 ? 10*eps() : pi-0.001
		end
		G = Geometry(thetaI, theta_e, phi)
		Y[i] = value(H,G)
	end
	return X,Y
end

function plot_primary_plane(H::ScatteringLaw, thetaI::Real, N::Integer)
	X,Y = get_primary_plane(H,thetaI,N)
	plot.plot(X,Y)
	plot.xlim(-90,90)
	plot.axvline(0,color="black")
	plot.show()
end

function plot_primary_plane{T<:ScatteringLaw}(A::Array{T,1}, thetaI::Real, N::Integer)
	for H in A
		X,Y = get_primary_plane(H,thetaI,N)
		plot.plot(X,Y)
	end
	plot.xlim(-90,90)
	plot.axvline(0,color="black")
	plot.show()
end

# Make a plot of the emergent plane
# (theta_i varies from -85 to +85 with theta_e constant)
function plot_emergent_plane(H::ScatteringLaw, theta_e::Real, N::Integer)
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
function plot_disk(H::ScatteringLaw, alpha::Real, N::Integer)
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


function plot_phase(H::ScatteringLaw, N::Integer)
    Geometries = Geometry[]
    Thetas = linspace(0.0, 80.0, N) * pi/180
    
    for i = 1:N
        theta_i = Thetas[i]
        for j = 1:N
            theta_e = Thetas[j]
            kmax = fld(N^2-N,2)
            Phis = linspace(0.0, 180.0, kmax) * pi/180
            for k = 1:kmax
                phi = Phis[k]
                push!(Geometries, Geometry(theta_i, theta_e, phi))
            end
        end
    end
    Alphas = [ScatteringLaws.phase_angle(g) for g in Geometries]
    Values = [value(H,g) for g in Geometries]
    plot.scatter(Alphas, Values)
    plot.show()
    
end


end # module
