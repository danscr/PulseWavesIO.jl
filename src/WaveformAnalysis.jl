function GaussianDecomposition(wave::Vector{UInt8}; smooth::Bool=true, width::Number=1, rel_threshold::Number=0.22, min_peakprominence::Number=1, abs_threshold::Number=4, sigma_min::Float64=2.) # the absolute threshold is based on Riegl Q780 LUT data, where the first 4 entries (= digital numbers of the signal) are -infinity dB, i.e. below the detection threshold
  wave = Float64.(wave .- minimum(wave)) # If min is not set to 0 the optimization will find one infinite width component that we don't want
  wave[ wave .<= abs_threshold] .= 0.
  ### Smoothing
  if smooth
    wave_smooth = ImageFiltering.imfilter(wave, ImageFiltering.Kernel.gaussian((width,)))
  else
    wave_smooth = wave
  end
  ### Peak detection
  peaks, amplitudes = Peaks.findmaxima(wave_smooth)
  _, peakproms = Peaks.peakproms!(peaks, wave_smooth; minprom=min_peakprominence) #see Peaks.jl documentation https://github.com/halleysfifthinc/Peaks.jl
  # peak prominence > 1 eliminates most noise peaks, but 1D Gaussian smoothing is the safer bet
  maxIntensity = maximum(wave)
  peaks = filter(peak -> (wave[peak] > rel_threshold * maxIntensity) & (wave[peak] > abs_threshold), peaks)
  n_peaks = length(peaks)
  @assert(n_peaks > 0, "No peaks found")
  ### Gaussian parameter fitting
  p0 = vcat([[wave[peak], peak, 8.] for peak in peaks]...) # initial parameter guess needs to be 1D for some reason... The third element is the standard deviation
  p_lower = vcat([[0., 0., sigma_min] for peak in peaks]...) # set lower boundary for standard deviation parameter
  LMfit = LsqFit.curve_fit(GaussianSum, 1:length(wave), wave, p0; lower=p_lower, autodiff=:finite) #fit Gaussian components with Levenberg-Marquart method, and calculate the Jacobian with finite difference method
  #@assert(LMfit.converged, "Gaussian decomposition failed. NLS did not converge.")
  params = reshape(LMfit.param, (3,n_peaks)) # reshape to 3xN (amplitude, mu, sigma)xN
  # drop peaks that were fitted with zero amplitude
  params = params[:,findall((params[1,:] .!= 0.) .& (params[1,:] .<= 255))]
  @assert(all(params .> 0.), "Non-positive parameters detected.")
  bias = mean(LMfit.resid)
  rmse = sqrt(mean(LMfit.resid .^ 2))
  params, bias, rmse, LMfit.converged
end

function IterativelySmoothedGaussianDecomposition(wave::Vector{UInt8}, width::Number=1, increment::Float64=1.1, width_threshold::Number=15; kwargs...)
  while true
    p, bias, rmse, converged = GaussianDecomposition(wave; width=width, kwargs...)
    if converged
      return p, bias, rmse, width
    else
      width *= increment
      if width > width_threshold
	throw("Iterative smoothing exceeded blurring width threshold")
      end
    end
  end
end

function CheckDecomposition(waves, nmin, nmax)
  for i in nmin:nmax
    wave = waves[i]
    for j in 1:length(wave)
      try
	p, bias, rmse, w = PulseWavesIO.IterativelySmoothedGaussianDecomposition(wave[j].Samples; smooth=true, width=1, rel_threshold=.1, abs_threshold=4)
        f= Figure(resolution=(1600,600))
	ax=Axis(f[1,1:5], limits=(nothing, (0,150)))
	lines!(1:.01:length(wave[j].Samples), GaussianSum(1:.01:length(wave[j].Samples), p[:]), color=:firebrick, label="Decomposition")
        lines!(wave[j].Samples, color=:royalblue, linestyle="--", label="Original waveform (with noise)"); hlines!([4], color=:black, linestyle="..", label="Absolute noise threshold")
	#pars = reshape(p, (3,div(length(p),3)+1))
	for k in 1:size(p)[2]
	  lines!(repeat([p[2,k]],2), [0,p[1,k]], linewidth=1, color=:black)
	  fwhm = p[3,k] .* 2.355
	  lines!([p[2,k] .- fwhm ./ 2, p[2,k] .+ fwhm ./ 2],repeat([p[1,k]],2) ./ 2, linewidth=1, color=:black)
	end
	axislegend()
	text!(0,120,text="Bias = $(bias)\nRMSE = $(rmse)\nGaussian blur width = $w\nParameters : $(p)")
        return f#GLMakie.save("dev/system_pulse_decomposition/wavedecomposition_$(i)_$(j).png", f)
      catch e
	println("wave $(i) segment $(j) failed: $(e)")
	f= Figure(resolution=(1600,600))
	ax=Axis(f[1,1:5], limits=(nothing, (0,255)))
	#lines!(1:.01:length(wave[j].Samples), PulseWavesIO.GaussianSum(1:.01:length(wave[j].Samples), p[:]), color=:firebrick, label="Decomposition")
        lines!(wave[j].Samples, color=:royalblue, linewidth=5, linestyle="--", label="Original waveform (with noise)"); hlines!([4], color=:black, linewidth=5, linestyle="..", label="Absolute noise threshold")
	axislegend()
	text!(30,120,text="Decomposition failed!", fontsize=30)
        return f#GLMakie.save("dev/system_pulse_decomposition/wavedecomposition_$(i)_$(j).png", f)
      end
    end
  end
end


function GaussianPDFWithAmplitude(x::Int64, A::Float64, mu::Float64, sigma::Float64)
  A * exp(-(x - mu)^2 / (2 * sigma^2))
end

function GaussianSum(x::UnitRange{Int64}, p::Vector{Float64})
  n = Int(length(p) / 3)
  #sum([GaussianPDFWithAmplitude.(x, p[i,1], p[i,2], p[i,3]) for i in 1:size(p)[1]])
  sum([GaussianPDFWithAmplitude.(x, p[(i-1)*3 + 1], p[(i-1)*3 + 2], p[(i-1)*3 + 3]) for i in 1:n])
  #sum([GaussianPDFWithAmplitude.(x, p[1], p[2], p[3]), GaussianPDFWithAmplitude.(x,p[4],p[5],p[6])])
end

function AdditiveGaussianDistributionPDF(x::Number, distributions::Vector{Distributions.Normal{Float64}}, scale_factors::Union{Float64,Vector{Float64}}; upper::Number)
  if x > upper
    return NaN
  end
  pdfs = pdf.(distributions, x)
  sum(pdfs .* scale_factors)
end


"""
RayGroundIntersection( ray, dem, dem_header, threshold=0.5 )
ray is a 2x3 matrix (1x3 origin, 1x3 direction vectors)
dem is a 2D array
dem_header is a Dict object with entries "YLLCORNER", "XLLCORNER", "NODATA_VALUE" and "CELLSIZE"
threshold is the height above ground that is to be considered the classification threshold (in meters). Defaults to 0.5.
resolution is the step size that is used to find the intersection between DEM and ray (in sampling units). Defaults to 1.

This function returns the number of units along the direction vector of the ray are required to travel until the ground (within the threshold) is reached.
"""
function RayGroundIntersection(ray::Matrix{Float64}, dem::Dict, threshold::Float64=1.5, resolution::Float64=1.) # If this method is too slow, try Cohen et al 1993: https://doi.org/10.1111/1467-8659.1230363
  dem_mask = dem[:DEM] .!= dem[:header]["NODATA_VALUE"]
  maxElevation = maximum(dem[:DEM][dem_mask])
  λᵢ = (maxElevation + threshold - ray[1,3]) / ray[2,3] # number of sampling units until the ray reaches the maxElevation + threshold
  rᵢ = ray[1,:] .+ λᵢ .* ray[2,:]
  zᵢ = getDEMvalue(rᵢ[1], rᵢ[2], dem)
  while rᵢ[3] > zᵢ + threshold
    λᵢ += resolution
    rᵢ = ray[1,:] .+ λᵢ .* ray[2,:]
    zᵢ = getDEMvalue(rᵢ[1], rᵢ[2], dem)
  end
  # refine estimate of λ
  λ = (zᵢ + threshold - ray[1,3]) / ray[2,3]
  λ #sampling units along ray's direction vector until z-coordinate is within the ground threshold
end

function GetPulseDescriptor(PulseDescriptorIndex::Integer, header::PulseWavesHeader)
  vlrRecordIDs = map(x -> x.RecordID, header.VariableLengthRecords)
	header.VariableLengthRecords[only(findall(vlrRecordIDs .== PulseDescriptorIndex + 200000))]
end

function GetLUTs(header::PulseWavesHeader)
  vlrRecordIDs = map(x -> x.RecordID, header.VariableLengthRecords)
	LUTIdcs = findall(300000 .<= vlrRecordIDs .< 300500)
	lowChannelIdx = contains.(map(x -> x.Data.LUT.Description, header.VariableLengthRecords[LUTIdcs]), "low channel")
	highChannelIdx = contains.(map(x -> x.Data.LUT.Description, header.VariableLengthRecords[LUTIdcs]), "high channel")
	@assert(all(xor.(lowChannelIdx, highChannelIdx)))
	LUTDict = Dict(:lowChannel => only(header.VariableLengthRecords[LUTIdcs[lowChannelIdx]]).Data.LUT.Entries, :highChannel => only(header.VariableLengthRecords[LUTIdcs[highChannelIdx]]).Data.LUT.Entries)
	LUTDict
end
	
function WaveToCoverFractionDistribution(wave::Vector{WaveRecord}, groundintersection::Float64, header::PulseWavesHeader, PulseDescriptorIndex::Integer, ground_to_canopy_reflectance_ratio::Float64=1.)
  #get params from all wave segments
  #attenuation correction
  #deconvolution with system pulse
  #mixed model probability density function
  PulseDescriptor = GetPulseDescriptor(PulseDescriptorIndex, header)
  LUTs = GetLUTs(header)
  n_samplings = PulseDescriptor.Data.Composition.NumberOfSamplings	
  amplitudes = Float64[]
  mus = Float64[]
  sigmas = Float64[]
  @assert(wave[1].NumberOfSamples == 28, "Unknown system pulse sampling encountered.")
  system_pulse_params, _, _, _ = IterativelySmoothedGaussianDecomposition(wave[1].Samples; smooth=true, width=1, rel_threshold=.1, abs_threshold=4)
  @assert(size(system_pulse_params)[2] == 1, "More than one pulse detected in system pulse record.")
  system_pulse_sigma = system_pulse_params[3,1]
  k = 1
  for j in 1:n_samplings
    SamplingRecord = PulseDescriptor.Data.Sampling[j]
    n_segments = SamplingRecord.NumberOfSegments
    for i in 1:n_segments#2:length(wave)#3xN (amplitude, mu, sigma)xN
      if SamplingRecord.Type == 2
        params, bias, rmse, actual_smoothing_width = IterativelySmoothedGaussianDecomposition(wave[k].Samples; smooth=true, width=1, rel_threshold=.1, abs_threshold=4)
        #params)[2]
        if SamplingRecord.Channel == 1
	  amplitude = map(x -> dBToPowerRatio(InterpolateLUT(LUTs[:lowChannel],x)), params[1,:])
        elseif SamplingRecord.Channel == 0
	  amplitude = map(x -> dBToPowerRatio(InterpolateLUT(LUTs[:highChannel],x)), params[1,:])
        else
	  throw("Unknown Sampling channel encountered: $(SamplingRecord.Channel) in PulseDescriptor $(PulseDescriptor)")
        end
        append!(amplitudes, amplitude)
        append!(mus, params[2,:] .+ wave[k].DurationFromAnchor)
        append!(sigmas, params[3,:]) # pick up from here!
      end
      k += 1
    end
  end
  #sort hits along the trajectory
  sortidcs = sortperm(mus)
  amplitudes, mus, sigmas = amplitudes[sortidcs], mus[sortidcs], sigmas[sortidcs]
  #deconvolution with system pulse is trivial with Gaussians
  sigmas[-.5 .< sigmas .- system_pulse_sigma .< 0] .= 1.01 * system_pulse_sigma
  sigmas = sqrt.(sigmas.^2 .- system_pulse_sigma^2)
  #calculate total energy per hit
  energy = amplitudes .* sigmas
  # find out which returns are ground hits 
  groundhits = mus .> groundintersection
  # modify the energy of ground hits
  energy[groundhits] .= energy[groundhits] * ground_to_canopy_reflectance_ratio
  # get coverfraction for each individual hit
  coverfraction = energy ./ sum(energy)
  #cum_coverfraction = cumsum(coverfraction)
  # ocassionally the maximum of the cumsum is 1 + 1e-16 or so, causing problems down the line
  #@assert(sum(cum_coverfraction .> 1) <= 1)
  #if any(cum_coverfraction .> 1)
  #cum_coverfraction = cum_coverfraction ./ maximum(cum_coverfraction)
  #end
  # gap fraction is the cumulative energy, normalized to one
  #gapfraction = 1 .- cum_coverfraction
  # second and later returns are corrected by gap fraction up to the return, the first return remains the same
  # but this information is not needed
  #attenuation_corrected_energy = energy ./ vcat([1.],gapfraction[1:end-1])
  # calculate the gap fraction for each hit, i.e. the gap fraction of the remaining visible footprint
  #if length(gapfraction) == 1
  #  gapfraction_i = gapfraction
  #else
  #  gapfraction_i = vcat(gapfraction[1], gapfraction[2:end] ./ gapfraction[1:end-1])
  #end
  #λ = - log.(gapfraction_i)
  #@assert(all(isinf.(λ[groundhits])), "Finite ground hit attenuation encountered.")
  #distributions = [Normal(mu, sigma) for (mu, sigma) in zip(mus, sigmas)]

  cf_distribution = MixtureModel(Normal,[(mu, sigma) for (mu, sigma) in zip(mus, sigmas)], coverfraction)
  # construct mixture model using the canopy hit gap fractions, and an upper boundary at the ground threshold
  #MM = MixtureModel(Normal, [(mu, sigma) for (mu, sigma) in zip(mus[.!groundhits], sigmas[.!groundhits])], 1 .- attenuation_corrected_gapfractions[.!groundhits], upper=groundintersection)
  #distributions[.!groundhits], λ[.!groundhits], groundintersection, cf_distribution#, coverfraction, gapfraction, groundhits, amplitudes, mus, groundintersection, sigmas
  cf_distribution
  # compare wavedorm to deconvoluted foliage profile
  # ground to canopy reflectance method by Armston et al 2013
end

function RayPointIntersection(ray::Matrix{Float64}, xyz::Vector{Float64})
  rxyz = xyz .- ray[1,:]
  rl = LinearAlgebra.normalize(ray[2,:])
  normalized_distance = rl * LinearAlgebra.dot(rl, rxyz)
  r = (normalized_distance ./ ray[2,:])[1]
  closest_point = ray[1,:] .+ normalized_distance
  t = LinearAlgebra.norm(xyz - closest_point)
  r, t
end

function Ray_AABB_Intersection(ray::Matrix{Float64}, aabb::Matrix{Float64})
    direction=ray[2,:]
    dir_fraction = zeros(Float64,3)
    dir_fraction[direction .== 0] .= Inf
    dir_fraction[direction .!= 0] .= 1. ./ direction[direction .!= 0]
    t1 = (aabb[1,1] - ray[1,1]) * dir_fraction[1]
    t2 = (aabb[2,1] - ray[1,1]) * dir_fraction[1]
    t3 = (aabb[1,2] - ray[1,2]) * dir_fraction[2]
    t4 = (aabb[2,2] - ray[1,2]) * dir_fraction[2]
    t5 = (aabb[1,3] - ray[1,3]) * dir_fraction[3]
    t6 = (aabb[2,3] - ray[1,3]) * dir_fraction[3]
    tmin = max(-Inf, min(t1,t2), min(t3,t4), min(t5,t6))
    tmax = min(Inf, max(t1,t2), max(t3,t4), max(t5,t6)) 
                                                        
    if tmax < 0                                         
        return nothing                                  
    end                                                 
    if tmin > tmax                                      
        return nothing                                  
    end                                                 
    t = minimum([x for x in [tmin, tmax] if x >= 0])    
    point = ray[1,:] + (ray[2,:] * t)                   
    return point                                        
end


#function ConstructPointGrid(size, resolution, rotation, origin)

