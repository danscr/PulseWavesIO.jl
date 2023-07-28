using JLD2, GLMakie, CSV, DataFrames, Distributions, FastGaussQuadrature, LinearAlgebra, ProgressMeter#, LsqFit, ImageFiltering, Peaks 
using Alert
import Distributed.pmap
using Revise
using PulseWavesIO

jldopen("pulsewavesfile.jld2") do io
  #global header, pulses, rays, waves, pdi, dem = (io[i] for i in keys(io))
  global dem = io["dem"]
end


#distributions, lambda, gi, cfdist = WaveToCoverFractionDistribution(waves[7], rays[7,:,:], dem, header, pdi[7])
function GaussLegendreIntegral(f, a, b, nodes, weights)
  dot(weights, f.((a+b)/2 .+ (b-a) ./ 2 .* nodes) * (b-a)/2)
end


function FreePath(cfdist, r, s, nodes, weights)
  a, b = r-s/2, r+s/2
  cf_0 = cdf(cfdist, a)
  freepath = GaussLegendreIntegral(x -> 1 - cdf(cfdist, x), a, b, nodes, weights) / (1 - cf_0) 
  δ_cf = cdf(cfdist, b) - cf_0
  #println(cf_0, cdf(cfdist,b))
  freepath, 1 - cf_0, δ_cf
end

function WaveInformation(wave::Vector{PulseWavesIO.WaveRecord}, groundintersection::Float64, header::PulseWavesIO.PulseWavesHeader, pdi::Integer, nodes::Vector{Float64}, weights::Vector{Float64})
  try
    local cfdist = WaveToCoverFractionDistribution(wave, groundintersection, header, pdi)
    return f(r,s=5) = FreePath(cfdist, r, s, nodes, weights)
  catch e
    #throw("uhetonas")
    if isa(e, String)
      if occursin("Iterative smoothing", e)
        return nothing
      else
	println(wave)
	rethrow(e)
      end
    #elseif isa(e, AssertionError)
    #  return nothing
  elseif isa(e,AssertionError)
    return nothing#string(e)
  else
      rethrow(e)
    end
  end
end
lateral_weight(t::Float64,beam_footprint_size::Float64=0.15) = pdf.(Normal(0,1), t / beam_footprint_size)

# construct sampling point grid
plots = CSV.read("C:/Users\\schraid1\\Downloads\\forest_structure.csv", DataFrame)
plotnumber = 11
xy = Vector(plots[plotnumber,[:x_utm,:y_utm]]);
z0 = PulseWavesIO.getDEMvalue(xy[1],xy[2],dem);
resolution = .5;
aabb = [xy[1]-12.5 xy[2]-12.5 z0-20; xy[1]+12.5 xy[2]+12.5 z0+50];
xyz = [[x,y,z] for x in xy[1]-12.5:resolution:xy[1]+12.5 for y in xy[2]-12.5:resolution:xy[2]+12.5 for z in z0-20:resolution:z0+50];
xyz = Matrix(hcat(xyz...)');

path = "E:/pulse_waves/" 
#path = "dev/" 
filenames = path .* filter(x -> endswith(x,"pls"),readdir(path));
nodes, weights = gausslegendre(100);

#flightlinenumber = 1
fl_atts = []
fl_n = []


@time begin
for flightlinenumber in 5:length(filenames)
  @time begin
  println("Processing flight line $flightlinenumber.")
  #println("Reading PulseWaves data...")
  local header, pulses, wv_header, waves, avlr
  try
  header, pulses, wv_header, waves, avlr = load(filenames[flightlinenumber]);
  catch e
    println(e)
    println("Skipping flight line.")
    continue
  end
  rays = PulsesToRays(pulses, header);
  pdi = map(x -> x.PulseDescriptorIndex, pulses);
  # filter rays that cross the bbox
  #println("Calculating plot and pulse geometries...")
  in_bbox = .!isnothing.(map(i -> PulseWavesIO.Ray_AABB_Intersection(rays[i,:,:],aabb), 1:size(rays)[1]));
  if sum(in_bbox) <= 0
    println("No rays within the bounding box, skipping flight line.")
    continue
  end
  #@assert(sum(in_bbox) > 0, "No rays within the bounding box.")
  ### the code below here needs to be expanded for all grid points
  # for all ray-point combinations, find the distance to the closest point on the ray
  beam_subset = findall(in_bbox);
  groundIntersections = pmap(i -> PulseWavesIO.RayGroundIntersection(rays[i,:,:], dem), beam_subset);
  
  #println("Extracting waveform information...")
  wave_info = pmap(i -> WaveInformation(waves[beam_subset[i]], groundIntersections[i], header,pdi[beam_subset[i]], nodes, weights), 1:length(beam_subset));
  n_pts = size(xyz)[1];
  atts = Vector{Float64}(undef, n_pts);
  n_beams = Vector{Int64}(undef, n_pts);
  #s = 2.;
  R = 5.;
  progress = Progress(n_pts, dt=1.0)
  #println("Beginnig spatial point sampling...")
  Threads.@threads for j in 1:n_pts
  #for j in 1:n_pts
    #@time begin
    rts = map(i -> PulseWavesIO.RayPointIntersection(rays[i,:,:], xyz[j,:]), beam_subset);
    mask_ground = map(i -> rts[i][1] <= groundIntersections[i], 1:length(rts));
    #all(mask_ground);
    ts = map(i -> rts[i][2], 1:length(rts));
    mask_lateral = ts .< R;#map(i -> rts[i][2] <= R, 1:length(rts));
    nearby = findall(mask_ground .& mask_lateral .& .!isnothing.(wave_info));
    if length(nearby) > 0
      s = max.(2.,ts)#2*sqrt.(R^2 .- ts[nearby] .^ 2)
      freepaths_gaps0 = map(i -> wave_info[nearby[i]](rts[nearby[i]][1],s[i]), 1:length(nearby));
      freepaths = map(i -> freepaths_gaps0[i][1], 1:length(nearby));
      gaps0 = map(i -> freepaths_gaps0[i][2], 1:length(nearby));
      interceptance = map(i -> freepaths_gaps0[i][3], 1:length(nearby));
      interceptance ./= gaps0
      t_weight = lateral_weight.(ts[nearby]);#map(i -> lateral_weight(rts[i][2]), nearby);
      unoccluded = gaps0 .> 0.01
      weight = t_weight[unoccluded] .* gaps0[unoccluded]#ones(length(gaps0[unoccluded]))#t_weight .* gaps0
      weight ./= sum(weight)
      #gaps0_w = ones(length(gaps0))#gaps0 ./ sum(gaps0)
      #t_weight = ones(length(gaps0))#t_weight ./ sum(t_weight)
      #sphere_radius = sqrt.(ts[nearby] .^2 .+ (s/2)^2)
      #sphere_volume = 4/3*π* (sphere_radius .^3)
      attenuation = sum(interceptance[unoccluded] ./ freepaths[unoccluded] ./ 0.14985 .* weight)# .* 2 .* sphere_volume) / sum(sphere_volume)# / sum(t_weight .* gaps0_w)# / mean(freepaths ./ s)
      atts[j] =  attenuation
      n_beams[j] =  length(gaps0)
    else
      atts[j] = NaN
      n_beams[j] = 0
    end
    next!(progress)
  end
  append!(fl_atts, [atts])
  append!(fl_n,[n_beams])
  alert("Flight line $(flightlinenumber) is processed.")
  println("Flight line $(flightlinenumber) is processed. LAI = $(sum(atts[.!isnan.(atts) .& .!isinf.(atts)]) * resolution^3 / .5 / 25^2)")
end
end
end

jldsave("dev/working_results_res50cm_variableS_weighted.jld2"; fl_atts, fl_n)

#function gpu_processing!(atts, n_beams, xyz, rays, groundIntersections, beam_subset, header, pdi)
#  index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#  stride = gridDim().x * blockDim().x
#  waveinfo = WaveInformation(waves[beam_subset[i]], groundIntersections[i], header,pdi[beam_subset[i]])
#  for j in index:stride:n_pts
#    #@time begin
#    rts = map(i -> PulseWavesIO.RayPointIntersection(rays[i,:,:], xyz[j,:]), beam_subset);
#    mask_ground = map(i -> rts[i][1] <= groundIntersections[i], 1:length(rts));
#    all(mask_ground);
#    mask_lateral = map(i -> rts[i][2] <= 1, 1:length(rts));
#    nearby = findall((mask_ground .& mask_lateral) .& .!isnothing.(wave_info));
#    if length(nearby) > 0
#      freepaths_gaps0 = map(i -> wave_info[i](rts[i][1],s), nearby);
#      freepaths = map(i -> freepaths_gaps0[i][1], 1:length(nearby));
#      gaps0 = map(i -> freepaths_gaps0[i][2], 1:length(nearby));
#      interceptance = map(i -> freepaths_gaps0[i][3], 1:length(nearby));
#      t_weight = map(i -> lateral_weight(rts[i][2]), nearby);
#      gaps0_w = gaps0 ./ sum(gaps0)
#      t_weight = t_weight ./ sum(t_weight)
#      attenuation = sum(interceptance .* t_weight .* gaps0_w) / sum(freepaths .* t_weight .* gaps0_w)
#      atts[j] =  attenuation
#      n_beams[j] =  length(gaps0)
#    else
#      atts[j] = NaN
#      n_beams[j] = 0
#    end
#    return nothing
#  end
#end
#
#function WaveInfoToAttenuation(wave_info_i, rts_i,s)
#  freepath, gap0, interceptance = wave_info_i(rts_i[1],s)
#  t_weight = lateral_weight(rts_i[2])
#  gap0 ./= sum(gap0)
#  t_weight ./= sum(t_weight)
#  sum(interceptance .* t_weight .* gap0) / sum(freepaths .* t_weight .* gap0)
#end
#
#@showprogress for j in 1:n_pts
#  #@time begin
#  @time rts = map(i -> PulseWavesIO.RayPointIntersection(rays[i,:,:], xyz[j,:]), beam_subset);
#  @time mask_ground = map(i -> rts[i][1] <= groundIntersections[i], 1:length(rts));
#  @time mask_lateral = map(i -> rts[i][2] <= 1, 1:length(rts));
#  @time nearby = findall((mask_ground .& mask_lateral) .& .!isnothing.(wave_info));
#  if length(nearby) > 0
#    @time freepaths_gaps0 = map(i -> wave_info[i](rts[i][1],s), nearby);
#    @time freepaths = map(i -> freepaths_gaps0[i][1], 1:length(nearby));
#    @time gaps0 = map(i -> freepaths_gaps0[i][2], 1:length(nearby));
#    @time interceptance = map(i -> freepaths_gaps0[i][3], 1:length(nearby));
#    @time t_weight = map(i -> lateral_weight(rts[i][2]), nearby);
#    @time gaps0_w = gaps0 ./ sum(gaps0);
#    @time t_weight = t_weight ./ sum(t_weight);
#    @time attenuation = sum(interceptance .* t_weight .* gaps0_w) / sum(freepaths .* t_weight .* gaps0_w);
#    atts[j] =  attenuation
#    n_beams[j] =  length(gaps0)
#  else
#    atts[j] = NaN
#    n_beams[j] = 0
#  end
#end


