using JLD2, CSV, DataFrames, Distributions, FastGaussQuadrature, LinearAlgebra, ProgressMeter#, LsqFit, ImageFiltering, Peaks 
import Distributed.pmap
#using Revise
using PulseWavesIO

SU_to_m = .14985490630017206 #Sampling units to meters (approximately one lightnanosecond, measured from 100 pulses in the data)


plots = CSV.read("scratch/forest_structure.csv", DataFrame)
plotnumber = parse(Int64,ARGS[1])
resolution = parse(Float64,ARGS[2])
flightlinenumber = parse(Int64,ARGS[3])
outpath = "scratch/$(plots[plotnumber,:plot_ID])_res$(resolution)m_voxelgrid_flightline$(flightlinenumber).jld2"

"""
Returns two scalars at entry and exit to an axis-aligned bounding box along a ray.
"""
function Ray_AABB_intersection(ray, aabb)
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
return tmin, tmax
end


function InterceptanceTrilinearWeights(cfdist::MixtureModel, nodes::Vector{Float64}, weights::Vector{Float64}, ray::Matrix{Float64}, xyzi::Vector{Float64}, resolution::Float64)
ray_aabb_intersection = Ray_AABB_intersection(ray, hcat(xyzi .- resolution/2, xyzi .+ resolution/2)') #voxel grid -> half resolution compared to grid spacing
if !isnothing(ray_aabb_intersection)
a, b = ray_aabb_intersection
cf_0 = cdf(cfdist, a)
distance_to_xyz(x) = abs.(ray[1,:] .+ ray[2,:] .* x .- xyzi)
w(x) = any(distance_to_xyz(x) .> resolution ./ 2) ? 0. : 1.  ### now with voxel weights, the following is inverse distance weight: #prod(1 .- distance_to_xyz(x) ./ resolution)
integral_spacing = 1e-3
#delta_cf = pdf.(cfdist, a:integral_spacing:b) * integral_spacing .* w.(a:integral_spacing:b) * SU_to_m
delta_cf = pdf.(cfdist, a:integral_spacing:b) .* w.(a:integral_spacing:b) .* integral_spacing
#delta_cf = GaussLegendreIntegral(x -> pdf(cfdist,x) * w(x), a, b, nodes, weights)
#W = GaussLegendreIntegral(x -> w(x), a,b,nodes,weights)
W = sum(w.(a:integral_spacing:b) .* integral_spacing * SU_to_m)
attenuation = sum(delta_cf) / (1-cf_0) / W
return nothing, 1-cf_0, attenuation, W
else
return nothing
end
end

function WaveInformation_trilinear(ray::Matrix{Float64}, xyzi::Vector{Float64}, resolution::Float64, wave::Vector{PulseWavesIO.WaveRecord}, groundintersection::Float64, header::PulseWavesIO.PulseWavesHeader, pdi::Integer, nodes::Vector{Float64}, weights::Vector{Float64})
try
local cfdist = WaveToCoverFractionDistribution(wave, groundintersection, header, pdi)
return InterceptanceTrilinearWeights(cfdist, nodes, weights, ray, xyzi, resolution)
catch e
#throw("uhetonas")
if isa(e, String)
if occursin("Iterative smoothing", e)
return nothing
else
println(wave)
#rethrow(e)
end
#elseif isa(e, AssertionError)
#  return nothing
elseif isa(e,AssertionError)
return nothing#string(e)
else
println(e)
end
end
end


#distributions, lambda, gi, cfdist = WaveToCoverFractionDistribution(waves[7], rays[7,:,:], dem, header, pdi[7])
function GaussLegendreIntegral(f, a, b, nodes, weights)
dot(weights, f.((a+b)/2 .+ (b-a) ./ 2 .* nodes) * (b-a)/2)
end

if isfile(outpath)
println("Output file already exists, skipping...")
else
jldopen("scratch/pulsewavesfile.jld2") do io
#global header, pulses, rays, waves, pdi, dem = (io[i] for i in keys(io))
global dem = io["dem"]
end
# construct sampling point grid
#plotnumber = 11
xy = Vector(plots[plotnumber,[:x_utm,:y_utm]]);
z0 = PulseWavesIO.getDEMvalue(xy[1],xy[2],dem);
#resolution = .2;
println("Processing plot $(plotnumber) ($(plots[plotnumber,:plot_ID])) at resolution $(resolution).")
aabb = [xy[1]-12.5 xy[2]-12.5 z0-20; xy[1]+12.5 xy[2]+12.5 z0+50];
xyz = [[x,y,z] for x in xy[1]-12.5:resolution:xy[1]+12.5 for y in xy[2]-12.5:resolution:xy[2]+12.5 for z in z0-20:resolution:z0+50];
xyz = Matrix(hcat(xyz...)');

path = "scratch/pulse_waves/" 
#path = "dev/" 
filenames = path .* filter(x -> startswith(x,"PulseWaves") & endswith(x,"pls"),readdir(path));
global gauss_nodes, gauss_weights = gausslegendre(100);

#flightlinenumber = 1
fl_atts = []
fl_n = []
fl_w = []

for i in 1
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
beam_subset = findall(in_bbox);
groundIntersections = pmap(i -> PulseWavesIO.RayGroundIntersection(rays[i,:,:], dem), beam_subset);

#println("Extracting waveform information...")
#old weights##wave_info = map(i -> WaveInformation(waves[beam_subset[i]], groundIntersections[i], header,pdi[beam_subset[i]], gauss_nodes, gauss_weights), 1:length(beam_subset));
n_pts = size(xyz)[1];
atts = Vector{Float64}(undef, n_pts);
w_beams = Vector{Float64}(undef, n_pts);
n_beams = Vector{Int64}(undef, n_pts);
#s = 2.;
global R = resolution*sqrt(3.1);
#progress = Progress(n_pts, dt=1.0)
#println("Beginnig spatial point sampling...")
Threads.@threads for j in 1:n_pts
#for j in 1:n_pts
#@time begin
rts = map(i -> PulseWavesIO.RayPointIntersection(rays[i,:,:], xyz[j,:]), beam_subset);
### weighting scheme here

###
mask_ground = map(i -> rts[i][1] <= groundIntersections[i], 1:length(rts));
#all(mask_ground);
ts = map(i -> rts[i][2], 1:length(rts));
mask_lateral = ts .<= R;#map(i -> rts[i][2] <= R, 1:length(rts));
nearby = findall(mask_ground .& mask_lateral);
if length(nearby) > 0
freepaths_gaps0 = map(i -> WaveInformation_trilinear(rays[beam_subset[i],:,:],xyz[j,:], resolution, waves[beam_subset[i]], groundIntersections[i], header, pdi[beam_subset[i]], gauss_nodes, gauss_weights), nearby)
#nearby = nearby[.!isnothing(freepaths_gaps0)]
filter!(x -> !isnothing(x), freepaths_gaps0)
if length(freepaths_gaps0) > 0
#waveform_processing_converged = freepaths_gaps0 .!= nothing
#freepaths_gaps0 = freepaths_gaps0[waveform_processing_converged]
#freepaths_gaps0 = map(i -> wave_info[i], 1:length(nearby));
#freepaths = map(i -> freepaths_gaps0[i][1], 1:length(nearby));
gaps0 = map(i -> i[2], freepaths_gaps0);
attenuation = map(i -> i[3], freepaths_gaps0);
weights = map(i -> i[4],freepaths_gaps0);
#t_weight = lateral_weight.(ts[nearby]);#map(i -> lateral_weight(rts[i][2]), nearby);
unoccluded = gaps0 .> 0.01
if any(unoccluded)
weight = weights[unoccluded] .* gaps0[unoccluded]
weight = weight ./ sum(weight)
attenuation = sum(attenuation[unoccluded] .* weight)
atts[j] =  attenuation
#println(sum(unoccluded))
n_beams[j] =  sum(unoccluded)
w_beams[j] = sum(gaps0[unoccluded])
else
atts[j] = NaN
n_beams[j] = 0
w_beams[j] = 0
end
else
atts[j] = NaN
n_beams[j] = 0
w_beams[j] = 0
end
else
atts[j] = NaN
n_beams[j] = 0
w_beams[j] = 0
end
#next!(progress)
end
append!(fl_atts, [atts])
append!(fl_n,[n_beams])
append!(fl_w,[w_beams])
#alert("Flight line $(flightlinenumber) is processed.")
println("Flight line $(flightlinenumber) is processed. LAI = $(sum(atts[.!isnan.(atts) .& .!isinf.(atts)]) * resolution^3 / .5 / 25^2)")
end
end

jldsave(outpath ; fl_atts, fl_n, fl_w)
end #end of else from if isfile(outpath)
println("Workflow complete.")
