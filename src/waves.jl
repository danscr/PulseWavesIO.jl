struct WavesHeader
    fileSignature ::AbstractString
    Compression ::UInt32
    Reserved    ::AbstractString
end

function Base.read(io::IO, ::Type{WavesHeader})
    fileSignature = readstring(io,16)
    @assert(fileSignature == "PulseWavesWaves", "Waves file does not contain valide file signature.")
    Compression = read(io,UInt32)
    @assert(Compression == 0, "Compressed waves not implemented.")
    Reserved = readstring(io,40)

    WavesHeader(
		fileSignature,
		Compression,
		Reserved
		)
end

function Base.write(io::IO, wv_header::WavesHeader)
  write(io, wv_header.fileSignature)
  write(io, wv_header.Compression)
  write(io, wv_header.Reserved)
end

"""
Wave structure is highly flexible, only one instance is implemented for now.
"""
mutable struct WaveRecord
  DurationFromAnchor::Union{UInt32,Float32}
  NumberOfSamples::UInt16
  Samples::Union{Vector{UInt8},Vector{Float32}}
end

function Base.show(io::IO, w::WaveRecord)
  println(io,string("\nDurationFromAnchor = ",w.DurationFromAnchor))
  println(io,string("NumberOfSamples    = ",w.NumberOfSamples))
  println(io,string("Samples            = ",w.Samples))
end

function readWave(io::IO, p::PulseRecord, h::PulseWavesHeader, convert_to_power_ratio=true)
  vlr = only(filter(x -> x.RecordID - 200000 == p.PulseDescriptorIndex, h.VariableLengthRecords))
  n_samplings = vlr.Data.Composition.NumberOfSamplings
  seek(io, p.OffsetToWaves)
  dat = []
  #SamplingTypes = []
  for sampling in 1:n_samplings
    samplingRecord = vlr.Data.Sampling[sampling]
    n_segments = samplingRecord.NumberOfSegments
    LUT_idx = samplingRecord.LookupTableIndex
    LUT = only(filter(x -> x.RecordID - 300000 == LUT_idx, h.VariableLengthRecords)).Data.LUT
    @assert(LUT.UnitOfMeasurement == 1, "Unexpected lookup table encountered. Lookup tables are only used for intensity correction (Unit of measurement == 1), but other was encountered (undefined or range correction)")
    for segment in 1:n_segments
      wv = read(io, WaveRecord)
      if convert_to_power_ratio
	wv.Samples = dBToPowerRatio.(LUT.Entries[wv.Samples .+ 1]) #LUT is in dB, we want the power ratio, i.e. the intensity
      end
      #however, scaling in this way introduces non-linearity and the peaks may be exaggerated when analysing the entire waveform
      wv.DurationFromAnchor = Float32(samplingRecord.OffsetForDurationFromAnchor + samplingRecord.ScaleForDurationFromAnchor * wv.DurationFromAnchor)
      append!(dat,[wv])
      #append!(SamplingTypes, samplingRecord.Type)
    end
  end
  Vector{WaveRecord}(dat)#, p.PulseDescriptorIndex# Vector{UInt8}(SamplingTypes)
end

function Base.read(io::IO, ::Type{WaveRecord})
  DurationFromAnchor = read(io, UInt32)
  NumberOfSamples = read(io, UInt16)
  Samples = Vector{UInt8}(undef, NumberOfSamples)
  for i ∈ 1:NumberOfSamples
    Samples[i] = read(io, UInt8)
  end

  WaveRecord(
	      DurationFromAnchor,
	      NumberOfSamples,
	      Samples
	     )
end


function WaveformToGaussianMixture(wave::Vector{WaveRecord}, PulseDescriptorIndex::Int, header::PulseWavesHeader)
  vlr = only(filter(x -> x.RecordID == PulseDescriptorIndex + 200000, header.VariableLengthRecords))
  n_samplings = vlr.Data.Composition.NumberOfSamplings
  luts = filter(x -> 300000 <= x.RecordID < 300100, header.VariableLengthRecords)
  lowIdx = map(x -> x.RecordID == 300001, luts)
  luts = map(x -> x.Data.LUT, luts)
  amplitudes = Float64[]
  mus = Float64[]
  sigmas = Float64[]
  
  k = 0
  for i in 1:n_samplings
    SamplingRecord = vlr.Data.Sampling[i]
    n_segments = SamplingRecord.NumberOfSegments
    for j in 1:n_segments
      k += 1
      if SamplingRecord.Type == 2
	dc = rcopy(wflidar.decom_adaptive(wave[k].Samples, smooth=true, thres=0.22, width=3))
	n_peaks = size(dc[3])[1]
	amplitudes = dc[3][:,2]
	if SamplingRecord.Channel == 1
	  lutIdx = lowIdx
	elseif SamplingRecord.Channel == 0
	  lutIdx = .!loxIdx
	else
	  throw("Unknown sampling channel encountered: $(SamplingRecord.Channel)")
	end
	lut = dBToPowerRatio.(only(luts[lutIdx]).Data.LUT.Entries)
	append!(amplitudes, map(x -> InterpolateLUT(lut, x), amplitudes))
	append!(mus, dc[3][:,3])
	append!(sigmas, dc[3][:,4])
      end
    end
  end
  #create mixture model
end

function mergeWaveSegments(wave::Vector{WaveRecord}, vlr::PulseWavesVariableLengthRecord)
  n_samplings = vlr.Data.Composition.NumberOfSamplings
  MergedWave = zeros(UInt8, 0)
  LowPowerChannel = zeros(Bool,0)
  k = 0
  for i in 1:n_samplings
    SamplingRecord = vlr.Data.Sampling[i]
    n_segments = SamplingRecord.NumberOfSegments
    for j in 1:n_segments
      k += 1
      if SamplingRecord.Type == 2
	if j != n_segments
	  dt = Int(round(wave[k+1].DurationFromAnchor - wave[k].DurationFromAnchor))
	  segment = zeros(UInt8, dt-1)
	  segment[1:wave[k].NumberOfSamples] .= wave[k].Samples
	  MergedWave = vcat(MergedWave, segment)
	else
	  MergedWave = vcat(MergedWave, wave[k].Samples)
	end
      end
    end
    if SamplingRecord.Channel == 1
      LowPowerChannel = vcat(LowPowerChannel, ones(Bool, length(MergedWave)))
    elseif SamplingRecord.Channel == 0
      LowPowerChannel = vcat(LowPowerChannel, zeros(Bool, length(MergedWave)))
    end
  end
  MergedWave, LowPowerChannel
  #return one long vector where the gaps between segments are zeros
end

function mergeWaves(pulses::Vector{PulseRecord}, waves::Vector{WaveRecord}, header::PulseWavesHeader)
  vlrRecordIDs = map(x -> x.RecordID, header.VariableLengthRecords)
  PulseDescriptorIndices = map(x -> x.PulseDescriptorIndex, pulses)
  vlrIdcs = map(x -> x[1], findall(vlrrids .== PulseDescriptorIndices' .+ 200000))
end
