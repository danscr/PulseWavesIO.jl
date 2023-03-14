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

function readWave(io::IO, p::PulseRecord, h::PulseWavesHeader)
  vlr = only(filter(x -> x.RecordID - 200000 == p.PulseDescriptorIndex, h.VariableLengthRecords))
  n_samplings = vlr.Data.Composition.NumberOfSamplings
  seek(io, p.OffsetToWaves)
  dat = []
  SamplingTypes = []
  for sampling in 1:n_samplings
    samplingRecord = vlr.Data.Sampling[sampling]
    n_segments = samplingRecord.NumberOfSegments
    LUT_idx = samplingRecord.LookupTableIndex
    LUT = only(filter(x -> x.RecordID - 300000 == LUT_idx, h.VariableLengthRecords)).Data.LUT
    @assert(LUT.UnitOfMeasurement == 1, "Unexpected lookup table encountered. Lookup tables are only used for intensity correction (Unit of measurement == 1), but other was encountered (undefined or range correction)")
    for segment in 1:n_segments
      wv = read(io, WaveRecord)
      wv.Samples = LUT.Entries[wv.Samples .+ 1]
      wv.DurationFromAnchor = Float32(samplingRecord.OffsetForDurationFromAnchor + samplingRecord.ScaleForDurationFromAnchor * wv.DurationFromAnchor)
      append!(dat,[wv])
      append!(SamplingTypes, samplingRecord.Type)
    end
  end
  Vector{WaveRecord}(dat), Vector{UInt8}(SamplingTypes)
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
