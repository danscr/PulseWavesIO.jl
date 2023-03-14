struct PulseWavesLookupTableHeader
  Size::UInt32
  Reserved::UInt32
  NumberTables::UInt32
  Description::AbstractString
end

function Base.read(io::IO, ::Type{PulseWavesLookupTableHeader})
  Size = read(io,UInt32)
  Reserved = read(io,UInt32)
  NumberTables = read(io,UInt32)
  Description = readstring(io,64)

  PulseWavesLookupTableHeader(Size,
			      Reserved,
			      NumberTables,
			      Description
			      )
end

struct PulseWavesLookupTableRecord
  Size::UInt32
  Reserved::UInt32
  NumberEntries::UInt32
  UnitOfMeasurement::UInt16
  DataType::UInt8
  Options::UInt8
  Compression::UInt32
  Description::AbstractString
  Entries::Vector{Float32}
end

function Base.read(io::IO, ::Type{PulseWavesLookupTableRecord})
  Size = read(io,UInt32)
  Reserved = read(io, UInt32)
  NumberEntries = read(io,UInt32)
  UnitOfMeasurement = read(io,UInt16)
  DataType = read(io,UInt8)
  Options = read(io,UInt8)
  Compression = read(io,UInt32)
  Description = readstring(io,64)
  Entries = Vector{Float32}(undef,NumberEntries)
  for i ∈ 1:NumberEntries
    Entries[i] = read(io,Float32)
  end

  PulseWavesLookupTableRecord(Size,
			      Reserved,
			      NumberEntries,
			      UnitOfMeasurement,
			      DataType,
			      Options,
			      Compression,
			      Description,
			      Entries
			      )
end

struct PulseWavesLookupTable
  Header::PulseWavesLookupTableHeader
  LUT::PulseWavesLookupTableRecord
end

function Base.read(io::IO, ::Type{PulseWavesLookupTable})
  LUTHeader = read(io, PulseWavesLookupTableHeader)
  LUT = read(io, PulseWavesLookupTableRecord)

  PulseWavesLookupTable(LUTHeader, LUT)
end


struct PulseWavesCompositionRecord
  Size::UInt32
  Reserved::UInt32
  OpticalCenterToAnchorPoint::Int32
  NumberOfExtraWaveBytes::UInt16
  NumberOfSamplings::UInt16
  SampleUnits::Float32
  Compression::UInt32
  ScannerIndex::UInt32
  Description::AbstractString
end

function Base.show(io::IO, cr::PulseWavesCompositionRecord)
  println(io, string("\n\tSize                       = ",cr.Size))
  println(io, string("\tReserved                   = ",cr.Reserved))
  println(io, string("\tOpticalCenterToAnchorPoint = ",cr.OpticalCenterToAnchorPoint))
  println(io, string("\tNumberOfExtraWaveBytes     = ",cr.NumberOfExtraWaveBytes    ))
  println(io, string("\tNumberOfSamplings          = ",cr.NumberOfSamplings         ))
  println(io, string("\tSampleUnits                = ",cr.SampleUnits               ))
  println(io, string("\tCompression                = ",cr.Compression               ))
  println(io, string("\tScannerIndex               = ",cr.ScannerIndex              ))
  println(io, string("\tDescription                = ",cr.Description               ,"\n"))
end

function Base.read(io::IO, ::Type{PulseWavesCompositionRecord})
  Size = read(io, UInt32)
  Reserved = read(io, UInt32)
  OpticalCenterToAnchorPoint = read(io, Int32)
  NumberOfExtraWaveBytes = read(io, UInt16)
  NumberOfSamplings = read(io, UInt16)
  SampleUnits = read(io, Float32)
  Compression = read(io, UInt32)
  ScannerIndex = read(io, UInt32)
  Description = readstring(io, 64)

  @assert(NumberOfExtraWaveBytes == 0, "Extra Wave bytes not yet implemented.")

  PulseWavesCompositionRecord(Size,
		    Reserved,
		    OpticalCenterToAnchorPoint,
		    NumberOfExtraWaveBytes,
		    NumberOfSamplings,
		    SampleUnits,
		    Compression,
		    ScannerIndex,
		    Description
		    )
end

struct PulseWavesSamplingRecord
  Size::UInt32
  Reserved::UInt32 #must be zero
  Type::UInt8 #1=outgoing waveform, 2=returning waveform
  Channel::UInt8 #0 for single sensor, up to h-1 (with h channels)
  Unused::UInt8 #must be zero
  BitsForDurationFromAnchor::UInt8 #"start" of waveform D
  ScaleForDurationFromAnchor::Float32
  OffsetForDurationFromAnchor::Float32 # d = scale * D + offset
  BitsForNumberOfSegments::UInt8
  BitsForNumberOfSamples::UInt8 #non-zero= variable sampling
  NumberOfSegments::UInt16
  NumberOfSamples::UInt32 #non-zero = fixed sampling
  BitsPerSample::UInt16 
  LookupTableIndex::UInt16
  SampleUnits::Float32
  Compression::UInt32
  Description::AbstractString
end

function Base.show(io::IO, sr::PulseWavesSamplingRecord)
  println(io,string("\n\t\tSize                        = ",sr.Size))
  println(io,string("\t\tReserved                    = ",sr.Reserved))
  println(io,string("\t\tType                        = ",sr.Type))
  println(io,string("\t\tChannel                     = ",sr.Channel))
  println(io,string("\t\tUnused                      = ",sr.Unused))
  println(io,string("\t\tBitsForDurationFromAnchor   = ",sr.BitsForDurationFromAnchor))
  println(io,string("\t\tScaleForDurationFromAnchor  = ",sr.ScaleForDurationFromAnchor))
  println(io,string("\t\tOffsetForDurationFromAnchor = ",sr.OffsetForDurationFromAnchor))
  println(io,string("\t\tBitsForNumberOfSegments     = ",sr.BitsForNumberOfSegments))
  println(io,string("\t\tBitsForNumberOfSamples      = ",sr.BitsForNumberOfSamples))
  println(io,string("\t\tNumberOfSegments            = ",sr.NumberOfSegments))
  println(io,string("\t\tNumberOfSamples             = ",sr.NumberOfSamples))
  println(io,string("\t\tBitsPerSample               = ",sr.BitsPerSample))
  println(io,string("\t\tLookupTableIndex            = ",sr.LookupTableIndex))
  println(io,string("\t\tSampleUnits                 = ",sr.SampleUnits))
  println(io,string("\t\tCompression                 = ",sr.Compression))
  println(io,string("\t\tDescription                 = ",sr.Description,"\n"))
end

function Base.read(io::IO, ::Type{PulseWavesSamplingRecord})
  Size = read(io, UInt32)
  Reserved = read(io, UInt32)
  PWType = read(io, UInt8)
  PWChannel = read(io, UInt8)
  Unused = read(io, UInt8)
  BitsForDurationFromAnchor = read(io, UInt8)
  ScaleForDurationFromAnchor = read(io, Float32)
  OffsetForDurationFromAnchor = read(io, Float32)
  BitsForNumberOfSegments = read(io, UInt8)
  BitsForNumberOfSamples = read(io, UInt8)
  NumberOfSegments = read(io, UInt16)
  NumberOfSamples = read(io, UInt32)
  BitsPerSample = read(io, UInt16)
  LookupTableIndex = read(io, UInt16)
  SampleUnits = read(io, Float32)
  Compression = read(io, UInt32)
  Description = readstring(io, 64)
  @assert(BitsForDurationFromAnchor == 32, "Bits for duration from anchor != 32, this has not been implemented yet.")
  @assert(NumberOfSamples == 0, "Fixed number of samples is not implemented.")
  @assert(BitsPerSample == 8, "Only 8 bits per sample are implemented.")
  @assert(Compression == 0, "Compression in Wave sampling not implemented.")

  PulseWavesSamplingRecord(Size                        ,
			    Reserved                   ,
  	      		    PWType                     ,
  	      		    PWChannel                  ,
  	      		    Unused                     ,
  	      		    BitsForDurationFromAnchor  ,
  	      		    ScaleForDurationFromAnchor ,
  	      		    OffsetForDurationFromAnchor,
  	      		    BitsForNumberOfSegments    ,
  	      		    BitsForNumberOfSamples     ,
  	      		    NumberOfSegments           ,
  	      		    NumberOfSamples            ,
  	      		    BitsPerSample              ,
  	      		    LookupTableIndex           ,
  	      		    SampleUnits                ,
  	      		    Compression                ,
  	      		    Description
			    )
end

struct PulseWavesPulseDescriptor
  Composition::PulseWavesCompositionRecord
  Sampling::Vector{PulseWavesSamplingRecord}
end

function Base.read(io::IO, ::Type{PulseWavesPulseDescriptor})
  Composition = read(io, PulseWavesCompositionRecord)
  Sampling = Vector{PulseWavesSamplingRecord}(undef,Composition.NumberOfSamplings)
  for i ∈ 1:Composition.NumberOfSamplings
      Sampling[i] = read(io, PulseWavesSamplingRecord)
  end
  
  PulseWavesPulseDescriptor(Composition, Sampling)
end

struct PulseWavesVariableLengthRecord
    UserID      ::AbstractString
    RecordID    ::Int32
    Reserved    ::UInt32
    RecordLengthAfterHeader ::Int64
    Description ::AbstractString
    Data#::Union{PulseWavesPulseDescriptor,PulseWavesLookupTable,Vector{UInt8}}
end

function Base.show(io::IO, vlr::PulseWavesVariableLengthRecord)
  println(io,string("\nUserID                   = ",vlr.UserID      ))
  println(io,string("RecordID                 = ",vlr.RecordID    ))
  println(io,string("Reserved                 = ",vlr.Reserved    ))
  println(io,string("RecordLengthAfterHeader  = ",vlr.RecordLengthAfterHeader ))
  println(io,string("Description              = ",vlr.Description ))
  println(io,string("Data                     = ",vlr.Data,"\n"))
end

function Base.read(io::IO, ::Type{PulseWavesVariableLengthRecord})
    UserID = readstring(io, 16)
    RecordID = read(io, Int32)
    Reserved = read(io, UInt32)
    RecordLengthAfterHeader = read(io, UInt64)
    Description = readstring(io, 64)
    Data = read(io, RecordLengthAfterHeader)#ReadVLRData(io, RecordID, RecordLengthAfterHeader)
    #if RecordID == 100001
    #  throw("Not implemented")
    if 200000 <= RecordID < 200255
      Data = read(IOBuffer(Data), PulseWavesPulseDescriptor) #TODO: here be bugs. If the IOBuffer isn't used but instead the IO stream directrly is read for a PulseDescriptor, it throws an OutOfMemory error
    elseif 300000 <= RecordID < 300255
      Data = read(IOBuffer(Data), PulseWavesLookupTable)
    end
    #elseif RecordID == 34735
    #  throw("Not implemented")
    #elseif RecordID == 34736
    #  throw("Not implemented")
    #elseif RecordID == 34737
    #  throw("Not implemented")
    #end

    PulseWavesVariableLengthRecord(
                                   UserID,
                                   RecordID,
                                   Reserved,
                                   RecordLengthAfterHeader,
                                   Description,
                                   Data
                                  )
end

function Base.write(io::IO, VLR::PulseWavesVariableLengthRecord)
    writestring(io, VLR.UserID, 16)
    write(io, VLR.RecordID)
    write(io, VLR.Reserved)
    @assert(sizeof(VLR.Data) == VLR.RecordLengthAfterHeader)
    write(io, VLR.RecordLengthAfterHeader)
    write(io, VLR.Description)
    write(io, VLR.Data)
    nothing
end
