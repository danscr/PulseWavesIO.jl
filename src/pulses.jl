struct PulseRecord
    GPSTimestampT::Int64
    OffsetToWaves::Int64
    AnchorX::Int32
    AnchorY::Int32
    AnchorZ::Int32
    TargetX::Int32
    TargetY::Int32
    TargetZ::Int32
    FirstReturningSample::Int16
    LastReturningSample::Int16
    PulseDescriptorIndex::Int8
    OtherData::Int8
    #Reserved::UInt8
    #EdgeOfScanLine::Bool
    #ScanDirection::Bool
    #MirrorFacet::Int8
    Intensity::UInt8
    Classification::UInt8
end

function ReadAncillaryPulseData(dat::Int8)
  data = bitstring(dat)
  @assert(length(data) == 8)
  Reserved = parse(UInt8, data[1:4],base=2)
  EdgeOfScanLine = parse(Bool, data[5])
  ScanDirection = parse(Bool, data[6])
  MirrorFacet = parse(Int8, data[7:8],base=2)
  Reserved, EdgeOfScanLine, ScanDirection, MirrorFacet
end

function WriteAncillaryPulseData(Reserved::UInt8, EdgeOfScanLine::Bool, ScanDirection::Bool, MirrorFacet::Int8)
  @assert(Reserved == 0, "Error; Pulse encountered where Reserved byte is non-zero")
  Reserved_bs = bitstring(Reserved)[5:end]
  EdgeOfScanLine_bs = bitstring(EdgeOfScanLine)[end]
  ScanDirection_bs = bitstring(ScanDirection)[end]
  @assert(MirrorFacet <= 3, "Overflow error with Mirror Facet")
  MirrorFacet_bs = bitstring(MirrorFacet)[7:end]
  outString = Reserved_bs * EdgeOfScanLine_bs * ScanDirection_bs * MirrorFacet_bs
  parse(Int8, outString)
end

function Base.read(io::IO, ::Type{PulseRecord})
  GPSTimestampT = read(io,Int64)
  OffsetToWaves = read(io,Int64)
  AnchorX = read(io,Int32)
  AnchorY = read(io,Int32)
  AnchorZ = read(io,Int32)
  TargetX = read(io,Int32)
  TargetY = read(io,Int32)
  TargetZ = read(io,Int32)
  FirstReturningSample = read(io,Int16)
  LastReturningSample = read(io,Int16)
  PulseDescriptorIndex = read(io,Int8)
  OtherData = read(io,Int8)
  #Reserved, EdgeOfScanLine, ScanDirection, MirrorFacet = ReadAncillaryPulseData(read(io,Int8))
  Intensity = read(io,UInt8)
  Classification = read(io,UInt8)

  PulseRecord(
    GPSTimestampT,
    OffsetToWaves,
    AnchorX,
    AnchorY,
    AnchorZ,
    TargetX,
    TargetY,
    TargetZ,
    FirstReturningSample,
    LastReturningSample,
    PulseDescriptorIndex,
    OtherData,
    #Reserved,
    #EdgeOfScanLine,
    #ScanDirection,
    #MirrorFacet,
    Intensity,
    Classification)
end

function Base.write(io::IO, p::PulseRecord)
  write(io, p.GPSTimestampT)
  write(io, p.OffsetToWaves)
  write(io, p.AnchorX)
  write(io, p.AnchorY)
  write(io, p.AnchorZ)
  write(io, p.TargetX)
  write(io, p.TargetY)
  write(io, p.TargetZ)
  write(io, p.FirstReturningSample)
  write(io, p.LastReturningSample)
  write(io, p.PulseDescriptorIndex)
  write(io, p.OtherData)#WriteAncillaryPulseData(p.Reserved, p.EdgeOfScanLine, p.ScanDirection, p.MirrorFacet))
  write(io, p.Intensity)
  write(io, p.Classification)
end

function Base.show(io::IO, p::PulseRecord)
  println(io,string("GPSTimestampT        = ",p.GPSTimestampT))
  println(io,string("OffsetToWaves        = ",p.OffsetToWaves))
  println(io,string("AnchorX              = ",p.AnchorX))
  println(io,string("AnchorY              = ",p.AnchorY))
  println(io,string("AnchorZ              = ",p.AnchorZ))
  println(io,string("TargetX              = ",p.TargetX))
  println(io,string("TargetY              = ",p.TargetY))
  println(io,string("TargetZ              = ",p.TargetZ))
  println(io,string("FirstReturningSample = ",p.FirstReturningSample))
  println(io,string("LastReturningSample  = ",p.LastReturningSample))
  println(io,string("PulseDescriptorIndex = ",p.PulseDescriptorIndex))
  println(io,string("Other Data		  = ",bitstring(p.OtherData)))
  #println(io,string("Reserved             = ",p.Reserved))
  #println(io,string("EdgeOfScanline       = ",p.EdgeOfScanLine))
  #println(io,string("ScanDirection        = ",p.ScanDirection))
  #println(io,string("MirrorFacet          = ",p.MirrorFacet))
  println(io,string("Intensity            = ",p.Intensity))
  println(io,string("Classification       = ",p.Classification))
end

function PulsesToRays(p::Vector{PulseRecord}, h::PulseWavesHeader)
    X0 = h.XOffset
    X1 = h.XScaleFactor
    Y0 = h.YOffset
    Y1 = h.YScaleFactor
    Z0 = h.ZOffset
    Z1 = h.ZScaleFactor
    rays = Array{Float64}(undef,(length(p),2,3))
    for i ∈ 1:length(p)
        rays[i,1,1] = X0 + X1 * p[i].AnchorX
        rays[i,1,2] = Y0 + Y1 * p[i].AnchorY
        rays[i,1,3] = Z0 + Z1 * p[i].AnchorZ
        rays[i,2,1] = X1 * (p[i].TargetX - p[i].AnchorX) / 1000 #dx
        rays[i,2,2] = Y1 * (p[i].TargetY - p[i].AnchorY) / 1000 #dy
        rays[i,2,3] = Z1 * (p[i].TargetZ - p[i].AnchorZ) / 1000 #dz
    end
    rays
end

function RescalePulseTimestamps(p::Vector{PulseRecord}, h::PulseWavesHeader)
    T0 = h.TOffset
    T1 = h.TScaleFactor
    Times = Vector{Float64}(undef,length(p))
    for i ∈ 1:length(p)
        Times[i] = T0 + T1 * p[i].GPSTimestampT
    end
    Times
end
