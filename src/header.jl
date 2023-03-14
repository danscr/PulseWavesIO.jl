"""
PulseWaves file header according to the specifications in
https://github.com/PulseWaves/Specification/blob/master/specification.rst
Several attributes are left out for now, these are TODO
"""
struct PulseWavesHeader
    fileSignature       ::AbstractString
    globalParameters    ::UInt32
    FileSourceID        ::UInt32
    ProjectGUID_1       ::UInt32
    ProjectGUID_2       ::UInt16
    ProjectGUID_3       ::UInt16
    ProjectGUID_4       ::AbstractString
    SystemIdentifier    ::AbstractString
    GeneratingSoftware  ::AbstractString
    FileCreationDayOfYear   ::UInt16
    FileCreationYear        ::UInt16
    VersionMajor        ::UInt8
    VersionMinor        ::UInt8
    HeaderSize	    	::UInt16
    OffsetToPulseData  	::Int64
    NumberOfPulses  	::Int64
    PulseFormat  	::UInt32
    PulseAttributes  	::UInt32
    PulseSize  	    	::UInt32
    PulseCompression    ::UInt32
    Reserved            ::Int64
    NumberOfVariableLengthRecords	::UInt32
    NumberOfAppendedVariableLengthRecords ::Int32
    TScaleFactor  	::Float64
    TOffset  		::Float64
    MinT  	    	::Int64
    MaxT  	    	::Int64
    XScaleFactor  	::Float64
    YScaleFactor  	::Float64
    ZScaleFactor  	::Float64
    XOffset  	 	::Float64
    YOffset 	 	::Float64
    ZOffset 	 	::Float64
    MinX 	 	    ::Float64
    MaxX 	 	    ::Float64
    MinY 	 	    ::Float64
    MaxY 	 	    ::Float64
    MinZ 	 	    ::Float64
    MaxZ 	 	    ::Float64
    VariableLengthRecords ::Vector{PulseWavesVariableLengthRecord}
    UserDefinedBytes::Vector{UInt8}
end

function Base.show(io::IO, h::PulseWavesHeader)
 println(io, string("\tfileSignature	= ",h.fileSignature 		))
 println(io, string("\tglobalParameters	= ",h.globalParameters 		))
 println(io, string("\tFileSourceID	= ",h.FileSourceID 		))
 println(io, string("\tProjectGUID_1	= ",h.ProjectGUID_1 		))
 println(io, string("\tProjectGUID_2	= ",h.ProjectGUID_2 		))
 println(io, string("\tProjectGUID_3	= ",h.ProjectGUID_3 		))
 println(io, string("\tProjectGUID_4	= ",h.ProjectGUID_4 		))
 println(io, string("\tSystemIdentifier	= ",h.SystemIdentifier 		))
 println(io, string("\tGeneratingSoftware	= ",h.GeneratingSoftware 		))
 println(io, string("\tFileCreationDayOfYear	= ",h.FileCreationDayOfYear		))
 println(io, string("\tFileCreationYear	= ",h.FileCreationYear 		))
 println(io, string("\tVersionMajor	= ",h.VersionMajor 		))
 println(io, string("\tVersionMinor	= ",h.VersionMinor 		))
 println(io, string("\tHeaderSize	= ",h.HeaderSize	 		))
 println(io, string("\tOffsetToPulseData	= ",h.OffsetToPulseData 		))
 println(io, string("\tNumberOfPulses	= ",h.NumberOfPulses 		))
 println(io, string("\tPulseFormat	= ",h.PulseFormat 		))
 println(io, string("\tPulseAttributes	= ",h.PulseAttributes 		))
 println(io, string("\tPulseSize	= ",h.PulseSize 	 		))
 println(io, string("\tPulseCompression	= ",h.PulseCompression 		))
 println(io, string("\tReserved	= ",h.Reserved 		))
 println(io, string("\tNumberOfVariableLengthRecords	= ",h.NumberOfVariableLengthRecords	))
 println(io, string("\tNumberOfAppendedVariableLengthRecords	= ",h.NumberOfAppendedVariableLengthRecords))
 println(io, string("\tTScaleFactor	= ",h.TScaleFactor 		))
 println(io, string("\tTOffset	= ",h.TOffset 	 		))
 println(io, string("\tMinT	= ",h.MinT 	 		))
 println(io, string("\tMaxT	= ",h.MaxT 	 		))
 println(io, string("\tXScaleFactor	= ",h.XScaleFactor 		))
 println(io, string("\tYScaleFactor	= ",h.YScaleFactor 		))
 println(io, string("\tZScaleFactor	= ",h.ZScaleFactor 		))
 println(io, string("\tXOffset	= ",h.XOffset 	 		))
 println(io, string("\tYOffset	= ",h.YOffset 	 		))
 println(io, string("\tZOffset	= ",h.ZOffset 	 		))
 println(io, string("\tMinX	= ",h.MinX 	 		))
 println(io, string("\tMaxX	= ",h.MaxX 	 		))
 println(io, string("\tMinY	= ",h.MinY 	 		))
 println(io, string("\tMaxY	= ",h.MaxY 	 		))
 println(io, string("\tMinZ	= ",h.MinZ 	 		))
 println(io, string("\tMaxZ	= ",h.MaxZ 	 		))
 #println(io, string("\tVariableLengthRecords	= ",h.VariableLengthRecords		))
 println(io, string("\tUserDefinedBytes	= ",h.UserDefinedBytes 		))
  println(io, string("WARNING: Variable Length Records are not shown, please print them separately!"))
end




function Base.read(io::IO, ::Type{PulseWavesHeader})
    fileSignature =             readstring(io, 16)
    @assert(fileSignature == "PulseWavesPulse", "File is not a valid PulseWaves (.pls) file! File signature does not match specifications.")
    globalParameters =          read(io,UInt32)
    FileSourceID =	            read(io,UInt32)
    ProjectGUID_1 =             read(io, UInt32)
    ProjectGUID_2 =             read(io, UInt16)
    ProjectGUID_3 =             read(io, UInt16)
    ProjectGUID_4 =             readstring(io,8)
    SystemIdentifier = 	        readstring(io,64)
    GeneratingSoftware = 	    readstring(io,64)
    FileCreationDayOfYear =     read(io,UInt16)
    FileCreationYear = 	        read(io,UInt16)
    VersionMajor =	            read(io,UInt8)
    VersionMinor = 	            read(io,UInt8)
    HeaderSize =                read(io,Int16)
    OffsetToPulseData =     	read(io,Int64)
    NumberOfPulses =        	read(io,Int64)
    PulseFormat = 	            read(io,UInt32)
    PulseAttributes = 	        read(io,UInt32)
    PulseSize = 		        read(io,UInt32)
    PulseCompression = 	        read(io,UInt32)
    Reserved = 		            read(io,Int64)
    NumberOfVariableLengthRecords = read(io,Int32)
    NumberOfAppendedVariableLengthRecords = read(io,Int32)
    TScaleFactor = 	            read(io,Float64)
    TOffset = 		            read(io,Float64)
    MinT = 		                read(io,Int64)
    MaxT = 		                read(io,Int64)
    XScaleFactor = 	            read(io,Float64)
    YScaleFactor = 	            read(io,Float64)
    ZScaleFactor = 	            read(io,Float64)
    XOffset = 	 	            read(io,Float64)
    YOffset =	 	            read(io,Float64)
    ZOffset =	 	            read(io,Float64)
    MinX =	 	                read(io,Float64)
    MaxX =	 	                read(io,Float64)
    MinY =	 	                read(io,Float64)
    MaxY =	 	                read(io,Float64)
    MinZ =	 	                read(io,Float64)
    MaxZ =	 	                read(io,Float64)

    VariableLengthRecords = [read(io, PulseWavesVariableLengthRecord) for i ∈ 1:NumberOfVariableLengthRecords]

    UserDefinedBytes = read(io, OffsetToPulseData - position(io))

    PulseWavesHeader(
                     fileSignature,
                     globalParameters,
                     FileSourceID,
                     ProjectGUID_1,
                     ProjectGUID_2,
                     ProjectGUID_3,
                     ProjectGUID_4,
                     SystemIdentifier,
                     GeneratingSoftware,
                     FileCreationDayOfYear,
                     FileCreationYear,
                     VersionMajor,
                     VersionMinor,
		             HeaderSize,
  		             OffsetToPulseData,
  		             NumberOfPulses,
  		             PulseFormat,
  		             PulseAttributes,
  		             PulseSize,
                     PulseCompression,
                     Reserved,
  		             NumberOfVariableLengthRecords,
  		             NumberOfAppendedVariableLengthRecords,
  		             TScaleFactor,
  		             TOffset,
  		             MinT,
  		             MaxT,
  		             XScaleFactor,
  		             YScaleFactor,
  		             ZScaleFactor,
  		             XOffset,
  		             YOffset,
  		             ZOffset,
  		             MinX,
  		             MaxX,
  		             MinY,
  		             MaxY,
  		             MinZ,
  		             MaxZ,
                     VariableLengthRecords,
                     UserDefinedBytes
  		    )
end


function Base.write(io::IO, h::PulseWavesHeader)
    writestring(io, h.fileSignature, 16)
    write(io, h.globalParameters)
    write(io, h.FileSourceID)
    write(io, h.ProjectGUID_1)
    write(io, h.ProjectGUID_2)
    write(io, h.ProjectGUID_3)
    writestring(io, h.ProjectGUID_4, 8)
    writestring(io, h.SystemIdentifier, 64)
    writestring(io, h.GeneratingSoftware, 64)
    write(io, h.FileCreationDayOfYear)
    write(io, h.FileCreationYear)
    write(io, h.VersionMajor)
    write(io, h.VersionMinor)
    write(io, h.HeaderSize)
    write(io, h.OffsetToPulseData)
    write(io, h.NumberOfPulses)
    write(io, h.PulseFormat)
    write(io, h.PulseAttributes)
    write(io, h.PulseSize)
    write(io, h.PulseCompression)
    write(io, h.Reserved)
    @assert length(h.VariableLengthRecords) == h.NumberOfVariableLengthRecords
    write(io, h.NumberOfVariableLengthRecords)
    write(io, h.NumberOfAppendedVariableLengthRecords)
    write(io, h.TScaleFactor)
    write(io, h.TOffset)
    write(io, h.MinT)
    write(io, h.MaxT)
    write(io, h.XScaleFactor)
    write(io, h.YScaleFactor)
    write(io, h.ZScaleFactor)
    write(io, h.XOffset)
    write(io, h.YOffset)
    write(io, h.ZOffset)
    write(io, h.MinX)
    write(io, h.MaxX)
    write(io, h.MinY)
    write(io, h.MaxY)
    write(io, h.MinZ)
    write(io, h.MaxZ)
    for i ∈ 1:h.NumberOfVariableLengthRecords
        write(io, h.VariableLengthRecords[i])
    end
    write(io, h.UserDefinedBytes)
    nothing
end
