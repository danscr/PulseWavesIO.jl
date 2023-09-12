function load(f_pls::File{format"PLS"})
    open(f_pls) do s_pls
        f_wvs = replace(filename(f_pls), "$(file_extension(f_pls))" => ".wvs")
        @assert(isfile(f_wvs), "No Waves (.wvs) file was found in the given path. Waves file must be in same directory as Pulses file.")
        open(f_wvs) do s_wvs
	    #read header
    	    header = read(s_pls, PulseWavesHeader)
    	    n = header.NumberOfPulses
    	    #Number of pulses can be -1, in which case it needs to be inferred
    	    @assert(n > 0, "Unspecified number of pulses not yet implemented.")

    	    #Read pulses
    	    pulses = Vector{PulseRecord}(undef, n)
    	    @assert(position(s_pls) == header.OffsetToPulseData, "Unexpected (possibly user-defined) bytes between header and pulse records.")
    	    for i in 1:n
    	        pulses[i] = read(s_pls, PulseRecord)
    	    end

    	    #Sanity checks to ensure the .pls file is read correctly
    	    pos = position(s_pls)
    	    seekend(s_pls)
    	    @assert(header.NumberOfAppendedVariableLengthRecords == 0, "Appended variable length records not yet implemented.")
    	    # Even if n_AVLR==0, there has to be one mandatory AVLR with only a header that has 96 bytes
    	    @assert(position(s_pls) - pos == 96, "More than the single mandatory appended variable length record suspected.")
    	    #Start reading the mandatory AVLR header to make sure everything is correct
    	    seek(s_pls, pos)
	    MandatoryAppendedVariableLengthRecord = read(s_pls,96)
	    @assert(only(reinterpret(UInt64, MandatoryAppendedVariableLengthRecord[25:32])) == 0, "Mandatory appended variable length record has non-zero length after header! Something is wrong...") #record length after header must be 0
    	    @assert(eof(s_pls),"EOF of Pulse file not reached.")

    	    # Segments can only be fixed size for now, the 8 lines below are checking that only fixed segmenting is used
    	    pulseDescriptorIndices = unique(map(x -> x.PulseDescriptorIndex, pulses))
    	    pulseDescriptors = filter(x -> x.RecordID - 200000 in pulseDescriptorIndices, header.VariableLengthRecords)
    	    for i in 1:length(pulseDescriptors)
    	      local SamplingRecords = pulseDescriptors[i].Data.Sampling
    	      #for j in 1:length(SamplingRecords)
    	      #  @assert(SamplingRecords[j].BitsForNumberOfSegments == 0, "This file contains pulses with variable segmenting. This is not yet implemented.")
    	      #end
    	    end

    	    # ######
    	    # Read waves
    	    # ######
    	    wv_header = read(s_wvs, WavesHeader)
    	    waves = Vector{Vector{WaveRecord}}(undef, n)
	    #WavePulseDescriptorIndices[i] = Vector{Int8}(undef,n)
    	    #samplingType = Vector{Vector{UInt8}}(undef, n)
    	    for i ∈ 1:n
	      waves[i] = readWave(s_wvs, pulses[i], header, false)
    	    end
	header, pulses, wv_header, waves, MandatoryAppendedVariableLengthRecord
	end
    end
end


function save(f_pls::File{format"PLS"}, header::PulseWavesHeader, pulses::Vector{PulseRecord}, wv_header::WavesHeader, waves::Vector{Vector{WaveRecord}}, AVLR::Vector{UInt8})
  open(f_pls, "w") do s_pls
    write(s_pls, header)
    for i ∈ 1:length(pulses)#header.NumberOfPulses
      write(s_pls, pulses[i])
    end
    write(s_pls, AVLR)
  end
  f_wvs = replace(filename(f_pls), "$(file_extension(f_pls))" => ".wvs")
  open(f_wvs, "w") do s_wvs
    write(s_wvs, wv_header)
    # writing waves not yet implemented
    #for i ∈ 1:length(waves)
    #  for j ∈ 1:length(waves[i])
    #    write(s_wvs, waves[i][j])
    #  end
    #end
  end
  nothing
end

