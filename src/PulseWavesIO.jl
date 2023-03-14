module PulseWavesIO
using FileIO


export
    # Types
    PulseWavesHeader,
    PulseRecord,
    PulseWavesVariableLengthRecord,
    WavesHeader,
    WavesRecord

    # Functions
    include("meta.jl")
    include("Util.jl")
    include("vlr.jl")
    include("header.jl")
    include("pulses.jl")
    include("waves.jl")
    include("fileio.jl")

function __init__()
    # these should eventually go in
    # https://github.com/JuliaIO/FileIO.jl/blob/master/src/registry.jl
    add_format(format"PLS", "PulseWavesPulse", ".pls", [:LidarLeafArea])
end

end
