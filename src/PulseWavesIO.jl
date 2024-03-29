module PulseWavesIO
using FileIO

#for gaussian decomp
import LinearAlgebra
using Distributions
using Statistics
import Peaks
import LsqFit
import ImageFiltering

export
    # Types
    PulseWavesHeader,
    PulseRecord,
    PulseWavesVariableLengthRecord,
    WavesHeader,
    WavesRecord,
    PulsesToRays,
    RescalePulseTimestamps,
    WaveToCoverFractionDistribution

    # Functions
    include("meta.jl")
    include("Util.jl")
    include("vlr.jl")
    include("header.jl")
    include("pulses.jl")
    include("waves.jl")
    include("fileio.jl")
    include("WaveformAnalysis.jl")

function __init__()
    # these should eventually go in
    # https://github.com/JuliaIO/FileIO.jl/blob/master/src/registry.jl
    add_format(format"PLS", "PulseWavesPulse", ".pls", [:PulseWavesIO])
end

end
