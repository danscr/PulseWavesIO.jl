using PulseWavesIO, FileIO
using Test

# top level tests
@testset "PulseWavesIO.jl" begin
  header, pulses, wv_header, waves, samplingtype, avlr = load("test.pls")
  @test typeof(header) == LidarLeafArea.PulseWavesHeader
  r = PulsesToRays(pulses,header)
  t = RescalePulseTimestamps(pulses, header)
  @test typeof(pulses) == Vector{PulseRecord}
  @test length(pulses) == 15
  @test pulses[1] == PulseRecord(129863735377, 60, 23500619, 80005128, 126118, 23504847, 80003736, 111808, 8213, 8251, 0x00, 0x00)
end
