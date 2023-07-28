function readstring(io, nb::Integer)
    bytes = read(io, nb)
    # strip possible null bytes
    lastchar = findlast(bytes .!= 0)
    if lastchar == nothing
        return ""
    else
        return String(bytes[1:lastchar])
    end
end

function writestring(io, str::AbstractString, nb::Integer)
    n = length(str)
    npad = nb - n
    if npad < 0
        error("string too long")
    elseif npad == 0
        write(io, str)
    else
        writestr = string(str * "\0"^npad)
        write(io, writestr)
    end
end

function dBToPowerRatio(x::AbstractFloat)
  10^(x/10)
end

function InterpolateLUT(lut::Vector{Float32}, idx::AbstractFloat)::Float64
  i0 = Int(floor(idx))
  i1 = Int(ceil(idx))
  if i0 == 0
    return lut[i1]
  elseif i0 != i1
    return lut[i0] + (lut[i1] - lut[i0]) * (idx - i0) / (i1 - i0)
  else
    return lut[i0]
  end
end

#ensure to rotate the dem 90 deg. right!
function getDEMvalue(x,y,dem::Dict)
  xi = Int(round(x - dem[:header]["XLLCORNER"]))
  yi = Int(round(y - dem[:header]["YLLCORNER"]))
  dem[:DEM][xi,yi]
end

function RayZenithAngle(ray::Matrix{Float64})
  hxy = sqrt(ray[2,1]^2 + ray[2,2]^2)
  α = atan(ray[2,3], hxy)
  α + π/2
end

