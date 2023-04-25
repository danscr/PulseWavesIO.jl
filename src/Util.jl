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
