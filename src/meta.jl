""""
From LasIO.jl.
Generate read (unpack) method for structs.
Approximately evaluates to:
function Base.read(io::IO, t::Type{LasPoint0})
    (LasPoint0)(read(io, Int32), read(io, Int32),...)
"""
function generate_read(T::Type)
    fc = fieldcount(T)
    types = [fieldtype(T, i) for i = 1:fc]

    # Create unpack function expression
    function_expression = :(function Base.read(io::IO, t::Type{$T}) end)

    # Create Type call expression and add parameters
    type_expression = :(($T)())
    for t in types
        read_expression = :(read(io, $t))
        append!(type_expression.args, 0)  # dummy with known length
        type_expression.args[end] = read_expression
    end

    # Replace empty function body with Type call
    function_expression.args[2] = type_expression

    eval(function_expression)
end

"""
From LasIO.jl.
Generate write (pack) method for structs.
Approximately evaluates to:
function Base.write(io::IO, T::LasPoint0)
    write(io, T.x)
    write(io, T.y)
    write(io, T.z)
    ...)
"""
function generate_write(T::Type)
    # Create pack function expression
    function_expression = :(function Base.write(io::IO, T::$T) end)

    body_expression = quote end
    for t in fieldnames(T)
        append!(body_expression.args, 0)  # dummy with known length
        write_expression = :(write(io, T.$t))
        body_expression.args[end] = write_expression
    end

    # Return nothing at the end
    append!(body_expression.args, 0)  # dummy with known length
    body_expression.args[end] = :(nothing)

    # Replace empty function body with write calls
    function_expression.args[2] = body_expression

    eval(function_expression)
end

function generate_io(T::Type)
    generate_read(T)
    generate_write(T)
end

"""
From LasIO.jl.
Generate IO expressions macro.
Approximately evaluates to:
generate_io(LasPoint0)
call order
gen_io -> generate_io -> generate_read/write
"""
macro gen_io(typ::Expr)
    T = typ.args[2]
    if Meta.isexpr(T, :(<:))
        T = T.args[1]
    end
    if Meta.isexpr(T, :curly)
        T = T.args[1]
    end

    ret = Expr(:toplevel, :(Base.@__doc__ $(typ)))
    push!(ret.args, :(generate_io($T)))
    return esc(ret)
end


