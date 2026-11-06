module TransLines
# Simulates a transmission line
export TransLine,filter!,transmit!

fastRate::Float64=0.5
slowRate::Float64=1/8/11

mutable struct TransLine
  fast	::Float64
  slow	::Float64
  TransLine()=new(0,0)
end

function filter!(tl::TransLine,sample::AbstractFloat)
  tl.fast=tl.fast*(1-fastRate)+sample*fastRate
  tl.slow=tl.slow*(1-slowRate)+sample*slowRate
  tl.fast-tl.slow
end

"""
    transmit!(tl::TransLine,stream::Vector{UInt8}[,nbits::Integer])

`stream` is little-endian and is a stream of codewords consisting of bits.
Returns a stream of `Int8` samples, eight samples per bit. If `nbits` is not
specified, it is the number of bits in `stream`.
"""
function transmit!(tl::TransLine,stream::Vector{UInt8},nbits::Integer)
  ret=Int8[]
  for i in 8:nbits+7
    b=Float64((stream[i÷8]>>(i%8))&1)
    for j in 1:8
      sample=filter!(tl,b)
      push!(ret,floor(sample*128))
    end
  end
  ret
end

function transmit!(tl::TransLine,stream::Vector{UInt8})
  transmit!(tl,stream,length(stream)*8)
end

end
