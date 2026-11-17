module TransLines
# Simulates a transmission line
using OffsetArrays,Random
export TransLine,filter!,transmit!,Receiver,receiveBits!
export Noisemaker,noise!

fastRate::Float64=0.2
slowRate::Float64=1/8/1331

rng=RandomDevice()

mutable struct Noisemaker
  next	::Int
  int	::Int
  Noisemaker()=new(0,0)
end

const fibo::OffsetVector{Int16}=OffsetVector(
  [13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711]
  ,-1) # start with 13 because 8 is a power of 2

function delay(n::Integer)
  2^(n÷16+4)+fibo[n%16]
end

function noise!(nm::Noisemaker)
  nm.next-=1
  if nm.next<0
    nm.next+=delay(rand(rng,UInt8))
    nm.int=rand(rng,UInt8)
  end
  if nm.next<2
    real((7im/8-1/2)^(nm.int+nm.next))*64
  else
    0.0
  end
end

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
  ret=Int16[]
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

recvRate::Float64=1/8/11

mutable struct Receiver
  avg	::Float64
  buf	::Vector{Float64}
  phase	::Int
  pherr	::Int
  Receiver()=new(0,Float64[],16,0)
end

function receiveBits!(rcv::Receiver,stream::Vector{<:Integer})
  ret=Int8[] # pack 8 to a byte later
  for sample in stream
    sample+=0.5
    push!(rcv.buf,sample-rcv.avg)
    rcv.avg=sample*recvRate+rcv.avg*(1-recvRate)
    if length(rcv.buf)>16
      popfirst!(rcv.buf)
    end
    rcv.phase-=1
    if rcv.phase==0
      push!(ret,rcv.buf[3]+rcv.buf[5]+rcv.buf[7]>0)
      # TODO should it convert the three samples to bits and vote?
      if rcv.buf[8]*rcv.buf[10]<0
	rcv.pherr+=(rcv.buf[8]*rcv.buf[9]<0) ? 1 : -1
      end
      if rcv.buf[7]*rcv.buf[3]<0
	rcv.pherr+=(rcv.buf[7]*rcv.buf[5]<0) ? 1 : -1
      end
      if rcv.pherr>10
	rcv.phase=7
	rcv.pherr=0
      elseif rcv.pherr<-10
	rcv.phase=9
	rcv.pherr=0
      else
	rcv.phase=8
      end
    end
  end
  ret
end

end
