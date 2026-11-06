module TransLines
# Simulates a transmission line
export TransLine,filter!

fastRate::Float64=0.5
slowRate::Float64=1/8/11

mutable struct TransLine
  fast	::Float64
  slow	::Float64
  TransLine()=new(-2,-2)
end

function filter!(tl::TransLine,sample::AbstractFloat)
  tl.fast=tl.fast*(1-fastRate)+sample*fastRate
  tl.slow=tl.slow*(1-slowRate)+sample*slowRate
  tl.fast-tl.slow
end

end
