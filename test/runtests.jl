###########################################################################
#   Copyright (C) 2025 by Pierre Abbat                                    #
#   phma@bezitopo.org                                                     #
#   This file is part of Mumzel.                                          #
#                                                                         #
#   Mumzel is free software; you can redistribute it and/or modify        #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 3 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   Mumzel is distributed in the hope that it will be useful,             #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with Mumzel; if not, see <http://www.gnu.org/licenses/>.        #
###########################################################################

using MumzelCode,Test,Printf
using MumzelCode:c43434,c33435,c43425,c33525,perm60,perm20,permute,invZel

@test MumzelCode.invLetter[0x65]==0x3a

@test cycleType(Codeword([0,1,2,3,4,0]))==Codeword([1,1,1,1,1,0])
@test cycleType(Codeword([0,4,1,2,3,0]))==Codeword([1,4,4,4,4,0])
@test cycleType(Codeword([0,4,3,2,1,0]))==Codeword([1,2,6,6,2,0])

function permutedAllDifferent(c::Codeword,v::AbstractVector{<:Integer})
  s=Set{Codeword}()
  for i in v
    push!(s,permute(c,i))
  end
  length(s)==length(v)
end

@test permutedAllDifferent(c43425,perm60)
@test permutedAllDifferent(c33525,perm60[0:29])
@test permutedAllDifferent(c33435,perm20)
@test permutedAllDifferent(c43434,perm20[0:9])

function testEncode1(n::Unsigned)
  cw=encode(n,32)
  int=codewordInt(cw)
  if count_ones(int)!=18
    @printf "%08x->%09x has %d bits\n" n int count_ones(int)
  end
  dec,bits,codeBits,upsideDown=decode(cw,false)
  if dec!=n
    @printf "%08x->%09x misdecoded as %08x\n" n int dec
  end
  count_ones(int)==18 && dec==n && bits==32 && codeBits==36
end

function testEncode()
  for i in 0x0:0xffff
    if !testEncode1(i*0x69969669)
      return false
    end
  end
  true
end

@test testEncode()

function testDistinctMumPart(z::Integer)
  zelPart=invZel[z]
  i=0
  j=0
  mumDict=Dict{UInt64,UInt32}()
  for mumPart in 0x0:0x1ffff
    code=codewordInt(encode(mumPart<<14+zelPart,32))
    if haskey(mumDict,code)
      if i==j*j
	@printf "%x has the same code as %x\n" mumDict[code] mumPart
	j+=1
      end
      i+=1
    end
    mumDict[code]=mumPart
  end
  i==0
end

@test testDistinctMumPart(0o11111)
@test testDistinctMumPart(0o22222)
@test testDistinctMumPart(0o33333)
@test testDistinctMumPart(0o44444)
@test testDistinctMumPart(0o01234)
@test testDistinctMumPart(0o12345)
@test testDistinctMumPart(0o23456)
