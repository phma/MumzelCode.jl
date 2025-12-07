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

function testEncode1(n::Unsigned,b::Integer)
  cw=encode(n,b)
  int=codewordInt(cw)
  if count_ones(int)!=18
    @printf "%08x->%09x has %d bits\n" n int count_ones(int)
  end
  dec,bits,codeBits,upsideDown=decode(cw,false)
  if dec!=n
    @printf "%08x->%09x misdecoded as %08x\n" n int dec
  end
  count_ones(int)==18 && dec==n && bits==b && codeBits==36
end

function testEncode()
  for i in 0x0:0xffff
    if !testEncode1(i*0x69969669,32)
      return false
    end
  end
  for i in 0x0:0xffff
    if !testEncode1((i*0x131ec09)%0x1eefe01,25)
      return false
    end
  end
  for i in 0x0:0xffff
    if !testEncode1(i,16)
      return false
    end
  end
  for i in 0x0:0xff
    if !testEncode1(i,8)
      return false
    end
  end
  for i in 0x0:0x7f
    if !testEncode1(i,7)
      return false
    end
  end
  for i in 0x0:0x3f
    if !testEncode1(i,6)
      return false
    end
  end
  for i in 0x0:0x1f
    if !testEncode1(i,5)
      return false
    end
  end
  for i in 0x0:0xf
    if !testEncode1(i,4)
      return false
    end
  end
  for i in 0x0:0x7
    if !testEncode1(i,3)
      return false
    end
  end
  for i in 0x0:0x3
    if !testEncode1(i,2)
      return false
    end
  end
  for i in 0x0:0x1
    if !testEncode1(i,1)
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
  for mumPart in 0x0:0x205ff
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

function testSyncWord()
  stream=UInt8[]
  for i in 1:5 # It takes 8 syncwords to synchronize and 1 more to flip.
    appendCodeword!(stream,encode(0x1,0))
    appendCodeword!(stream,encode(0x0,0))
  end
  appendCodeword!(stream,encode(0x636e7973,32))
  for i in 0:72
    dec=decodeStream(stream,i,true)
    if dec[end][1]!=0x636e7973
      return false
    end
    dec=decodeStream(stream,i,false)
    if dec[end][1]!=0x636e7973
      return false
    end
  end
  true
end

@test testSyncWord()

function testInvalidCodeword(cw::Codeword)
  dec,bits,codeBits,upsideDown=decode(cw,false)
  bits<0
end

function testInvalidCodeword(cw::Vector{<:Integer})
  testInvalidCodeword(Codeword(cw))
end

@test testInvalidCodeword([0x1f,0x54,0x7f,0x7f,0x55,0])
# This has a run of 16 bits. Real Mumzel has no run longer than 11.
@test testInvalidCodeword([0x46,0x36,0x17,0x4e,0x2c,0])
# partB is 249; the maximum valid is 223, except for syncwords.
@test testInvalidCodeword([0x58,0x4a,0x16,0x35,0x69,1])
# This is row 120 of the 125 rows of pattern 43434, but only 112 rows are used.
@test testInvalidCodeword([0x16,0x69,0x4a,0x6d,0x12,1])
# This is row 265 of the 270 rows of pattern 43425, but only 256 rows are used.
