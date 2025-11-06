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

module MumzelCode
using OffsetArrays,StaticArrays
export Codeword,permcode,permoct,cycleType

const letter=OffsetVector(
# 0101010 1010100 1010001 1000101 0010101 1100000 1000001 0000011
# 1100100 0110010 0011001 1010010 0110001 0110000 0100001 1000010
# 1001100 0100110 0010011 0100101 1000110 0011000 1000100 0010001
# 1101000 0110100 0011010 0001101 0101001 0001100 1001000 0100100
# 1011000 0101100 0010110 0001011 1001010 0000110 0010010 0001001
# 1000011 1100010 1001001 0100011 1100001 1010000 0010100 0000101
# 0000111 0001110 0011100 0111000 1110000 0101000 0100010 0001010
# Rows 0 to 6 consist of 5 numbers with 3 bits set and 3 numbers with 2 bits set,
# followed by their complements. The row number comes from the zel code, which
# is a base-7 code for 14 bits of the number (0-16383). The column comes from
# other parts of the number. Row 7 is not used for encoding; if it appears in
# decoding, the codeword is invalid or a framing error. A framing error of the
# idle codes appears as two adjacent row-0 codes opposite a row-7 code.
[
  0x2a,0x54,0x51,0x45,0x15, 0x60,0x41,0x03, 0x55,0x2b,0x2e,0x3a,0x6a, 0x1f,0x3e,0x7c,
  0x64,0x32,0x19,0x52,0x31, 0x30,0x21,0x42, 0x1b,0x4d,0x66,0x2d,0x4e, 0x4f,0x5e,0x3d,
  0x4c,0x26,0x13,0x25,0x46, 0x18,0x44,0x11, 0x33,0x59,0x6c,0x5a,0x39, 0x67,0x3b,0x6e,
  0x68,0x34,0x1a,0x0d,0x29, 0x0c,0x48,0x24, 0x17,0x4b,0x65,0x72,0x56, 0x73,0x37,0x5b,
  0x58,0x2c,0x16,0x0b,0x4a, 0x06,0x12,0x09, 0x27,0x53,0x69,0x74,0x35, 0x79,0x6d,0x76,
  0x43,0x62,0x49,0x23,0x61, 0x50,0x14,0x05, 0x3c,0x1d,0x36,0x5c,0x1e, 0x2f,0x6b,0x7a,
  0x07,0x0e,0x1c,0x38,0x70, 0x28,0x22,0x0a, 0x78,0x71,0x63,0x47,0x0f, 0x57,0x5d,0x75,
  0x01,0x02,0x04,0x08,0x10,0x20,0x40, 0x00, 0x7e,0x7d,0x7b,0x77,0x6f,0x5f,0x3f, 0x7f
],-1)

function invertLetter()
  inv=OffsetVector(fill(0xff,128),-1)
  for i in 0x00:0x7f
    inv[letter[i]]=i
  end
  inv
end

const invLetter=invertLetter()

# [1:5] are the letters, [6] is the sign bit
# The bytes can be the actual bit patters of letters (0x00-0xff), indices of
# letters (0x00-0xff), permutations (0x0-0x4), or bit counts (0x2-0x5).
# Permutations are written in reverse order (43210 is the identity).
Codeword=SVector{6,UInt8}

# Make the zel code table. Zel codes are 5-digit base-7 numbers.
# All data zel codes have at most two zeros and do not consist
# entirely of fives and sixes.
# Codes with four or five zeros are reserved for idle channel codes.
function makezel()
  zel=OffsetVector(fill(0xffff,16384),-1)
  invZel=OffsetVector(fill(0xffff,32768),-1)
  n=0
  l=OffsetVector([0,0,0,0,0],-1)
  for i in 0:7^5-1
    r=i
    m=m2=0
    for j in 0:4
      l[j]=r%7
      r÷=7
      m+=l[j]==0
      m2+=l[j]>=5
    end
    if m<3 && m2<5
      r=0
      for j in 4:-1:0
	r=8*r+l[j]
      end
      invZel[r]=n
      zel[n]=r
      n+=1
    end
  end
  @assert n==16384
  zel,invZel
end

const zel,invZel=makezel()

function permute(cword::Codeword,perm::Integer)
  mcword=MVector(cword)
  for i in 0:9
    if (perm>>i)&1==1
      mcword[i%5+1],mcword[(i+2)%5+1]=mcword[(i+2)%5+1],mcword[i%5+1]
    end
  end
  SVector(mcword)
end

"""
    permcode(cword::Codeword)

Given the number of 1-bits in each letter of a codeword, returns a 10-bit number
which can be looked up in a table to find how to undo the permutation and unpack
the bits. If a letter has <2 1-bits, returns a negative number. If a letter has
>5 ones, returns garbage.
"""
function permcode(cword::Codeword)
  pc=0
  for i in 1:5
    pc|=(cword[i]-2)<<(2*i-2)
  end
  pc
end

"""
    permoct(cword::Codeword)

Given a permutation of 0-4, returns a 15-bit number which can be printed in octal.
"""
function permoct(cword::Codeword)
  pc=0
  for i in 1:5
    pc|=cword[i]<<(3*i-3)
  end
  pc
end

"""
    cycleType(cword::Codeword)

Given a permutation of 43210 (`cword[6]` is ignored), returns a `Codeword` in
which each byte in `[1:5]` is the cycle length of the byte in `cword`, except
if there are two 2-cycles, in which case one is distinguished as 6.
"""
function cycleType(cword::Codeword)
  cycle=[0,0,0,0,0,0]
  for i in 1:5
    if cycle[i]==0
      j=i
      k=0
      while (j!=i || k==0) && k<128
	j=cword[j]+1
	k+=1
      end
      if k==2
	for l in 1:5
	  if cycle[l]==k
	    k+=4
	  end
	end
      end
      j=i
      for l in 1:k
	cycle[j]=k
	j=cword[j]+1
      end
    end
  end
  Codeword(cycle)
end

end # module MumzelCode
