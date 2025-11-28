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
using OffsetArrays,StaticArrays,Printf
export Codeword,permcode,permoct,cycleType,makeperms,makeperms2,permstr,outComb
export halfEncode,encode,codewordInt,perm20Inverses,perm60Inverses
export halfDecode,decode

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
# The bytes can be the actual bit patterns of letters (0x00-0xff), indices of
# letters (0x00-0xff), permutations (0x0-0x4), or bit counts (0x2-0x5).
# Permutations are written in reverse order (43210 is the identity).
Codeword=SVector{6,UInt8}

const c43434=Codeword([4,3,4,3,4,0])
const c34343=Codeword([3,4,3,4,3,1])
const c33435=Codeword([5,3,4,3,3,0])
const c44342=Codeword([2,4,3,4,4,1])
const c43425=Codeword([5,2,4,3,4,0])
const c34352=Codeword([2,5,3,4,3,1])
const c33525=Codeword([5,2,5,3,3,0])
const c44252=Codeword([2,5,2,4,4,1])
const id=Codeword([0,1,2,3,4,5])

function codewordInt(c::Codeword)
  acc=UInt64(0)
  for i in 1:6
    acc|=UInt64(c[i])<<(7*(i-1))
  end
  acc
end

# Make the zel code table. Zel codes are 5-digit base-7 numbers.
# All data zel codes have at most two zeros and do not consist
# entirely of fives and sixes.
# Codes with four or five zeros are reserved for idle channel codes.
# Codes with a seven (which can't appear in a base-7 number) opposite two zeros,
# with the other two digits being anything but 7, are framing errors of the
# idle codes and are in the inverse zel table as 0x60nn, where nn is 25 (37)
# if the receiver should skip a bit and 23 (35) if it should repeat a bit.
# For *7*00, 0*7*0, and 00*7*, it sets the counter so that the next frame has
# the 7, which is at least 6 ones or 6 zeros in a row, centered on the parity
# bit.
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
  for i in 0:48
    word=0o00700+(i%7)<<9+(i÷7)<<3
    invZel[word]=0x6012
    word=(word<<3)&0o77770|(word>>12)&7
    invZel[word]=0x6019
    word=(word<<3)&0o77770|(word>>12)&7
    invZel[word]=0x6023
    word=(word<<3)&0o77770|(word>>12)&7
    invZel[word]=0x6025
    word=(word<<3)&0o77770|(word>>12)&7
    invZel[word]=0x602f
  end
  invZel[0]=0x7024 # bit 12 set means in sync
  zel,invZel
end

const zel,invZel=makezel()

# Permutations
#
# Permutations are written big-endian, e.g. 32104 is [4,0,1,2,3,0] (the sign bit
# in a Codeword is ignored).

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
    pc|=UInt16(cword[i])<<(3*i-3)
  end
  pc
end

function permstr(perm::Integer)
  str=""
  sv=permute(Codeword([0,1,2,3,4,5]),perm)
  for i in 5:-1:1
    str*='0'+sv[i]
  end
  str
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

function makeperms()
  ptable=OffsetVector(fill(0x0000,1024),-1)
  invPerm=OffsetVector(fill(0x0000,1024),-1)
  for i in 0x000:0x3ff
    a=permute(id,i)
    ptable[i]=permoct(a)
  end
  for j in 0x000:0x3ff
    for i in 0x000:0x3ff
      if ptable[i]==ptable[j] && count_ones(i)>=count_ones(j) && i!=j
	ptable[i]=j|32768
      end
    end
  end
  for i in 0x000:0x3ff
    if ptable[i]<32768
      a=permute(id,i)
      for j in 0x000:0x3ff
	if permute(a,j)==id
	  invPerm[j]=i
	end
      end
    end
  end
  (ptable,invPerm)
end

"""
    ptable::OffsetArray(::Vector{UInt16}, 0:1023)

Given the ten-bit code for a permutation, look up the 15-bit permoct code if it's
canonical, and another ten-bit code for the same permutation, with the high bit
set, if it isn't. Repeated lookup of codes with the high bit set will eventually
find the canonical code.
"""
const (ptable,invPerm)=makeperms()

"""
    invPerm::OffsetArray(::Vector{UInt16}, 0:1023)

Given the ten-bit code for a permutation, look up the ten-bit code for its
inverse. Looking up twice gets the canonical code for the original permutation.
"""
invPerm

"""
    octinx(oct::Integer)

Given the 15-bit octal code for a permutation, find its 10-bit code.
"""
function octinx(oct::Integer)
  for i in 0:1023
    if ptable[i]==oct
      return i
    end
  end
  return -1
end

function makeperms2()
  rot=0xdb # index of the permutation 04321, a rotation
  perm=OffsetMatrix(fill(0xfff,24,5),-1,-1)
  # There are 10 single swaps and 15 double swaps, 26 total involutions
  # including the identity. 5 of these are reflections. Of these, one is put
  # at perm[1,0] (0x039, 12340), and the others in perm[1,1:4]. That leaves
  # 22 involutions. The identity is put at perm[0,0], the multiplications by
  # 2 and 3 at perm[3,0] and perm[2,0] respectively, and the rest of the
  # involutions in perm[4:23,0]. All other permutations are computed by
  # rotating these 24.
  perm[:,0]=
  [ 0x000, 0x039, 0x0d0, 0x016
  , 0x04a, 0x001, 0x2d0, 0x006
  , 0x094, 0x002, 0x02d, 0x00c
  , 0x029, 0x004, 0x05a, 0x018
  , 0x052, 0x008, 0x0b4, 0x011
  , 0x025, 0x010, 0x168, 0x003
  ]
  for i in 0:23
    for j in 0:3
      perm[i,j+1]=octinx(permoct(permute(permute(id,perm[i,j]),rot)))
    end
  end
  perm
end

const perm=makeperms2()

"""
    outComb()

Output a table of the effects of permuting the four kinds of codewords.
"""
function outComb()
  for i in 0:119
    a=permoct(permute(id,perm[i%24,i÷24]))
    b=permoct(permute(c43434,perm[i%24,i÷24]))
    c=permoct(permute(c33435,perm[i%24,i÷24]))
    d=permoct(permute(c43425,perm[i%24,i÷24]))
    e=permoct(permute(c33525,perm[i%24,i÷24]))
    @printf "%3d %03x %05o %05o %05o %05o %05o\n" i perm[i%24,i÷24] a b c d e
  end
end

# Bits per letter, permutations thereof, and number of unpermuted words
# (not counting zel codes):
# 43434 10 5^5=3125	31250 125*250 100-16f 00-df 25088
# 33435 20 5^4*3=1875   37500 150*250 170-202 00-df 32928
# 43425 60 5^3*3^2=1125 67500 270*250 000-0ff 00-df 57344
# 33525 30 5^2*3^3=675  20250  81*250 000-202 e0-ff 16480
# Total		       156500 626*250
# Needed to encode 4 B 131072 (2^17; 14 bits are encoded by zel, 1 bit is encoded by flipping)
#
# Alternating group (60 members, starred are used for 33525):
# *43210 *04321 *10432 *21043 *32104
# *34120 *03412 *20341 *12034 *41203
# *02314 *40231 *14023 *31402 *23140
# *34201 *13420 *01342 *20134 *42013
# *03124 *40312 *24031 *12403 *31240
# *42301 *14230 *01423 *30142 *23014
#  10324  41032  24103  32410  03241
#  42130  04213  30421  13042  21304
#  10243  31024  43102  24310  02431
#  41320  04132  20413  32041  13204
#  02143  30214  43021  14302  21430
#  01234  40123  34012  23401  12340
# 2^nx+b group (20 members, starred are used for 43434):
# *43210 *04321 *10432 *21043 *32104
# *31420 *03142 *20314 *42031 *14203
#  12340  01234  40123  34012  23401
#  24130  02413  30241  13024  41302

const perm60=OffsetVector(
  [ 0x000, 0x0db, 0x055, 0x249, 0x03f
  , 0x2d0, 0x005, 0x09a, 0x074, 0x048
  , 0x02d, 0x00a, 0x059, 0x017, 0x090
  , 0x05a, 0x014, 0x093, 0x02e, 0x009
  , 0x0b4, 0x028, 0x126, 0x01d, 0x012
  , 0x168, 0x050, 0x04d, 0x03a, 0x024
  , 0x03c, 0x003, 0x0d8, 0x056, 0x24a
  , 0x02b, 0x125, 0x01e, 0x011, 0x06c
  , 0x018, 0x036, 0x06a, 0x0d2, 0x00f
  , 0x069, 0x0d1, 0x00c, 0x01b, 0x035
  , 0x099, 0x0ac, 0x04b, 0x078, 0x006
  , 0x027, 0x09c, 0x053, 0x04e, 0x039
  ],-1)

const perm20=OffsetVector(
  [ 0x000, 0x0db, 0x055, 0x249, 0x03f
  , 0x016, 0x091, 0x02c, 0x00b, 0x058
  , 0x039, 0x027, 0x09c, 0x053, 0x04e
  , 0x0d0, 0x00d, 0x01a, 0x034, 0x068
  ],-1)

const invPerm60=map(x->invPerm[x],perm60)
const invPerm20=map(x->invPerm[x],perm20)

function makeLetterPerm()
  letterPerm=OffsetVector(fill(0xff,1024),-1)
  for (i,c) in pairs(perm20[0:9])
    p=permute(c43434,c)
    letterPerm[permcode(p)]=i-1 # slicing an OffsetVector removes the offset
    letterPerm[1023-permcode(p)]=i-1
  end
  for (i,c) in pairs(perm20)
    p=permute(c33435,c)
    letterPerm[permcode(p)]=i+64
    letterPerm[1023-permcode(p)]=i+64
  end
  for (i,c) in pairs(perm60)
    p=permute(c43425,c)
    letterPerm[permcode(p)]=i+128
    letterPerm[1023-permcode(p)]=i+128
  end
  for (i,c) in pairs(perm60[0:29])
    p=permute(c33525,c)
    letterPerm[permcode(p)]=i+191
    letterPerm[1023-permcode(p)]=i+191
  end
  letterPerm
end

function perm60Inverses()
  for i in 0:59
    for j in 0:59
      if permute(permute(id,perm60[i]),perm60[j])==id
	print('*')
      else
	print(' ')
      end
    end
    println()
  end
end

function perm20Inverses()
  for i in 0:19
    for j in 0:19
      if permute(permute(id,perm20[i]),perm20[j])==id
	print('*')
      else
	print(' ')
      end
    end
    println()
  end
end

"""
    letterPerm::OffsetVector{UInt8}

Given the number returned by `permcode` (0-1023), returns an index to `perm20`
or `perm60`. The maximum valid value is 0xdd, so if the top three bits are all 1,
the permcode is invalid.
"""
const letterPerm=makeLetterPerm()

# Encoding

function halfEncode(sign::Integer,partA::Integer,partB::Integer,zelPart::Integer)
  # partA ranges from 0 to 625 and selects the kind of letters and permutations
  # and the 5 or 3 of the last three letters before permutation.
  # partB ranges from 0 to 249 and selects the 5 of the first two letters
  # and the 10 of the permutation.
  # zelPart ranges from 0 to 16384 and selects the 7 of all five letters.
  # If zelPart=16384, that means the 7s are all 0, i.e. an idle code.
  # sign inverts all the bits.
  #
  # The idle codes are half-encoded as follows:
  # Idle code 0 (which ends with 0) is 111111000101010101010101010101100000.
  # Sign bit=1
  # 0000011 1010101 0101010 1010101 0011111
  # 2-2-0 4-0-0 3-0-0 4-0-0 5-0-0
  # 00000 is not a valid zel code (all 0s is used only for idle).
  # 24345 is the 59th permutation (of 0-59, 12340) of 43425. Undoing the
  # permutation gives 4-0-0 3-0-0 4-0-0 2-2-0 5-0-0.
  # So partA=506 (last of the 6 groups of 45, plus 6 to get 020), partB=225
  # (the 9 of 59, times 25), and zelPart=16384.
  # Idle code 1 is 000000110101010101010101010110011111.
  # Sign bit=0
  # 0000011 0101010 1010101 0101011 0011111
  # 2-2-0 3-0-0 4-0-0 4-1-0 5-0-0
  # 23445 is the 16th permutation (of 0-59, 13420) of 43425. Undoing the
  # permutation gives 4-0-0 3-0-0 4-1-0 2-2-0 5-0-0.
  # So partA=335 (group 1 starts at 320, plus 15 to get 120), partB=150
  # (the 6 of 16, times 25), and zelPart=16384.
  perminx=permtype=0
  zelCode=0x0000
  cw=[0,0,0,0,0,0]
  @assert partA<625
  if partA<125
    perminx=0
    permtype=10 # pattern 43434
    cw[5]=8 # selects 4 bits set letters
    cw[4]=0 # selects 3 bits set letters
    cw[3]=8+partA÷25
    cw[2]=0+(partA÷5)%5
    cw[1]=8+partA%5
  elseif partA<275
    partA-=125
    perminx=partA÷75
    permtype=20 # pattern 33435
    partA%=75
    cw[5]=0
    cw[4]=0
    cw[3]=8+partA÷15
    cw[2]=0+(partA÷3)%5
    cw[1]=13+partA%3 # 13 selects 5 bits set letters, of which there are 3 columns
  elseif partA<545
    partA-=275
    perminx=partA÷45
    permtype=60 # pattern 43425
    partA%=45
    cw[5]=8
    cw[4]=0
    cw[3]=8+partA÷9
    cw[2]=5+(partA÷3)%3 # 5 selects 2 bits set letters
    cw[1]=13+partA%3
  else
    partA-=545
    perminx=partA÷27
    permtype=30 # pattern 33525
    partA%=27
    cw[5]=0
    cw[4]=0
    cw[3]=13+partA÷9
    cw[2]=5+(partA÷3)%3
    cw[1]=13+partA%3
  end
  perminx=perminx*10+partB÷25
  #println(perminx)
  cw[5]+=(partB÷5)%5
  cw[4]+=partB%5
  if zelPart>=0 && zelPart<16384
    zelCode=zel[zelPart]
  else
    zelCode=0x0000
  end
  for i in 0:4
    cw[i+1]+=((zelCode>>(3*i))&7)<<4
    cw[i+1]=letter[cw[i+1]]⊻(sign*127)
  end
  cw[6]=sign
  cw=Codeword(cw)
  if permtype%3==0
    cw=permute(cw,perm60[perminx])
  else
    cw=permute(cw,perm20[perminx])
  end
  cw
end

function encode(n::Unsigned,bits::Integer)
  # Encodes the following plaintexts:
  # n		bits
  # 32-bit int	32
  # 24-bit int	24
  # 16-bit int	16
  # 8-bit int	8
  # 0 to 1-127	1-7 (TBD)
  # 0 or 1	0 (idle code)
  zelPart=Int(n&0x03fff)
  mumPart=Int((n&0x7fffffff)>>14)
  signBit=Int(n>>31);
  partA=626
  partB=250
  if bits==32
    if signBit!=0
      zelPart⊻=0x3fff	# flipping all bits of the input
      mumPart⊻=0x1ffff	# flips all bits of the output
    end
  elseif bits==24
    mumPart+=0x20000
  elseif bits==16
    mumPart+=0x20400
  elseif bits==8
    mumPart+=0x20404
  elseif bits==0
    mumPart=156500
  else
    mumPart=262144
  end
  if mumPart<156500
    if mumPart&255<224
      if mumPart>>8<256
	partA=(mumPart>>8)+275
      elseif mumPart>>8<0x170
	partA=(mumPart>>8)-256
      else
	partA=(mumPart>>8)+125-0x170
	@assert partA<545
      end
      partB=mumPart&0xff
    else
      partA=(mumPart>>8)÷7+545
      partB=((mumPart>>8)%7)*32+(mumPart&31)
      @assert partA<626
    end
  else
    if isodd(bits)
      partA=335
      partB=150
      signBit=0
    else
      partA=506
      partB=225
      signBit=1
    end
    zelPart=16384
  end
  halfEncode(signBit,partA,partB,zelPart)
end

# Decoding

function halfDecode(cw::Codeword)
  pc=permcode(map(UInt8∘count_ones,cw))
  perminx=letterPerm[pc]
  letterRows=map(x->invLetter[x]>>4,cw)
  noPermute=perminx>0xdf || findfirst(x->x==7,letterRows)<6 # the sign bit turns into 7
  if noPermute
    cwup=cw
  else
    if perminx<0x80
      cwup=permute(cw,invPerm20[perminx&0x3f])
    else
      cwup=permute(cw,invPerm60[perminx&0x3f])
    end
  end
  letterRows=map(x->invLetter[x]>>4,cwup)
  letterColumns=map(x->invLetter[x]&0xf,cwup)
  zelPart=invZel[permoct(letterRows)]
  partB=perminx&63%10*25+(letterColumns[5]&7)*5+(letterColumns[4]&7)
  # partB ranges from 0 to 249, but only 0-223 is valid,
  # except syncword 1, which is 225.
  #println(letterColumns,' ',perminx)
  if perminx<64 # pattern 43434
    partA=(letterColumns[3]&7)*25+(letterColumns[2]&7)*5+(letterColumns[1]&7)
  elseif perminx<128 # pattern 33435
    partA=(letterColumns[3]&7)*15+(letterColumns[2]&7)*3+(letterColumns[1]&7-5)+125
    partA+=(perminx-64)÷10*75
  elseif perminx<192 # pattern 43425
    partA=(letterColumns[3]&7)*9+(letterColumns[2]&7-5)*3+(letterColumns[1]&7-5)+275
    partA+=(perminx-128)÷10*45
  else # pattern 33525
    partA=(letterColumns[3]&7-5)*9+(letterColumns[2]&7-5)*3+(letterColumns[1]&7-5)+545
    partA+=(perminx-192)÷10*27
  end
  signBit=cw[6]
  if zelPart<0x4000 && sum(map(count_ones,cw))!=18
    zelPart=0xffff
  end
  (signBit,partA,partB,zelPart)
end

"""
    decode(cw::Codeword,flip::Bool)

Decodes a `Codeword`. The return value is `(n,bits,codeBits,upsideDown)` where

- `bits`==0: it's a syncword or framing error thereof
- `bits`∈[1,32]: `n` is data
- `bits`>32: `n` is a control code
- `bits`<0: it's an error.

`codeBits` is the number of bits till the next frame (normally 36); `upsideDown` should be xored into the next `flip`.
"""
function decode(cw::Codeword,flip::Bool)
  if flip
    cw=Codeword([cw[1]⊻0x7f,cw[2]⊻0x7f,cw[3]⊻0x7f,cw[4]⊻0x7f,cw[5]⊻0x7f,cw[6]⊻0x1])
  end
  (signBit,partA,partB,zelPart)=halfDecode(cw)
  codeBits=36 # wait 36 bits for next frame, unless it's a framing error
  n=0
  upsideDown=(partA==335 && partB==150 && signBit==1) ||
	     (partA==506 && partB==225 && signBit==0)
  if zelPart==0xffff # it's undecodable
    bits=64
  elseif zelPart>0x6000 # idle code, syncword (including framing error)
    bits=0
    codeBits=Int(zelPart&0xff)
    if partA==335 && partB==150
      n=1
    elseif partA==506 && partB==225
      n=0
    else
      n=65535
    end
  else
    if partA<545
      if partA<125
	mumPart=partA+0x100
      elseif partA<275
	mumPart=partA-125+0x170
      else
	mumPart=partA-275
      end
      mumPart=mumPart*256+partB
    else
      mumPart=(partA-545)*7+partB
    end
    n=mumPart<<14+zelPart
    if isodd(signBit)
      n⊻=0xffffffff
    end
    bits=32
  end
  (n,bits,codeBits,upsideDown)
end

end # module MumzelCode
