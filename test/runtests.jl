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

using MumzelCode,Test

@test MumzelCode.invLetter[0x65]==0x3a

@test cycleType(Codeword([0,1,2,3,4,0]))==Codeword([1,1,1,1,1,0])
@test cycleType(Codeword([0,4,1,2,3,0]))==Codeword([1,4,4,4,4,0])
@test cycleType(Codeword([0,4,3,2,1,0]))==Codeword([1,2,6,6,2,0])
