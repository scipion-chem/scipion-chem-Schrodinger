# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys, math
from schrodinger import structure

def getCloseRes(st, center, radius):
  closeResidues = []
  for at in st.atom:
    coords = at.xyz
    if math.dist(center, coords) <= radius:
      atRes = at.getResidue().__str__()
      if atRes not in closeResidues:
        closeResidues.append(atRes)

  return closeResidues

if __name__ == "__main__":

  maeFile, radius = sys.argv[1:3]
  center = eval(sys.argv[3])
  outFile = sys.argv[4]
  st = structure.StructureReader.read(maeFile)

  closeResidues, i = [], 0
  while len(closeResidues) == 0:
    closeResidues = getCloseRes(st, center, float(radius)+i)
    i += 1

  with open(outFile, 'w') as f:
    f.write(f'Residues in {maeFile} closer to {center} than {radius}\n')
    f.write(','.join(closeResidues))

