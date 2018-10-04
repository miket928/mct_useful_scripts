################################################################################
# This script is for Z-score normalization of atomic B-factors. It is useful for
# comparing structures determined from different data sets. For example, you may
# want to compare structures of the same protein determined at different 
# resolutions, +/- ligand, different mutants, etc., or you may want to compare 
# structures that have some other fundamental difference in data quality. Also 
# can be useful for comparing B-factors of homologous structures.

# USAGE: python calculate_Bz.py pdb_filename.pdb

# Output is a PDB file with Bz in the B-factor column.

# NOTE: Z-scores can be negative, while B-factors cannot. If coloring by 
#       B-factor, make sure to include negative values in the spectrum range.
################################################################################

import math 
import string
import sys
import numpy as np

pdb_filenname = str(sys.argv[1])
filename_prefix = pdb_filenname[:-4]
output_filename = "%s_Bz.pdb" % (filename_prefix)

inputfile = open(pdb_filenname, "r")
outputfile = open(output_filename, "w")


atom_lines_list = []
B_factors_list = []

for line in inputfile.readlines():
  if line.startswith("ATOM"):
    line = line.strip('\n')
    atom_lines_list.append(line)
    B_factors_list.append(float(line[60:66]))

B_factors_array = np.array(B_factors_list)

Average_B = np.mean(B_factors_list)
StdDev_B = np.std(B_factors_list)

B_z_list = []

for value in B_factors_array:
  B_z = (value - Average_B)/StdDev_B
  B_z_list.append(B_z)


for i in range(len(atom_lines_list)):
  atom_line = atom_lines_list[i]
  beginning = atom_line[0:60]
  end  = atom_line[66:80]
  new_atom_line = "%s%6.2f%s\n" % (beginning, B_z_list[i], end)
  outputfile.write(new_atom_line)

 



