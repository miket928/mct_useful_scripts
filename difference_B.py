import os
import sys
import math
import string

pdb1_filename = sys.argv[1]
pdb2_filename = sys.argv[2]

cwd = os.getcwd()

pdb1_fullpath = os.path.join(cwd, pdb1_filename)
pdb2_fullpath = os.path.join(cwd, pdb2_filename)

#remove alt confs

remove_alt_loc_command_pdb1 = "phenix.pdbtools %s remove_alt_confs=True convert_to_isotropic=True output.file_name=pdb1_no_altloc.pdb" % (pdb1_fullpath)
os.system(remove_alt_loc_command_pdb1)

remove_alt_loc_command_pdb2 = "phenix.pdbtools %s remove_alt_confs=True convert_to_isotropic=True output.file_name=pdb2_no_altloc.pdb" % (pdb2_fullpath)
os.system(remove_alt_loc_command_pdb2)

pdb1_no_alt_loc_fullpath = os.path.join(cwd, "pdb1_no_altloc.pdb")
pdb2_no_alt_loc_fullpath = os.path.join(cwd, "pdb2_no_altloc.pdb")

pdb1_input = open(pdb1_no_alt_loc_fullpath, 'r')
pdb2_input = open(pdb2_no_alt_loc_fullpath, 'r')

pdb_output_fullpath = os.path.join(cwd, "difference_B.pdb")
pdb_output = open(pdb_output_fullpath, 'w')

pdb1_atom_lines = []

for line1 in pdb1_input.readlines():
  if line1.startswith('ATOM'):
    pdb1_atom_lines.append(line1)

test_line = pdb1_atom_lines[1000]
test_resnum = test_line[22:25]
test_atom_name = test_line[12:15]
test_B = test_line[61:65]

print pdb1_atom_lines[1000]
print test_line[22:26]
print test_line[13:17]
print test_line[61:66]

pdb2_atom_lines = []

for line2 in pdb2_input.readlines():
  if line2.startswith('ATOM'):
    pdb2_atom_lines.append(line2)

for atomA_line in pdb1_atom_lines:
  residueA_number = atomA_line[22:26]
  atomA_name = atomA_line[13:17]
  B_factorA = float(atomA_line[61:66])
  for atomB_line in pdb2_atom_lines:
    residueB_number = atomB_line[22:26]
    atomB_name = atomB_line[13:17]
    B_factorB = float(atomB_line[61:66])
    if (residueA_number == residueB_number) and (atomA_name == atomB_name):
      difference_B = B_factorA - B_factorB
      output_atom_line = "%s%6.2f%s" % (atomA_line[:60], difference_B, atomA_line[66:])
      pdb_output.write(output_atom_line)
    else:
      pass

pdb1_input.close()
pdb2_input.close()
pdb_output.close()
