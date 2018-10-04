#uses some simple prody functions to calculate the Ca RMSD between two aligned structures
#you need to do the alignment in some other program, like pymol
#input filenames are hard coded - i should fix this

from prody import * 

structure1 = parsePDB("/mnt/c/Users/mctuc/Documents/Papers/MultiDataSet_PCA/finish_refine/group1_aligned.pdb")
structure2 = parsePDB("/mnt/c/Users/mctuc/Documents/Papers/MultiDataSet_PCA/finish_refine/group2_aligned.pdb")
n_res = 215
chains = ['A', 'B', 'C']

outfile = open("RMSD.csv", 'w')

for chain in chains:
  for i in range(215):
    selection = "chain "+chain+" ca resnum "+str(i+1)
    ca1 = structure1.select(selection)
    ca2 = structure2.select(selection)
    rmsd = calcRMSD(ca1, ca2)
    output_line = "%s, %i, %6.3f\n" % (chain, i+1, rmsd)
    outfile.write(output_line)
