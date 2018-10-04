# Requires the ProDy library
# Usage: python calculate_RMSD_for_replicate_refinements.py pdb_filename_prefix
# this script assumes 10 replicate refinement runs, but this can be modified in main()

import string
import math
import os
import numpy
import sys
from prody import *


def create_pdb_file_list( starting_directory, n_runs, pdb_filename_prefix ):

  pdb_file_list = []
  n_runs = int(n_runs)
  for i in range(n_runs):
    run_index = str(i+1)
    pdb_file_fullpath = "%s/run%s/%s_run%s_refine_001.pdb" % (starting_directory, run_index, pdb_filename_prefix, run_index)
    pdb_file_list.append(pdb_file_fullpath)

  return pdb_file_list


def initialize_prody_working_directory( prody_directory_name ):

  cwd = os.getcwd()
  prody_directory_fullpath = os.path.join(cwd, prody_directory_name)
  os.mkdir(prody_directory_fullpath)
  os.chdir(prody_directory_fullpath)
  
  return prody_directory_fullpath


def copy_pdbs_to_prody_directory( pdb_file_list, prody_directory_fullpath ):

  for pdb_file in pdb_file_list:
    copy_command = "cp %s %s" % (pdb_file, prody_directory_fullpath)


def collect_and_align_structures( pdb_file_list ):

  structure_group = Ensemble()
  for pdb_file in pdb_file_list:
    pdb_filename = pdb_file.split('/')[-1]
    coordset = parsePDB(pdb_filename)
    protein_atoms = coordset.select('protein')
    structure_group.addCoordset(protein_atoms)
  aligned_structures = structure_group.iterpose(rmsd=0.0001)

  return aligned_structures


def calculate_RMSD( aligned_structures ):

  #aligned_structures is a prody Ensemble()

  overall_RMSD = aligned_structures.getRMSDs()

  return overall_RMSD


def calculate_RMSD_per_residue_csv( aligned_structures, output_csv_prefix ):

  #aligned_structures is a prody Ensemble()

  cwd = os.getcwd()
  n_structures = aligned_structures.numCoordsets()
  n_residues = aligned_structures[0].numResidues()

  output_filename = "%s_per_residue_RMSD.txt" % (output_csv_prefix)
  output_file_path = os.path.join(cwd, output_filename)
  output_file = open(output_file_path, 'w')

  for i in len(residues):
    residue_ensemble = Ensemble()
    residue_number = i+1
    for j in len(n_structures):
      structure = aligned_structures[j]
      residue = structure['A', residue_number]
      residue_ensemble.addCoordset(residue)
    residue_RMSD = residue_ensemble.getRMSDs()

    output_line = "%s,%s\n" % ()
    output_file.write(output_line)

  output_file.close()


def main():

  starting_directory = os.getcwd()
  n_runs = 10
  pdb_filename_prefix = sys.argv[1]
  prody_directory_name = "ProDy_RMSD"

  pdb_file_list = create_pdb_file_list( starting_directory, n_runs, pdb_filename_prefix )
  prody_directory_fullpath = initialize_prody_working_directory( prody_directory_name )
  copy_pdbs_to_prody_directory( pdb_file_list, prody_directory_fullpath )
  aligned_structures = collect_and_align_structures( pdb_file_list )
  overall_rmsd = calculate_RMSD( aligned_structures )

  output_to_screen = """

*******************************************************************************

The all-atom RMSD (Angstrom) for %s input structures is:

%5.3f

******************************************************************************* 

""" % (n_runs, overall_rmsd)

  calculate_RMSD_per_residue_csv( aligned_structures, pdb_filename_prefix )

