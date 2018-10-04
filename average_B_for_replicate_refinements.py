import string
import os
import numpy
import sys


def create_pdb_file_list( starting_directory, n_runs, pdb_filename_prefix ):

  pdb_file_list = []
  n_runs = int(n_runs)
  for i in range(n_runs):
    run_index = str(i+1)
    pdb_file_fullpath = "%s/run%s/%s_run%s_refine_001.pdb" % (starting_directory, run_index, pdb_filename_prefix, run_index)
    pdb_file_list.append(pdb_file_fullpath)
  return pdb_file_list


def average_B_per_atom_as_csv_pdb( pdb_file_list, output_prefix, output_directory ):

  coordinate_set = pdb_file_list[0]
  coordinate_set_input = open(coordinate_set, 'r')
  lines_before_B_factor = []
  lines_after_B_factor =[]
  for coordinate_line in coordinate_set_input.readlines():
    if coordinate_line.startswith('ATOM'):
      lines_before_B_factor.append(coordinate_line[:60])
      lines_after_B_factor.append(coordinate_line[66:])

  all_b_factors = []

  for pdb_file in pdb_file_list:
    pdb_input = open(pdb_file, 'r')
    list_of_b_factors_by_atom = []
    for line in pdb_input.readlines():
      if line.startswith('ATOM'):
        b_factor = float(line[60:66])
        list_of_b_factors_by_atom.append(b_factor)
    all_b_factors.append(list_of_b_factors_by_atom)

  output_csv_filename = "%s_average_B_per_atom.csv" % (output_prefix)
  output_csv_fullpath = os.path.join(output_directory, output_csv_filename)
  output_csv = open(output_csv_fullpath, 'w')

  output_pdb_filename = "%s_average_B_per_atom.pdb" % (output_prefix)
  output_pdb_fullpath = os.path.join(output_directory, output_pdb_filename)
  output_pdb = open(output_pdb_fullpath, 'w')


  for i in range(len(all_b_factors[0])):
    atom_number = i+1
    b_factors_for_atom = []
    for n in range(len(all_b_factors)):
      b_factors_for_atom.append(all_b_factors[n][i])
    mean_b_factor_for_atom = numpy.mean(b_factors_for_atom)
    B_factor_stdv_for_atom = numpy.std(b_factors_for_atom)
    output_csv_line = "%s,%s,%s\n" % (atom_number, mean_b_factor_for_atom, B_factor_stdv_for_atom)
    output_csv.write(output_csv_line)
    output_pdb_line = "%s%6.2f%s" % (lines_before_B_factor[i], mean_b_factor_for_atom, lines_after_B_factor[i])
    output_pdb.write(output_pdb_line)

  output_csv.close()
  output_pdb.close()


def main():

  prefix = sys.argv[1]
  n = sys.argv[2]
  cwd = os.getcwd()

  start_message = "Averaging atomic B-factors across %s replicate refinement runs..." % (n)
  print start_message

  pdb_list = create_pdb_file_list(cwd, n, prefix)
  average_B_per_atom_as_csv_pdb(pdb_list, prefix, cwd)

  end_message = "Averaging Complete!! CSV and PDB files created."
  print end_message


main()
