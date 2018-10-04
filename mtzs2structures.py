#Script for copying a set of Rfree flags across multiple isomorphous data sets and
#perform rigid-body refinement to initiate rebuild.
#Usage: libtbx.python mtzs2structures.py [isomorphous.pdb] [r_free.mtz]
#You will need to source a phenix installation to run with libtbx.python
#Run from a directory containing all of the derivative MTZ files and the other inputs (native pdb file and mtz with R-free set)

import os
import sys
import math
import string
from iotbx.reflection_file_reader import any_reflection_file
import subprocess

#Take command line input (isomorphous-ish PDB file and an MTZ file with desired FreeR_flags)
#assess MTZ files
#make a directory for each MTZ
cwd = os.getcwd()
isomorphous_input_pdb = sys.argv[1]
isomorphous_input_pdb_fullpath = os.path.join(cwd, isomorphous_input_pdb)
r_free_mtz = sys.argv[2]
r_free_mtz_fullpath = os.path.join(cwd, r_free_mtz)

root_dir = os.getcwd()

original_mtz_files = []
original_mtz_filenames = []
dir_paths = []

for root, dirs, filenames in os.walk('./'):
  for filename in filenames:
    if filename.split('.')[-1] == 'mtz':
      fullpath = os.path.join(cwd, filename)
      original_mtz_files.append(fullpath)
      original_mtz_filenames.append(filename.split('.')[-2])

for i, mtz_filename in enumerate(original_mtz_filenames):
  dir_path = os.path.join(root_dir, mtz_filename)
  dir_paths.append(dir_path)
  os.mkdir(dir_path, 0755)
  

#Run phenix.xtriage
#retrieve unit cell parameters and space group
#create input for reflection_file_editor and run
#create new pdb with modified CRYST1
#set up phenix.refine for rigid body refinement and run

for n, path in enumerate(dir_paths):
  print path
  os.chdir(path)

  cwd = os.getcwd()

  # xtriage_cmd = 'phenix.xtriage ../'+str(original_mtz_filenames[n])+'.mtz'
  # os.system(xtriage_cmd)
  hkl_file = any_reflection_file("../"+str(original_mtz_filenames[n])+'.mtz')
  info = hkl_file.file_content()
  xtals =info.crystals()
  if len(xtals) > 1:
    xtal = xtals[0]
  else:
    xtal =xtals
  space_group = info.space_group_number()
  unit_cell = xtal.unit_cell_parameters()
  uc_a, uc_b, uc_c, uc_al, uc_be, uc_ga = unit_cell

  uc_str = str(unit_cell)[1:-1]

  # xtriage_outfile = open('./logfile.log', 'r')
  # # unit_cell_line = xtriage_outfile.readlines()[20]
  # for line in xtriage_outfile.readlines():
  #   print line
  # unit_cell = unit_cell_line.split('=')[-1]
  # uc_a = unit_cell_line.split(' ')[-6]
  # uc_b = unit_cell_line.split(' ')[-5]
  # uc_c = unit_cell_line.split(' ')[-4]
  # uc_al = unit_cell_line.split(' ')[-3]
  # uc_be = unit_cell_line.split(' ')[-2]
  # uc_ga = unit_cell_line.split(' ')[-1]

  # space_group_line = xtriage_outfile.readlines()[21]
  # space_group = space_group_line.split('"')[-2]


  input_refl_mtz_file = original_mtz_files[n]
  output_mtz_filename = original_mtz_filenames[n]+"_flags.mtz"
  # output_mtz_file = os.path.join(cwd, flags_mtz_filename)
  output_mtz_file = os.path.join(cwd, output_mtz_filename)

  RFE_def_file_fullpath = os.path.join(cwd, "reflection_file_editor.def")
  RFE_def_file = open(RFE_def_file_fullpath, 'w')
 
  RFE_def_text = """show_arrays = False
dry_run = False
verbose = True
mtz_file {
  crystal_symmetry {
    unit_cell = %s
    space_group = %s
    output_unit_cell = None
    output_space_group = None
    change_of_basis = None
    eliminate_invalid_indices = False
    expand_to_p1 = False
    disable_unit_cell_check = True
    disable_space_group_check = False
    eliminate_sys_absent = True
  }
  d_max = None
  d_min = None
  wavelength = None
  output_file = %s
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s
    labels = IOBS,SIGIOBS
    output_labels = IOBS SIGIOBS
    column_root_label = None
    d_min = None
    d_max = None
    output_as = *auto intensities amplitudes_fw amplitudes
    anomalous_data = *Auto merged anomalous
    force_type = *auto amplitudes intensities
    scale_max = None
    scale_factor = None
    remove_negatives = False
    massage_intensities = False
    filter_by_signal_to_noise = None
    add_b_iso = None
    add_b_aniso = 0 0 0 0 0 0
    shuffle_values = False
    reset_values_to = None
  }
  miller_array {
    file_name = %s
    labels = R-free-flags
    output_labels = FreeR_flag
    column_root_label = None
    d_min = None
    d_max = None
    output_as = *auto intensities amplitudes_fw amplitudes
    anomalous_data = *Auto merged anomalous
    force_type = *auto amplitudes intensities
    scale_max = None
    scale_factor = None
    remove_negatives = False
    massage_intensities = False
    filter_by_signal_to_noise = None
    add_b_iso = None
    add_b_aniso = 0 0 0 0 0 0
    shuffle_values = False
    reset_values_to = None
  }
  r_free_flags {
    generate = True
    force_generate = False
    new_label = "FreeR_flag"
    fraction = 0.1
    max_free = 2000
    lattice_symmetry_max_delta = 5
    use_lattice_symmetry = True
    use_dataman_shells = False
    n_shells = 20
    random_seed = None
    extend = True
    old_test_flag_value = None
    export_for_ccp4 = False
    preserve_input_values = True
    warn_if_all_same_value = True
    adjust_fraction = False
    d_eps = 0.0001
    relative_to_complete_set = False
    remediate_mismatches = False
  }
}

""" % (uc_str, space_group, output_mtz_file, input_refl_mtz_file, r_free_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "wrote RFE def file successfully!"

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])

  # prepare input pdb file
  #     pdbtools
  #     cryst1

  cryst1_sg = str(space_group).rjust(10)
  cryst1_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s    \n" % (uc_a, uc_b, uc_c, uc_al, uc_be, uc_ga, space_group)
  # cryst1_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s    \n" % (unit_cell, space_group)
  # cryst1_line = "CRYST1%s%s    \n" % (unit_cell,space_group)

  orig_pdb_file = open(isomorphous_input_pdb_fullpath, 'r')
  output_pdb_file_fullpath = os.path.join(cwd, "coords_in_current_cell.pdb")
  output_pdb_file = open(output_pdb_file_fullpath, 'w')

  output_pdb_lines = []

  for line in orig_pdb_file.readlines():
    if line.startswith('ATOM'):
      output_pdb_lines.append(line)

  output_pdb_file.write(cryst1_line)

  for atom in output_pdb_lines:
    output_pdb_file.write(atom)

  output_pdb_file.close()

  pdbtools_cmd = "phenix.pdbtools %s convert_to_isotropic=True set_b_iso=20.00 remove='not protein' output.file_name='coords_in_current_cell.pdb'" % (output_pdb_file_fullpath)
  os.system(pdbtools_cmd)

  phenix_refine_cmd = "phenix.refine %s coords_in_current_cell.pdb strategy=rigid_body --space_group %s --unit_cell %s" % (output_mtz_file, space_group, uc_str.replace(" ",""))
  os.system(phenix_refine_cmd)

  os.chdir('..')


