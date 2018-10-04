import string
import re
import os

gamma_list = ["0", "0.2", "0.4", "0.6", "0.8", "1.0"]
weight_list = ["3", "10", "30", "100", "300"]

parent_dir = "/netapp/home/mthomp/Proteasome/ALS_20160422/DEN_REFINE/"
pdb_filename = "Proteasome_c6_x1_refine_002.pdb"
mtz_filename = "Proteasome_c6_x1_free.mtz"

directories = []
qsub_lines = []

for gamma in gamma_list:
  for weight in weight_list:

    dir_name = "DEN_ref_gamma_%s_weight_%s" % (gamma, weight)
    os.mkdir(parent_dir+dir_name)
    os.mkdir(parent_dir+dir_name+"/err")
     
    os.system("cp ./"+pdb_filename+" ./"+dir_name)
    os.system("cp ./"+mtz_filename+" ./"+dir_name)

    sub_filename = parent_dir+dir_name+"/DEN_ref_gamma_%s_weight_%s.sub" % (gamma, weight)
    sub_file = open(sub_filename, "w")

    sub_file_text = """#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o %s%s
#$ -e %s%s/err
#$ -r y
#$ -j y
#$ -l mem_free=2G
#$ -l arch=linux-x64
#$ -l netapp=2G,scratch=2G
#$ -l h_rt=24:00:00

source /netapp/home/rwoldeyes/phenix/phenix-1.9-1692/phenix-1.9-1692/phenix_env.sh

phenix.den_refine %s%s/%s %s%s/%s refinement.den.gamma=%s refinement.den.weight=%s refinement.den.optimize=False

""" % (parent_dir, dir_name, parent_dir, dir_name, parent_dir, dir_name, pdb_filename, parent_dir, dir_name, mtz_filename, gamma, weight)

    sub_file.write(sub_file_text)

    qsub_line = "qsub "+sub_filename
    qsub_lines.append(qsub_line)
      
    directories.append(parent_dir+dir_name)      

shell_script = open("submit_DEN_ref_to_cluster.sh", "w")
shell_script.write("#!/bin/bash\n\n")
for i in range(len(qsub_lines)):
  shell_script.write("cd "+directories[i-1]+"\n")
  shell_script.write(qsub_lines[i-1]+"\n")
shell_script.write("cd ..")

os.system("chmod 777 submit_DEN_ref_to_cluster.sh")
