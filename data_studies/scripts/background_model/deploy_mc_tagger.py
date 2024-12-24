import numpy as np
import os
import string 

num = "bipo"
isotope = 1

tag_list = np.loadtxt("../../runlists/simulated_list.txt", dtype=int)
# tag_list = np.loadtxt("missing_runs.txt", dtype = int)
for irun in tag_list:

    # if irun not in tag_list:
    #     continue
    os.system(f"./bipo_tagger_mc {irun} {isotope}")
    # read in the templates and create the files to submit run-by-run analysis
    # with open("template.submit", "r") as infile:
    #     rawTextSubmit = string.Template(infile.read())
    # with open("template.sh", "r") as infile:
    #     rawTextSh = string.Template(infile.read())
    
    # # fill in the run number
    # outTextSubmit = rawTextSubmit.substitute(RUN_NUMBER = irun, ISOTOPE = isotope)
    # outTextSh = rawTextSh.substitute(RUN_NUMBER = irun, ISOTOPE = isotope)

    # with open(f"../condor/submit/{irun}_{isotope}.submit", "w") as outfile:
    #     outfile.write(outTextSubmit)
    # with open(f"../condor/sh/{irun}_{isotope}.sh", "w") as outfile:
    #     outfile.write(outTextSh)
    #     # make executable
    #     os.chmod(f"../condor/sh/{irun}_{isotope}.sh", 0o0777)
    
    # command = f"condor_submit -b bipo{isotope}_{num} ../condor/submit/{irun}_{isotope}.submit"
    # os.system(command)
