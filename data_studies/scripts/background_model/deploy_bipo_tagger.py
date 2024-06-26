import numpy as np
import os
import string 
import time

num = "bipo"
isotope = 212
batch_name = f"{isotope}_data_tags"
tag_list = np.loadtxt("../../runlists/bipo_runs.txt", dtype=int)
start_count = 1
end_count   = 1000
counter = 1
for irun in tag_list:
    
    with open("../../condor/bipo_tagger_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())
    with open("../../condor/bipo_tagger_template.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())

    # fill in the run number
    outTextSubmit = rawTextSubmit.substitute(RUN_NUMBER = irun, ISOTOPE = isotope)
    outTextSh = rawTextSh.substitute(RUN_NUMBER = irun, ISOTOPE = isotope)

    with open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/submit/{irun}_{isotope}.submit", "w") as outfile:
        outfile.write(outTextSubmit)
    with open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/sh/{irun}_{isotope}.sh", "w") as outfile:
        outfile.write(outTextSh)
    
    # make executable
    os.chmod(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/sh/{irun}_{isotope}.sh", 0o0777)

    command = f"condor_submit -b {batch_name} /data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/submit/{irun}_{isotope}.submit"
    os.system(command)

    time.sleep(2)
    

