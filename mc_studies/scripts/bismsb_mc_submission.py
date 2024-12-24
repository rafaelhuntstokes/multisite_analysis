"""
Script submits the bisMSB macros needed to make the bisMSB Asimov dataset. The macros
are updated to use Po-Wei's/Serena's new optical model.
"""

import string
import os 
import numpy as np
import argparse


path         = "/data/snoplus2/miniPROD_bismsb_newOptics_oldFitter7015" # output directory to place MC in
isotope_name = "B8_solar_numu"                                                     # isotope / macro name (be careful to make it match the macro name!)
evs_per_job  = 5000                                                     

# load up the MC runlist (one run conditions job per run)
runlist = np.loadtxt("../runlists/bismsb_mc_runlist.txt", dtype = int)

# create a macro, sh and submit file per job (1 per run)
for irun in runlist:
    # load in the default macro text 
    with open(f"../condor/bismsb_macros/{isotope_name}.mac", "r") as infile:
        rawTextMac = string.Template(infile.read())
    # load in the default executable .sh 
    with open(f"../condor/template_simulate_bismsb.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())
    with open(f"../condor/template_simulate_bismsb.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())

    # name the macro according to the run number
    macName = f"{isotope_name}_{irun}_bismsb"

    # substitute the macro name, outputs, run number etc. to the template files
    outTextMac = rawTextMac.substitute(OUTPUT_NTUPLE = f"{path}/ntuples/{macName}.ntuple.root", OUTPUT_RATDS = f"{path}/ratds/{macName}.root", NUM_EVENTS = evs_per_job)
    outTextSh = rawTextSh.substitute(RUN_CONDITIONS = irun, MACNAME = macName)
    outTextSubmit = rawTextSubmit.substitute(SHNAME = macName)

    # create the specific macro, sh and submit file for this job
    with open(f"../condor/macros/" + macName + ".mac", "w") as outfile:
        outfile.write(outTextMac)
    with open(f"../condor/sh/" + macName + ".sh", "w") as outfile:
        outfile.write(outTextSh)
        # make it executable 
        os.chmod(f"../condor/sh/" + macName + ".sh", 0o0777)
    with open(f"../condor/submit/" + macName + ".submit", "w") as outfile:
        outfile.write(outTextSubmit)

    # submit the job
    command = f"condor_submit -b {isotope_name}_bismsb ../condor/submit/{macName}.submit"
    os.system(command)