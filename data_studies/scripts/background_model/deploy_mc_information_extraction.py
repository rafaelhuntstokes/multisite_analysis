import numpy as np
import string
import os

def deploy():
    # load up directionality run list
    run_list = np.loadtxt("mc_list.txt", dtype = int)

    isotope     = "BiPo212"
    BATCH_NAME  = f"{isotope}_mc_extraction"
    SUBMIT_PATH = "/data/snoplus3/hunt-stokes/clean_multisite/condor" 

    # loop over the runs and split into each job every NUM_RUNS_PER_JOB
    for irun in run_list:

            # check if the ratds file for this run exists in the data drive
            if os.path.isfile(f"/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ratds/simulation{isotope}_{irun}.root") == False:
                continue
            
            # create the .sh file and .submit file for the job
            with open("extract_mc_information.sh", "r") as infile:
                rawTextSh = string.Template(infile.read())
            with open("extract_mc_information.submit", "r") as infile:
                rawTextSubmit = string.Template(infile.read())

            job_name      = f"{isotope}_mc_extraction_{irun}"
            outTextSh     = rawTextSh.substitute(RUN_NUMBER = irun, ISOTOPE = isotope)
            outTextSubmit = rawTextSubmit.substitute(NAME = job_name)

            with open(f"{SUBMIT_PATH}/sh/{job_name}.sh", "w") as outfile:
                outfile.write(outTextSh)
                os.chmod(f"{SUBMIT_PATH}/sh/{job_name}.sh", 0o0777)
            with open(f"{SUBMIT_PATH}/submit/{job_name}.submit", "w") as outfile:
                outfile.write(outTextSubmit)
            
            # and submit the job
            command = f"condor_submit -b {BATCH_NAME} {SUBMIT_PATH}/submit/{job_name}.submit"
            os.system(command)

def check_result():
    run_list = np.loadtxt("mc_list.txt", dtype = int)
    isotope  = "BiPo212"
    missing  = []
    for irun in run_list:

        # check if the ratds file for this run exists in the data drive
        if os.path.isfile(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/energy/energy_{irun}.npy") == False:
            missing.append(irun)
            continue

    for i in missing:
        print(i)
    print(f"Missing {len(missing)} runs.")
check_result()