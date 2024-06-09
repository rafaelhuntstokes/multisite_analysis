from htcondor import dags
import htcondor
import numpy as np
import os
import time

def submit_analysis(jobs_per_layer, fv_cut, z_cut):
    """
    Function submits analysis scripts as a sequence of layers, with each layer
    consisting of N nodes (the 'breadth'). This allows me to automate the
    submission of my analysis tagging etc. across the directionality list.
    """

    # load the run list - 1 job / node per run
    runlist = np.loadtxt("../runlists/quiet_period.txt", dtype = int)

    # define the workin directory for log files and suchlike
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"

    # work out the number of layers
    num_layers = int(len(runlist) / jobs_per_layer)

    # add final layer containing the remainder
    if len(runlist) % jobs_per_layer != 0:
        num_layers += 1

    print(f"Creating DAG with {num_layers} layers of {jobs_per_layer} jobs. Final layer has {len(runlist) % jobs_per_layer} jobs.")

    # create a .txt file containing input arguments for each job in each layer
    for ilayer in range(num_layers):
        
        with open(f"../condor/dag_arguments/layer_{ilayer}_args.txt", "w") as argsfile:
            
            # each line is the arguments for a given job
            if ilayer != num_layers - 1:
                for ijob in range(jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{runlist[run_idx]} {fv_cut} {z_cut}\n")
            else:
                for ijob in range(len(runlist) % jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{runlist[run_idx]} {fv_cut} {z_cut}\n")

    
    # create a .submit file for each job in each layer
    for ilayer in range(num_layers):

        with open(f"../condor/dag_submit/layer_{ilayer}.sub", "w") as layerfile:
            layerfile.write(f"executable = {working_dir}/condor/dag_analysis.sh\n")
            layerfile.write(f"output     = {working_dir}/condor/dag_output/layer_{ilayer}_$(Process).out\n")
            layerfile.write(f"error      = {working_dir}/condor/dag_error/layer_{ilayer}_$(Process).err\n")
            layerfile.write(f"log        = {working_dir}/condor/dag_log/layer_{ilayer}_$(Process).log\n")
            layerfile.write("max_retries = 3\n")
            layerfile.write(f"queue arguments from {working_dir}/condor/dag_arguments/layer_{ilayer}_args.txt\n")

    
    # create the .dag file and specify the dependencies between layers
    with open("../condor/analysis_dag.dag", "w") as dagfile:

        for ilayer in range(num_layers):


            dagfile.write(f"JOB BATCH_{ilayer} {working_dir}/condor/dag_submit/layer_{ilayer}.sub\n")
            
            if ilayer != 0:
                dagfile.write(f"PARENT BATCH_{ilayer - 1} CHILD BATCH_{ilayer}\n")

submit_analysis(50, 6000.0, -6000.0)