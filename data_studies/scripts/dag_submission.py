from htcondor import dags
import htcondor
import numpy as np
import os
import time
import string
import glob
import ROOT

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

def submit_reproc(jobs_per_layer):
    """
    This function creates a DAG that reprocesses all the events in extracted 7.0.8
    RATDS files and reconstructs them with RAT 7.0.15.
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
        
        with open(f"../condor/dag_arguments_reproc/layer_{ilayer}_args.txt", "w") as argsfile:
            
            # each line is the arguments for a given job
            if ilayer != num_layers - 1:
                for ijob in range(jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{working_dir}/condor/dag_sh_reproc/reprocess_{runlist[run_idx]}.sh\n")
                    
            else:
                for ijob in range(len(runlist) % jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{working_dir}/condor/dag_sh_reproc/reprocess_{runlist[run_idx]}.sh\n")


    # create a .submit file for each layer
    for ilayer in range(num_layers):

        with open(f"../condor/dag_submit_reproc/layer_{ilayer}.sub", "w") as layerfile:
            layerfile.write(f"executable = {working_dir}/condor/dag_reproc.sh\n")
            layerfile.write(f"output     = {working_dir}/condor/dag_output_reproc/layer_{ilayer}_$(Process).out\n")
            layerfile.write(f"error      = {working_dir}/condor/dag_error_reproc/layer_{ilayer}_$(Process).err\n")
            layerfile.write(f"log        = {working_dir}/condor/dag_log_reproc/layer_{ilayer}_$(Process).log\n")
            layerfile.write("max_retries = 10\n")
            layerfile.write("request_memory = 2000\n")
            layerfile.write("request_cpus = 1\n")
            layerfile.write(f"queue arguments from {working_dir}/condor/dag_arguments_reproc/layer_{ilayer}_args.txt\n")

    # create the .sh file that calls the macros for each job
    # open the templates
    with open("../condor/reprocessing_template.mac", "r") as infile:
        rawTextMacro  = string.Template(infile.read())
    with open("../condor/reprocessing_template.sh", "r") as infile:
        rawTextSh     = string.Template(infile.read())
    
    data_dir    = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/extracted_ratds"
    output_path = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/reprocessed_ratds_7.0.15"
    
    # loop over the runs with candidates identified in them
    num_submitted = 0
    no_data  = []
    for irun in runlist:
        print(irun)

        # glob all of the subrun files associated with this run
        flist = glob.glob(f"{data_dir}/{irun}*.root")
        print(f"Found {len(flist)} subruns for run {irun}.")
        if len(flist) == 0:
            no_data.append(irun)
            print(".sh file will be empty (no extraction to run) but made anyway so condor doesn't fail.")
        else:
            # we have some subruns here --> check if there's events in those subruns and only add the extraction macro to the sh file if events exist in it
            macro_string = ""
            for ifile in range(len(flist)):
                file = ROOT.TFile.Open(flist[ifile])
                num_evs = file.Get("T").GetEntries()

                if num_evs > 0:
                    
                    
                    # check how many entries there are in the subrun file --> if there aren't any, don't add the subrun!
                    name = f"{irun}_{ifile}"
                    
                    # macro for each subrun
                    outTextMacro  = rawTextMacro.substitute(IN = flist[ifile], OUT = f"{output_path}/{name}.root")
                    with open(f"../condor/dag_macros_reproc/reprocess_{name}.mac", "w") as outfile:
                        outfile.write(outTextMacro)
                    
                    # create a text string containing each macro name to be run and the rat call separated by new lines
                    macro_string += f"rat /data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/dag_macros_reproc/reprocess_{name}.mac\n"
                    
                    # one sh file and one submit file per run
                    outTextSh     = rawTextSh.substitute(MACROS = macro_string)
                    with open(f"../condor/dag_sh_reproc/reprocess_{irun}.sh", "w") as outfile:
                        outfile.write(outTextSh)
                        os.chmod(f"../condor/dag_sh_reproc/reprocess_{irun}.sh", 0o0777)
                else:
                    print(f"No events in {ifile}!")

    # create the .dag file and specify the dependencies between layers
    with open("../condor/reproc_dag.dag", "w") as dagfile:

        for ilayer in range(num_layers):


            dagfile.write(f"JOB BATCH_{ilayer} {working_dir}/condor/dag_submit_reproc/layer_{ilayer}.sub\n")
            
            if ilayer != 0:
                dagfile.write(f"PARENT BATCH_{ilayer - 1} CHILD BATCH_{ilayer}\n")

def submit_bipo_tagging(jobs_per_layer, out_dir):
    """
    Creates a DAG to submit BiPo tagging over entire directionality list (or whatever
    list). Created for the Tl208 constraint systematic study.
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
        
        with open(f"../condor/dag_arguments_bipo/layer_{ilayer}_args.txt", "w") as argsfile:
            
            # each line is the arguments for a given job
            if ilayer != num_layers - 1:
                for ijob in range(jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{runlist[run_idx]} {out_dir}\n")
            else:
                for ijob in range(len(runlist) % jobs_per_layer):
                    run_idx = ilayer * jobs_per_layer + ijob

                    argsfile.write(f"{runlist[run_idx]} {out_dir}\n")

    
    # create a .submit file for each job in each layer
    for ilayer in range(num_layers):

        with open(f"../condor/dag_submit_bipo/layer_{ilayer}.sub", "w") as layerfile:
            layerfile.write(f"executable = {working_dir}/condor/dag_bipo.sh\n")
            layerfile.write(f"output     = {working_dir}/condor/dag_output_bipo/layer_{ilayer}_$(Process).out\n")
            layerfile.write(f"error      = {working_dir}/condor/dag_error_bipo/layer_{ilayer}_$(Process).err\n")
            layerfile.write(f"log        = {working_dir}/condor/dag_log_bipo/layer_{ilayer}_$(Process).log\n")
            layerfile.write("max_retries = 3\n")
            layerfile.write(f"queue arguments from {working_dir}/condor/dag_arguments_bipo/layer_{ilayer}_args.txt\n")

    
    # create the .dag file and specify the dependencies between layers
    with open("../condor/bipo_dag.dag", "w") as dagfile:

        for ilayer in range(num_layers):


            dagfile.write(f"JOB BATCH_{ilayer} {working_dir}/condor/dag_submit_bipo/layer_{ilayer}.sub\n")
            
            if ilayer != 0:
                dagfile.write(f"PARENT BATCH_{ilayer - 1} CHILD BATCH_{ilayer}\n")

submit_bipo_tagging(50, "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/systematic_tagged212")
# submit_reproc(50)
# submit_analysis(50, 6000.0, -6000.0)