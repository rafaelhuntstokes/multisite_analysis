import os
import string
import numpy as np
import time
import glob

"""
Submission scripts handle submitting analysis and pdf making scripts to condor.
"""

def submit_analysis(fv_cut, z_cut):
    # define the inputs and the outputs
    # input_fpath  = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/bi214/bismuth214_extracted_ratds"
    # input_fpath  = "/data/snoplus3/inacio/dataTaggedBiPo214_20May2022"
    input_fpath = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ratds"
    output_fpath = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis"
    # output_fpath = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/bi214/bismuth214_data_discriminants"

    # path to save the individual sh and submit files to
    condor_path = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor"
    
    # load up the pdf run list
    # run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    # run_list = np.loadtxt("../runlists/failed_analysis.txt", dtype = int)
    run_list = np.loadtxt("../runlists/may20th_bipo_runs.txt", dtype = int)

    # load up the templates
    with open("../condor/analysis_template.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())
    with open("../condor/analysis_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())

    for irun in run_list:
        outTextSh = rawTextSh.substitute(RUN_NUMBER = irun, FV_CUT = fv_cut, Z_CUT = z_cut, INPUT = input_fpath, OUTPUT = output_fpath)
        
        name = f"analysis_{irun}"
        
        outTextSubmit = rawTextSubmit.substitute(SH_NAME = f"{condor_path}/sh/{name}.sh", LOG_NAME = name)

        # create the files
        with open(f"{condor_path}/sh/{name}.sh", "w") as outfile:
            outfile.write(outTextSh)
        os.chmod(f"{condor_path}/sh/{name}.sh", 0o0777)

        with open(f"{condor_path}/submit/{name}.submit", "w") as outfile:
            outfile.write(outTextSubmit)
        
        # submit this pdf extraction
        command = f"condor_submit -b data_{fv_cut}_{z_cut} {condor_path}/submit/{name}.submit"
        os.system(command)
        
def submit_solar_identification(fv_cut, z_cut):
    """
    Code submits a job per run on the directionality run list to identify solar candidates
    in nutples. 
    """
    # define the inputs and the outputs
    input_fpath  = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples"
    output_fpath = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/run_by_run_ntuples"

    # path to save the individual sh and submit files to
    condor_path = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor"
    
    # load up the pdf run list
    run_list = np.loadtxt("../runlists/subset_dir.txt", dtype = int)

    # load up the templates
    with open("../condor/identify_candidates_template.sh", "r") as infile:
        rawTextSh = string.Template(infile.read())
    with open("../condor/identify_candidates_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())
    
    run = 0
    for irun in run_list:
        run = irun
        outTextSh = rawTextSh.substitute(RUN_NUMBER = irun, FV_CUT = fv_cut, Z_CUT = z_cut, INPUT = input_fpath, OUTPUT = output_fpath)
        
        name = f"extraction_{irun}"
        
        outTextSubmit = rawTextSubmit.substitute(SH_NAME = f"{condor_path}/sh/{name}.sh", LOG_NAME = name)

        # create the files
        with open(f"{condor_path}/sh/{name}.sh", "w") as outfile:
            outfile.write(outTextSh)
        os.chmod(f"{condor_path}/sh/{name}.sh", 0o0777)

        with open(f"{condor_path}/submit/{name}.submit", "w") as outfile:
            outfile.write(outTextSubmit)
        
        # submit this pdf extraction
        command = f"condor_submit -b identify_{fv_cut}_{z_cut} {condor_path}/submit/{name}.submit"
        os.system(command)
        time.sleep(1)
    print("Submit up to and including run: ", run)

def submit_ratds_extraction(candidate_runlist, gtid_path, output_path, batch_name):
    """

    INPUTS : candidate_runlist - path to .txt file of runs containing candidate
             gtid_path         - path to folder etc. containing .txt gtids for each run in candidate_runlist
             output_path       - path to folder to save the extracted events
             batch_name        - label for the condor batch of jobs
    For every run identified with a candidate in it [candidate_runlist], create a macro, sh and submit job
    per subrun in that run.
    
    Input to the macro: subrun ratds file to search through
                        gtid .txt list of candidates in that run to search for
                        output path for extracted events
    """

    # load the runlist
    candidate_runlist = np.loadtxt(candidate_runlist, dtype= int)

    # open the templates
    with open("../condor/extract_candidates_template.mac", "r") as infile:
        rawTextMacro  = string.Template(infile.read())
    with open("../condor/extract_candidates_combined_template.sh", "r") as infile:
        rawTextSh     = string.Template(infile.read())
    with open("../condor/extract_candidates_template.submit", "r") as infile:
        rawTextSubmit = string.Template(infile.read())

    # path to snoplus data
    data_dir = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ratds"
    # loop over the runs with candidates identified in them
    num_submitted = 0
    bad_runs = [300000, 300001, 300002, 301140, 301148, 305479]
    for irun in candidate_runlist:
        print(irun)
        if irun in bad_runs:
            continue
        # define the path to the GTID list
        gtids = f"{gtid_path}/{irun}.txt"

        # glob all of the subrun files associated with this run
        # be careful to get analysis 20 OR analysis 20R runs
        # some runs are downloaded with both and we don't want
        # duplicated events in the extraction
        if irun < 307613:
            flist = glob.glob(f"{data_dir}/Analysis20R_r0000{irun}*.root")
        else:
            flist = glob.glob(f"{data_dir}/Analysis20_r0000{irun}*.root")
        print(flist)
        print(f"Found {len(flist)} subruns for run {irun}.")

        # loop over each subrun and create the macro, sh and submit file for the job
        macro_string = ""
        for isub in range(len(flist)):
            name = f"{irun}_{isub}"
            
            # macro for each subrun
            outTextMacro  = rawTextMacro.substitute(IN = flist[isub], LIST = gtids, OUT = f"{output_path}/{name}.root")
            with open(f"../condor/macros/extract_{name}.mac", "w") as outfile:
                outfile.write(outTextMacro)
            
            # create a text string containing each macro name to be run and the rat call separated by new lines
            macro_string += f"rat /data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/macros/extract_{name}.mac\n"
        
        # one sh file and one submit file per run
        outTextSh     = rawTextSh.substitute(MACROS = macro_string)
        outTextSubmit = rawTextSubmit.substitute(SH_NAME = f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/sh/extract_{irun}.sh", LOG_NAME = f"extract_{irun}")
        
        with open(f"../condor/sh/extract_{irun}.sh", "w") as outfile:
            outfile.write(outTextSh)
            os.chmod(f"../condor/sh/extract_{irun}.sh", 0o0777)
        with open(f"../condor/submit/extract_{irun}.submit", "w") as outfile:
            outfile.write(outTextSubmit)

        command = f"condor_submit -b extract_{batch_name} ../condor/submit/extract_{irun}.submit"
        os.system(command)
        num_submitted += 1
        time.sleep(2)
        if num_submitted == 500:
            print("submitted up to and including: ", irun)
            break
# submit_analysis(5250.0, -6000.0)
submit_ratds_extraction("../runlists/subset_dir.txt", "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/gtids", f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/extracted_ratds", "full_dataset_ratds_extraction")