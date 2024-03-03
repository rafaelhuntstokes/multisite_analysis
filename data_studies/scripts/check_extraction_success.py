import numpy as np
import os
import ROOT
import glob

"""
Utility scripts to check solar candidate extraction worked.

check_runs()
    Identify any failed condor jobs to be re-submitted, negative livetime runs
    ('bad' runs) to be excluded.

runs_containing_solar_candidate()
    Create a run list of runs containing solar candidates, to download the 
    ratds files of and run the extraction on.

    Also return the total number of solar candidates found.

    Also create a list of GTIDs of solar candidates per run for extraction 
    script.    
"""

def check_runs(run_list, data_dir):
    """
    Simple utility script opens the livetime .txt file created by solar extraction
    code for each run in the directionality list.

    If the .txt file exists:
        - check if livetime > 0
            if not:
                add to broken runs list
    if .txt file doesn't exist:
        - condor job probably failed (.txt creation last part of extraction code)
        - add run to failed run list for resubmission

    print list of failed jobs for resubmission
    print list of livetime <= 0 runs for ignoring (e.g. missing tables etc.)
    """
    
    
    bad_runs       = []             # negative livetime runs
    missing_runs   = []             # no livetime .txt file (condor failed)
    total_livetime = 0
    livetime_old   = 0
    now_missing = np.loadtxt("../runlists/now_missing.txt", dtype = int)
    counter    = 1
    run_list   = np.loadtxt(run_list, dtype = int)
    total_runs = len(run_list)
    difference = 0 
    for irun in run_list:

        print(f"Loaded Run {counter} / {total_runs}.")

        # check the file exists
        if os.path.exists(f"{data_dir}/{irun}.txt"): #and os.path.exists(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/livetime/{irun}.txt"):
            
            # file exists --> check if livetime is > 0
            time = np.loadtxt(f"{data_dir}/{irun}.txt")
            # time_old = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/livetime/{irun}.txt")

            
            if time < 0:
                # have some broken tables and the run is not good for use
                bad_runs.append(irun)
            else:
                # add run livetime to total livetime
                total_livetime += time
                # livetime_old += time_old
                # if time_old != time:
                    # print(f"Old: {time_old} | New: {time}")
                    # difference += abs(time - time_old)
                    # print("Difference: ", difference)

        else:
            # job probably failed on condor and needs to be resubmitted
            missing_runs.append(irun)
        counter += 1
    print(len(missing_runs))
    print("Missing Runs: \n", missing_runs)
    print("Bad Runs: \n", bad_runs)
    # print(f"Total Livetime: {total_livetime / (60 * 60 * 24)} days.")

    print(f"Total livetime:\nDays: {total_livetime / (60 * 60 * 24)}\nHours: {total_livetime / (60 * 60)}\n")
    # print(f"Total livetime OLD:\nDays: {livetime_old / (60 * 60 * 24)}\nHours: {livetime_old / (60 * 60)}\n")

def runs_containing_solar_candidate_ntuple(run_list, data_dir):
    """
    Identify runs containing a solar candidate in the ntuples.
    
    Output the total number of solar candidates found, and create a run_list
    of runs containing at least 1 solar candidate.
    
    This run list is then used to download corresponding ratds files.
    
    Finally, each run with a solar candidate has a corresponding .txt file
    containing the GTID(s) of the solar candidates, for use with the
    event extraction/pruning macros.
    """

    run_contains_solar = []
    total_candidates   = 0
    ignore             = [300000, 301140, 301148, 305479] # orphans, negative livetimes identified in the check_runs() func
    
    for irun in run_list:
        gtids = [] # to write to run specific extraction list
        
        # open the corresponding ntuple
        if os.path.exists(f"{data_dir}/{irun}.root") and irun not in ignore:
            # ntuple exists! open it and find how many entries it has
            file = ROOT.TFile.Open(f"{data_dir}/{irun}.root")

            tree      = file.Get("tagged_solar")
            num_solar = tree.GetEntries()

            if num_solar > 0:
                run_contains_solar.append(irun)
                total_candidates += num_solar

                # loop over the entries and extract the event's GTID, runID
                for ientry in tree:
                    gtids.append(ientry.gtid)
                
                # write the run by run specific gtid list file
                np.savetxt(f"../extracted_data/above_5MeV/run_by_run_ntuples/gtid_lists/{irun}.txt", gtids, fmt = '%d')
            file.Close()

    print("Runs containing solar candidates: \n", run_contains_solar)
    print(f"Total solar candidates: {total_candidates}")

    # create run list .txt file
    np.savetxt("../runlists/contains_solar_candidates.txt", run_contains_solar, fmt = '%d')

def runs_containing_solar_candidate_ratds(run_list):
    """
    Identify runs containing a solar candidate in the ntuples.
    
    Output the total number of solar candidates found, and create a run_list
    of runs containing at least 1 solar candidate.
    """

    missing_solar          = []
    file_no_solar          = []
    total_candidates       = 0
    ignore                 = [300000, 301140, 301148, 305479] # orphans, negative livetimes identified in the check_runs() func
    run_list = np.loadtxt(run_list, dtype = int)
    for irun in run_list:
        num_in_run = 0
        # glob all the subruns into a filelist to check
        flist = glob.glob(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/extracted_solar_ratds/{irun}*.root")
        
        # loop over each subrun
        for isub in flist:
            # open the corresponding ntuple
            if os.path.exists(isub):
                # ratds exists! open it and find how many entries it has
                file = ROOT.TFile.Open(isub)

                tree      = file.Get("T")
                num_solar = tree.GetEntries()

                if num_solar > 0:
                    num_in_run += 1
                    # print(num_solar)
                    total_candidates += num_solar
                else:
                    file_no_solar.append(isub)
                file.Close()
        print(f"Run {irun} Contains: {num_in_run}")
        if num_in_run == 0:
            missing_solar.append(irun)

    # print("Runs containing solar candidates: \n", files_containing_solar)
    print(f"Total solar candidates: {total_candidates}")
    print(f"Missing solar candidates: \n", missing_solar)
    
    # create run list .txt file
    np.savetxt("../runlists/missing_extracted_solar_candidates.txt", missing_solar, fmt = "%d")

    # create a list of the subrun files which do not contain solar candidate
    print(f"{len(file_no_solar)} Files to Remove: \n", file_no_solar)

    # now remove the files...
    for ifile in file_no_solar:
        os.system(f"rm {ifile}")
# run_list   = np.loadtxt("../runlists/directionality_list.txt", dtype = int)
# data_dir   = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/run_by_run_ntuples/livetimes"
# check_runs(run_list, data_dir)

def check_analysed_candidate_output():
    """
    For each run identified as having a solar candidate, ensure the output analysed
    TTree exists.
    """

    run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    data_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/solar_data_discriminants"
    missing  = []
    for irun in run_list:
        # open the corresponding ntuple
        print(irun)
        if os.path.exists(f"{data_dir}/{irun}.root") == False:
            # ntuple missing!
            missing.append(irun)
        else:
            
            file = ROOT.TFile.Open(f"{data_dir}/{irun}.root")
            
            ntuple = file.Get("5p0_plus")
            try:
                num = ntuple.GetEntries()
            except:
                missing.append(irun)
            if num == 0:
                missing.append
    
    print("Missing Runs: \n", missing)
    np.savetxt("../runlists/failed_analysis.txt", missing, fmt = '%d')


# data_dir   = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/solar_run_by_run_ntuples/livetimes"
data_dir = "/data/snoplus3/hunt-stokes/clean_multisite/bipo_212_cleanSelec3/completed"
# run_list = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/runlists/contains_solar_candidates.txt"
# runs_containing_solar_candidate_ratds(run_list)
# check_analysed_candidate_output()
run_list = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/runlists/directionality_list.txt"
check_runs(run_list, data_dir)