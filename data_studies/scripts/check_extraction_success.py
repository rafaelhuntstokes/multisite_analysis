import numpy as np
import os
import ROOT
import glob
import rat
from ROOT import RAT
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
                missing.append(irun)
    
    print("Missing Runs: \n", missing)
    np.savetxt("../runlists/failed_analysis.txt", missing, fmt = '%d')

def check_extraction_success():
    """
    This function is to check the ratds extraction of events has succeeded for every
    run in the run list. 
    
    The number of events in the extracted ratds files should equal the number of
    events in the corresponding ntuple files. If this is not the case, something
    has gone wrong - most likely the extraction job has failed on condor and 
    needs to be resubmitted.
    """

    runlist = np.loadtxt("../runlists/subset_dir.txt", dtype = int)

    missing_ntuple   = [] # if the ntuple doesn't exist for a given run
    bad_ratds        = [] # files missing / corrupted / wrong number of events
    missing_ratds    = []
    bad_runs         = [300000, 301140, 301148, 305479, 300121]
    # run in airplane mode
    RAT.DB.Get().SetAirplaneModeStatus(True)
    total_num_ntuple = 0
    total_num_ratds  = 0
    for irun in runlist:
        if irun in bad_runs:
            continue
        print(irun)
        # open the ntuple for this run
        try:
            file_ntuple  = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/output_{irun}.root")
            file_tree    = file_ntuple.Get("clean_6m")
        except:
            missing_ntuple.append(irun)
            file_ntuple.Close()
            continue

        # find number of events inside ROI to compare to RATDS
        num_evs_ntpl = file_tree.GetEntries("energy>=2.5&energy <= 5.0")
        total_num_ntuple += num_evs_ntpl
        # now glob together all the extracted ratds subrun files for this run
        ratds_flist = glob.glob(f"../extracted_data/full_analysis/extracted_ratds/{irun}*.root")
        if len(ratds_flist) == 0:
            # no extraction happened for this run
            print(f"No subruns found for {irun}!")
            bad_ratds.append(irun)
            missing_ratds.append(irun)
            continue
        
        # find the total extracted number of events in the ratds subrun files
        num_evs_ratds = 0
        for iratds in ratds_flist:
            # print(iratds)
            # ds = RAT.DU.DSReader("../extracted_data/full_analysis/extracted_ratds/301125_2.root")
            try:
                file = ROOT.TFile.Open(iratds)
                tree = file.Get("T")
                num_evs_ratds += tree.GetEntries()
            except:
                print("Problem loading file: ", iratds)
                file.Close()
                continue

        file.Close()
        file_ntuple.Close()
            # print(num_evs_ratds)
        total_num_ratds += num_evs_ratds
        if num_evs_ntpl != num_evs_ratds:
            print(f"Mismatch in run: {irun}!")
            print(num_evs_ntpl, num_evs_ratds)
            bad_ratds.append(irun)
            continue

    print("Mismatch numbers: ")
    for irun in bad_ratds:
        print(irun)

    print(f"Problem with {len(bad_ratds)} runs.")

    print("Runs with no ratds files present at all: ")
    for irun in missing_ratds:
        print(irun)
    print("\n")
    print(f"Total events in ntuple: {total_num_ntuple}\nTotal events ratds: {total_num_ratds}")

def check_analysis_success():
    """
    Function checks that the number of events per run in the multisite 
    discriminant analysis files matches the number of events in the .ntuple
    energy spectrum file.
    
    Returns a list of runs missing events / missing output files for 
    re-running and bug checking.
    """

    runlist = np.loadtxt("../runlists/data_analysis.txt", dtype = int)

    missing_ntuple   = [] # if the ntuple doesn't exist for a given run
    bad_ratds        = [] # files missing / corrupted / wrong number of events
    missing_ratds    = []
    bad_runs         = [300000, 301140, 301148, 305479, 300121]
    # run in airplane mode
    RAT.DB.Get().SetAirplaneModeStatus(True)
    total_num_ntuple = 0
    total_num_ratds  = 0
    for irun in runlist:
        if irun in bad_runs:
            continue
        print(irun)
        # open the ntuple for this run
        try:
            file_ntuple  = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/output_{irun}.root")
            file_tree    = file_ntuple.Get("clean_4p5m")
        except:
            missing_ntuple.append(irun)
            file_ntuple.Close()
            continue

        # find number of events inside ROI to compare to RATDS
        num_evs_ntpl = file_tree.GetEntries("energy>=2.5&energy <= 5.0")
        total_num_ntuple += num_evs_ntpl
        
        # find the total extracted number of events in the ratds subrun files
        num_evs_ratds = 0
        try:
            file = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/processed_dataset/{irun}.root")
            tree = file.Get("2p5_5p0")
            num_evs_ratds += tree.GetEntries()
        except:
            print("Problem loading file: ", file)
            bad_ratds.append(irun)
            file.Close()
            continue

        file.Close()
        file_ntuple.Close()

        total_num_ratds += num_evs_ratds

        print(f"Num in ntuple: {num_evs_ntpl}\nNum in dataset: {num_evs_ratds}")

        if num_evs_ntpl != num_evs_ratds:
            print(f"Mismatch in run: {irun}!")
            print(num_evs_ntpl, num_evs_ratds)
            bad_ratds.append(irun)
            continue

    print("Mismatch numbers: ")
    for irun in bad_ratds:
        print(irun)

    print(f"Problem with {len(bad_ratds)} runs.")

    print("Runs with no ratds files present at all: ")
    for irun in missing_ratds:
        print(irun)
    print("\n")
    print(f"Total events in ntuple: {total_num_ntuple}\nTotal events dataset: {total_num_ratds}")

def delete_full_ratds():
    """
    Dangerous function! Assuming you have verified the extraction of ratds files
    succeeded (see function above), you can delete the full ratds files to 
    free up space. This function loads a run list, loops through it and 
    deletes the corresponding ratds files.
    """

    runlist  = np.loadtxt("../runlists/batch4_2_missing.txt", dtype = int)
    bad_runs = [300000, 301140, 301148, 305479, 300121]
    
    for irun in runlist:
        print(irun)
        # keep rat ds with lone orphans as maybe I need them later if I fix these runs?
        if irun in bad_runs:
            continue

        # glob all ratds files for this run
        data_path = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ratds";

        # be mindful of the reprocessing of earlier data
        if irun < 307613:
            inFile1 = data_path + "/" + "Analysis20R_r0000"
        if irun >= 307613:
            inFile1 = data_path + "/" + "Analysis20_r0000"
        inFileName = f"{inFile1}{irun}*.ntuple.root"
        ratds_flist = glob.glob(inFileName)
        print(ratds_flist)
        for isubrun in ratds_flist:
            command = f"rm {isubrun}"
            os.system(command)
            # break
        # break





# data_dir   = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/above_5MeV/solar_run_by_run_ntuples/livetimes"
# data_dir = "/data/snoplus3/hunt-stokes/clean_multisite/bipo_212_cleanSelec3/completed"
# run_list = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/runlists/contains_solar_candidates.txt"
# runs_containing_solar_candidate_ratds(run_list)
# check_analysed_candidate_output()
# run_list = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/runlists/directionality_list.txt"
# check_runs(run_list, data_dir)
# check_extraction_success()
check_analysis_success()
# delete_full_ratds()