import ROOT
import numpy as np

"""
Script opens up the run-by-run BiPo2124 tagged NTUPLEs and creates a .txt file of GTIDs of Bi and Po events.
"""

def get_bi214_gtids():
    run_list = np.loadtxt("../runlists/quiet_period.txt", dtype = int)
    runs_present = []
    num_events   = 0
    for irun in run_list:
        
        BI_IDS = []        
        # open ntuple file
        try:
            TFILE = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/better_tagging/raw_tagging_output/{irun}.root")
        except:
            print(f"Missing run {irun}")
            continue

        # get the ntuple inside the file
        ntuple = TFILE.Get("tag_info")

        for ientry in ntuple:
            # for every tagged BiPo event, extract the GTID of the Bi214 and Po214 separately
            bi_id = ientry.gtid_bi
            BI_IDS.append(bi_id)
        num_events += len(BI_IDS)
        runs_present.append(irun)
        # save these gtid list for the run
        np.savetxt(f"../extracted_data/better_tagging/bi214_gtids/{irun}.txt", BI_IDS, fmt = "%i")

    np.savetxt("../runlists/better_tagged_bi214.txt", runs_present, fmt  = "%i")

    print(f"Found {len(runs_present)} runs totalling {num_events} events.")
def get_dataset_gtids():
    """
    Get all the GTIDS for the full dataset within the ROI for extraction from the full RATDS files.
    """

    run_list = np.loadtxt("../runlists/batch4_1_missing.txt", dtype = int)
    missing  = []
    for irun in run_list:
        gtids = []

        # open the ntuple files for the dataset
        try:
            file = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/output_{irun}.root")
            ntuple = file.Get("clean_6m") # get the events inside the ROI for all FVs (incase I want to change this later on) with coincidence cuts ON
        except:
            missing.append(irun)
            continue

        for ientry in ntuple:
            energy = ientry.energy
            gtid   = ientry.gtid
            if energy >= 2.5 and energy <= 5.0:
                gtids.append(gtid)
        print(f"Found {len(gtids)} for run {irun}.")
        # now have all the gtids in 6 m not removed by coincidence tags and falling inside ROI --> write them to a file
        np.savetxt(f"../extracted_data/full_analysis/gtids/{irun}.txt", gtids, fmt = "%i")

    print(f"Missing {len(missing)} runs: \n")
    print(missing)
# get_dataset_gtids()
get_bi214_gtids()