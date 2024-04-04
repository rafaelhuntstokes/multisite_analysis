import ROOT
import numpy as np

"""
Script opens up the run-by-run BiPo2124 tagged NTUPLEs and creates a .txt file of GTIDs of Bi and Po events.
"""

def get_bi214_gtids():
    run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    for irun in run_list:
        
        BI_IDS = []        
        # open ntuple file
        try:
            TFILE = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_bipo/tagged_ntuples/run_by_run214/output_{irun}.root")
        except:
            print(f"Missing run {irun}")
            continue

        # get the ntuple inside the file
        ntuple = TFILE.Get("ntuple214")

        for ientry in ntuple:
            # for every tagged BiPo event, extract the GTID of the Bi214 and Po214 separately
            bi_id = ientry.gtid2
            BI_IDS.append(bi_id)

        # save these gtid list for the run
        np.savetxt(f"../extracted_data/bi214/gtid_lists/{irun}.txt", BI_IDS, fmt = "%i")
        print('Completed ', irun)

def get_dataset_gtids():
    """
    Get all the GTIDS for the full dataset within the ROI for extraction from the full RATDS files.
    """

    run_list = np.loadtxt("../runlists/directionality_list.txt", dtype = int)
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
        
        # now have all the gtids in 6 m not removed by coincidence tags and falling inside ROI --> write them to a file
        np.savetxt(f"../extracted_data/full_analysis/gtids/{irun}.txt", gtids, fmt = "%i")

    print(f"Missing {len(missing)} runs: \n")
    print(missing)
get_dataset_gtids()
# get_bi214_gtids()