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

get_bi214_gtids()