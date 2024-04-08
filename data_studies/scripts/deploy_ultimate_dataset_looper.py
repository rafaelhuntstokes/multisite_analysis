import numpy as np
import os
import string 
import time
import ROOT
def deploy_looper():
    """
    This function runs ultimate_dataset_looper.cpp, which:
    
    1. for a given run, creates a TChain from ntuple subruns
    2. performs a general coincidence tag on the data, flagging coincident pairs
    3. extracts all the events in each 0.5 m FV inside energy ROI to TTree
    
    These TTrees are stored in /data/snoplus3/hunt-stokes/clean_multisite/data_general_coincidences
    """
    
    batch_name = f"extract_roi_ntuples"
    tag_list = np.loadtxt("../runlists/now_missing.txt", dtype=int)

    for irun in tag_list:
        
        with open("../condor/extract_roi_ntuples.submit", "r") as infile:
            rawTextSubmit = string.Template(infile.read())
        with open("../condor/extract_roi_ntuples.sh", "r") as infile:
            rawTextSh = string.Template(infile.read())

        # fill in the run number
        outTextSubmit = rawTextSubmit.substitute(RUN_NUMBER = irun)
        outTextSh = rawTextSh.substitute(RUN_NUMBER = irun)

        with open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/submit/{irun}.submit", "w") as outfile:
            outfile.write(outTextSubmit)
        with open(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/sh/{irun}.sh", "w") as outfile:
            outfile.write(outTextSh)
        
        # make executable
        os.chmod(f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/sh/{irun}.sh", 0o0777)

        command = f"condor_submit -b {batch_name} /data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/submit/{irun}.submit"
        os.system(command)

        time.sleep(2)

def verify_result():
    tag_list     = np.loadtxt("directionality_list.txt", dtype=int)
    missing_runs = []
    for irun in tag_list:

        # try and open the livetime file (created at end of looper script)
        try:
            livetime = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/livetime/{irun}.txt")
            # livetime = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs/livetime/{irun}.txt")
        except:
            # file doesn't exist ...
            missing_runs.append(irun)

    
    for imissing in missing_runs:
        print(imissing)
    print(f"\nMissing {len(missing_runs)} runs.")

def calc_total_livetime():

    tag_list       = np.loadtxt("directionality_list.txt", dtype=int)
    total_livetime = 0 # in seconds
    for irun in tag_list:

        try:
            livetime = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/livetime/{irun}.txt")
            # livetime = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs/livetime/{irun}.txt")
            if livetime < 0:
                print(irun, livetime)
                continue
        except:
            print(irun)
            continue

        total_livetime += livetime

    print(f"Total livetime:\nDays: {total_livetime / (60 * 60 * 24)}\nHours: {total_livetime / (60 * 60)}\n")

def quiet_data_spectrum():
    # run_list = np.loadtxt("livetime_list.txt", dtype = int)
    run_list = np.loadtxt("dirty_list.txt", dtype = int)
    undownloaded_runs = np.loadtxt("not_downloaded_runs.txt", dtype = int)
    # just do this for 4.5 m FV
    data_spectrum = []
    num_tagged_bipo212 = 0
    num_tagged_bipo214 = 0
    total_livetime     = 0
    missing_runs = []
    bad_livetime = []
    missing_bipo212 = []
    missing_bipo214 = []
    num_runs = 0
    for irun in run_list:

        if irun in undownloaded_runs:
            continue
        missing_file = False
        # open the full data spectrum file
        try:
            file = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs/output_{irun}.root")
            livetime = np.loadtxt(f"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs/livetime/{irun}.txt")
            ntuple = file.Get("full_4p5m")
            if livetime < 0:
                bad_livetime.append(irun)
                missing_file = True
            else:
                total_livetime += livetime
                

        except:
            missing_runs.append(irun)
            missing_file = True
        if missing_file == False:
            try:
                for entry in ntuple:
                    energy = entry.energy
                    itr    = entry.ITR
                    r      = entry.x**2 + entry.y**2 + entry.z**2

                    if r > 4500*4500:
                        continue
                    
                    if itr > 0.3 or itr < 0.13:
                        continue
                    
                    data_spectrum.append(energy)
            except:
                missing_runs.append(irun)
            # find the number of tagged BiPo214 and 212 in this run number
            try:
                file214 = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/bipo_214_cleanSelec/output_214_{irun}.root")
                ntuple = file214.Get("4p5m")
                num_tagged_bipo214 += ntuple.GetEntries()
            except:
                missing_bipo214.append(irun)
            try:
                file212 = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/clean_multisite/bipo_212_cleanSelec2/output_212_{irun}.root")
                ntuple = file212.Get("4p5m")
                num_tagged_bipo212 += ntuple.GetEntries()
                num_runs += 1
            except:
                missing_bipo212.append(irun)

    # save the energy spectrum of the data
    np.save("./clean_data_spectrum_ITR.npy", data_spectrum)
    print(num_runs)
    print(f"Extracted {len(data_spectrum)} events from 4.5 m clean dataset.")
    print(f"Total livetime: {total_livetime/(60*60*24)} days.")
    print(f"Tagged {num_tagged_bipo214} BiPo214 and {num_tagged_bipo212} BiPo212 in 4.5 m.")
    print("Missing Runs / Corrupted: \n")
    print(missing_runs)
    print("\nBad Livetime: \n")
    print(bad_livetime)
    print("\nMissing BiPo214 Information: \n")
    print(missing_bipo214)
    print("\nMissing BiPo212 Information: \n")
    print(missing_bipo212)

# verify_result()
deploy_looper()
# calc_total_livetime()
# quiet_data_spectrum()