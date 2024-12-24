import numpy as np
import matplotlib.pyplot as plt
import ROOT

"""
Script to work out what fraction of events falls into each energy ROI based
on number of events in PDFS/MC.
"""

run_list     = np.loadtxt("../runlists/full_test.txt", dtype = int)
isotope      = "Tl208"
fv_cut       = 4500
mc_path      = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/full_analysis3_"
energies     = []
multisites   = []
missing_runs = [] # runs that don't exist / corrupted but included in run list
print(f"{mc_path}{isotope}/")
for irun in run_list:

    try:
        # open the MC ntuple file
        file   = ROOT.TFile.Open(f"{mc_path}{isotope}/{irun}.root")
        ntuple = file.Get("2p5_5p0")
        ntuple.GetEntries()
    except:
        missing_runs.append(irun)
        continue

    # loop over every event inside file
    for ientry in ntuple:
        
        # check event falls inside FV cut
        x = ientry.x
        y = ientry.y
        z = ientry.z
        r = np.sqrt(x**2 + y**2 + z**2)
        if r > fv_cut:
            continue
        itr = ientry.itr
        if itr < 0.22 or itr > 0.3:
            continue
        posFOM = ientry.posFOM / ientry.posFOM_hits
        if posFOM < 13.7:
            continue
        # event passed FV so add info to PDFs
        energies.append(ientry.energy)
        multisites.append(ientry.dlogL)
print(f"Found {len(missing_runs)} missing runs:\n{missing_runs}")

# convert to numpy array and count
energies = np.array(energies)
print(f"Fraction in 2.5 --> 3.0: {len(energies[(energies >= 2.5) & (energies < 3.0)]) / len(energies)}")
print(f"Fraction in 3.0 --> 3.5: {len(energies[(energies >= 3.0) & (energies < 3.5)]) / len(energies)}")
print(f"Fraction in 3.5 --> 4.0: {len(energies[(energies >= 3.5) & (energies < 4.0)]) / len(energies)}")
print(f"Fraction in 4.0 --> 4.5: {len(energies[(energies >= 4.0) & (energies < 4.5)]) / len(energies)}")
print(f"Fraction in 4.5 --> 5.0: {len(energies[(energies >= 4.5) & (energies <= 5.0)]) / len(energies)}")

print(f"\nFraction in 2.5 --> 3.75: {len(energies[(energies >= 2.5) & (energies <= 3.75)]) / len(energies)}")
print(f"Fraction in 3.75 --> 5.0: {len(energies[(energies >= 3.75) & (energies <= 5.0)]) / len(energies)}")
print(f"Fraction in 3.5 --> 5.0: {len(energies[(energies >= 3.5) & (energies <= 5.0)]) / len(energies)}")