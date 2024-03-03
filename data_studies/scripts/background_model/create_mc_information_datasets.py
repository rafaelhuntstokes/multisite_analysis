import numpy as np
import glob

"""
Script takes the run-by-run arrays of ITR, energy, position and nhits obtained
from Tl-208 simulations and saved individually and adds them into single 
distributions of each variable.
"""

def adder(files):
    """
    For a given list of files, open them one by one and add together.
    Save the final combined array.
    """

    combined = []
    for ifile in files:
        contents = np.load(ifile).tolist()

        combined += contents
    
    return combined

# gather run by run files
isotope = "BiPo212"
path = f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information"
energy_files    = glob.glob(path + "/energy/*.npy")
itr_files       = glob.glob(path + "/itr/*.npy")
nhits_files     = glob.glob(path + "/scaled_nhits/*.npy")
pos_x_files     = glob.glob(path + "/position/posX*.npy")
pos_y_files     = glob.glob(path + "/position/posY*.npy")
pos_z_files     = glob.glob(path + "/position/posZ*.npy")
pos_r_files     = glob.glob(path + "/position/posR*.npy")
rho2_files      = glob.glob(path + "/position/rho2*.npy")
residuals_files = glob.glob(path + "/time_residuals/*.npy")

# add all the files together
combined_energy    = adder(energy_files)
combined_itr       = adder(itr_files)
combined_nhits     = adder(nhits_files)
combined_x         = adder(pos_x_files)
combined_y         = adder(pos_y_files)
combined_z         = adder(pos_z_files)
combined_r         = adder(pos_r_files)
combined_rho2      = adder(rho2_files)
combined_residuals = adder(residuals_files)

# save combined arrays
np.save(path + "/energy.npy", combined_energy)
np.save(path + "/itr.npy", combined_itr)
np.save(path + "/nhits_scaled.npy", combined_nhits)
np.save(path + "/posX.npy", combined_x)
np.save(path + "/posY.npy", combined_y)
np.save(path + "/posZ.npy", combined_z)
np.save(path + "/posR.npy", combined_r)
np.save(path + "/rho2.npy", combined_rho2)
np.save(path + "/time_residuals.npy", combined_residuals)