import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import rat
import matplotlib.font_manager as fm
from ROOT import RAT
import ROOT
import json

path2 = '/data/snoplus3/hunt-stokes/nuPhysPoster/scripts/Times_New_Roman_Normal.ttf'
prop_font = fm.FontProperties(fname=path2, size = 12)

plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams.update({'axes.unicode_minus' : False})

# set position of axis labels
plt.rcParams["xaxis.labellocation"] = 'right'
plt.rcParams["yaxis.labellocation"] = 'top'

# set global parameters with rcParams -- copy whole block below before plotting code to set plotting template globally

# set height and width of big markings on axis x
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.6
# set height and width of small markings on axis x
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1.6
# set height and width of big markings on axis y
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 1.6
# set height and width of small markings on axis y
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 1.6
# set thickness of axes
plt.rcParams['axes.linewidth'] = 1.6
# set plot background color
plt.rcParams['figure.facecolor'] = 'white'
# set plot aspect ratio -- change according to needs
plt.rcParams['figure.figsize'] = (8.5, 6.5)
# set padding (between ticks and axis label)
plt.rcParams['xtick.major.pad'] = '12' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '12'
# set padding (between plot and title)
plt.rcParams['axes.titlepad'] = 12
# set markings on axis to show on the inside of the plot, can change if needed
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# set ticks on both sides of x and y axes
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

histogram_style = {
    'histtype': 'step', 
    'color': 'blue',
    'alpha': 0.7,
    'linewidth': 2
}

scatter_style = {
    'marker': 's',
    'color': 'black',
    's': 25
}

errorbar_style = {
    'linestyle': 'None',
    'color': 'black',
    'capsize': 1.5
}

line_plot_style = {
    'linestyle': '-',
    'color': 'black',
    'linewidth': 2.5
}

"""
Script creates a comparison plot showing the agreement between the energy spectrum in data and the 
background model expected rates.

The MC spectrum gives the SHAPE, and the calculated number of events scales the histograms such that
the area = number of expected events.

The comparison is done within the ROI : 2 - 6 MeV.

Two versions of the plot are created: one before removing tagged coincidences, and one after.

Each version is made for a range of different FV cuts.

Event Types Considered:

    1. BiPo214
    2. BiPo212
    3. alpha-n
    4. solar-nue
    5. solar-numu
    6. Tl-210
    7. Tl208
"""

def apply_scaling(energy_spectrum, expected_number, binning):
    """
    Function takes in the unbinned energy spectrum from MC, bins and scales it such that the sum of
    the counts is equal to the expected number of events given by the background model calculations.

    The poisson error on these scaled counts is calculated.

    INPUTS : raw energy spectrum of MC background
             expected number of events given by background model
             energy spectrum binning to be applied

    RETURNS: scaled counts, poisson error on counts, midpoint of bins.
    """

    # bin the raw energy spectrum
    unscaled_counts, bins = np.histogram(energy_spectrum, bins = binning)

    # normalise
    unscaled_counts = unscaled_counts / sum(unscaled_counts)

    # scale the counts such that the sum of counts equals expected number of events
    scaled_counts = unscaled_counts  * expected_number

    # find the midpoint of the bins for plotting
    mid_points = binning + np.diff(binning)[0] / 2

    # calculate the poisson error on the counts of each bin # 
    ################## TO BE IMPLEMENTED ####################

    return unscaled_counts, scaled_counts, mid_points[:-1]

def create_background_model(FV, coincidence_cuts, itr_cuts, dist_type, bin_width = 0.05):
    """
    For a given FV, the function creates the background model vs data energy spectrum for the Tl208 vs B8 solar ROI.

    INPUTS: chosen FV [3, 3.5, 4, 4.5, 5, 5.5, 6], float
            apply coincidence_cuts [True, False], bool
            apply itr_cuts [True, False], bool           | note that ITR cuts were optimised in another analysis to accept events if: 0.18 < ITR < 0.3
            dist_type "energy", "nhits_clean" or "nhits_corrected", str
            bin_width, float, the bin width for energy spectra
    OUTPUT: None, but saves a background model graph.
    """
    ITR_LOW = 0.18
    ITR_HIGH = 0.3
    print("hi")
    # load the MC energy and radial distributions
    path = "/data/snoplus3/hunt-stokes/clean_multisite"
    energy_bi214     = np.load(path + "/BiPo214_mc_information/energy.npy")
    energy_po214     = np.load(path + "/Po214_mc_information/energy.npy")
    energy_bi212     = np.load(path + "/BiPo212_mc_information/energy.npy")
    energy_po212     = np.load(path + "/Po212_mc_information/energy.npy")
    energy_alphaN_p  = np.load(path + "/AlphaN_LAB_13C_mc_information/energy.npy") # prompt spectrum + pileup
    energy_alphaN_d  = np.load(path + "/gammas_2p2MeV_mc_information/energy.npy")  # delayed 2.2 MeV gamma spectrum
    energy_tl208     = np.load(path + "/Tl208_mc_information/energy.npy")
    energy_tl210     = np.load(path + "/Tl210_mc_information/energy.npy")
    energy_B8_nue    = np.load(path + "/B8_solar_nue_mc_information/energy.npy")
    energy_B8_numu   = np.load(path + "/B8_solar_numu_mc_information/energy.npy")
    energy_pa234m    = np.load(path + "/Pa234m_mc_information/energy.npy")

    radius_bi214     = np.load(path + "/BiPo214_mc_information/posR.npy")
    radius_po214     = np.load(path + "/Po214_mc_information/posR.npy")
    radius_bi212     = np.load(path + "/BiPo212_mc_information/posR.npy")
    radius_po212     = np.load(path + "/Po212_mc_information/posR.npy")
    radius_alphaN_p  = np.load(path + "/AlphaN_LAB_13C_mc_information/posR.npy") 
    radius_alphaN_d  = np.load(path + "/gammas_2p2MeV_mc_information/posR.npy")
    radius_tl208     = np.load(path + "/Tl208_mc_information/posR.npy")
    radius_tl210     = np.load(path + "/Tl210_mc_information/posR.npy")
    radius_B8_nue    = np.load(path + "/B8_solar_nue_mc_information/posR.npy")
    radius_B8_numu   = np.load(path + "/B8_solar_numu_mc_information/posR.npy")
    radius_pa234m    = np.load(path + "/Pa234m_mc_information/posR.npy")

    itr_bi214     = np.load(path + "/BiPo214_mc_information/itr.npy")
    itr_po214     = np.load(path + "/Po214_mc_information/itr.npy")
    itr_bi212     = np.load(path + "/BiPo212_mc_information/itr.npy")
    itr_po212     = np.load(path + "/Po212_mc_information/itr.npy")
    itr_alphaN_p  = np.load(path + "/AlphaN_LAB_13C_mc_information/itr.npy") 
    itr_alphaN_d  = np.load(path + "/gammas_2p2MeV_mc_information/itr.npy")
    itr_tl208     = np.load(path + "/Tl208_mc_information/itr.npy")
    itr_tl210     = np.load(path + "/Tl210_mc_information/itr.npy")
    itr_B8_nue    = np.load(path + "/B8_solar_nue_mc_information/itr.npy")
    itr_B8_numu   = np.load(path + "/B8_solar_numu_mc_information/itr.npy")
    itr_pa234m    = np.load(path + "/Pa234m_mc_information/itr.npy")

    nhitscorrected_bi214     = np.load(path + "/BiPo214_mc_information/nhits_scaled.npy")
    nhitscorrected_po214     = np.load(path + "/Po214_mc_information/nhits_scaled.npy")
    nhitscorrected_bi212     = np.load(path + "/BiPo212_mc_information/nhits_scaled.npy")
    nhitscorrected_po212     = np.load(path + "/Po212_mc_information/nhits_scaled.npy")
    nhitscorrected_alphaN_p  = np.load(path + "/AlphaN_LAB_13C_mc_information/nhits_scaled.npy")
    nhitscorrected_alphaN_d  = np.load(path + "/gammas_2p2MeV_mc_information/nhits_scaled.npy")
    nhitscorrected_tl208     = np.load(path + "/Tl208_mc_information/nhits_scaled.npy")
    nhitscorrected_tl210     = np.load(path + "/Tl210_mc_information/nhits_scaled.npy")
    nhitscorrected_B8_nue    = np.load(path + "/B8_solar_nue_mc_information/nhits_scaled.npy")
    nhitscorrected_B8_numu   = np.load(path + "/B8_solar_numu_mc_information/nhits_scaled.npy")
    nhitscorrected_pa234m    = np.load(path + "/Pa234m_mc_information/nhits_scaled.npy")

    nhitsclean_bi214     = np.load(path + "/BiPo214_mc_information/nhits_cleaned.npy")
    nhitsclean_bi212     = np.load(path + "/BiPo212_mc_information/nhits_cleaned.npy")
    nhitsclean_alphaN_p  = np.load(path + "/AlphaN_LAB_13C_mc_information/nhits_cleaned.npy")
    nhitsclean_alphaN_d  = np.load(path + "/gammas_2p2MeV_mc_information/nhits_cleaned.npy")
    nhitsclean_tl208     = np.load(path + "/Tl208_mc_information/nhits_cleaned.npy")
    nhitsclean_tl210     = np.load(path + "/Tl210_mc_information/nhits_cleaned.npy")
    nhitsclean_B8_nue    = np.load(path + "/B8_solar_nue_mc_information/nhits_cleaned.npy")
    nhitsclean_B8_numu   = np.load(path + "/B8_solar_numu_mc_information/nhits_cleaned.npy")
    nhitsclean_pa234m    = np.load(path + "/Pa234m_mc_information/nhits_cleaned.npy")

    # plt.hist(nhitsclean_bi214, bins = np.arange(0, 1000, 1), histtype = "step")
    # isolate MC events which fall within the FV cut
    if itr_cuts == True:
        idx_bi214     = np.where((radius_bi214 < FV) & (itr_bi214 > ITR_LOW) & (itr_bi214 < ITR_HIGH) & (energy_bi214 > 2) & (energy_bi214 < 6))
        idx_po214     = np.where((radius_po214 < FV) & (itr_po214 > ITR_LOW) & (itr_po214 < ITR_HIGH) & (energy_po214 > 2) & (energy_po214 < 6))
        idx_bi212     = np.where((radius_bi212 < FV) & (itr_bi212 > ITR_LOW) & (itr_bi212 < ITR_HIGH) & (energy_bi212 > 2) & (energy_bi212 < 6))
        idx_po212     = np.where((radius_po212 < FV) & (itr_po212 > ITR_LOW) & (itr_po212 < ITR_HIGH) & (energy_po212 > 2) & (energy_po212 < 6))
        idx_alphaN_p  = np.where((radius_alphaN_p < FV) & (itr_alphaN_p > ITR_LOW) & (itr_alphaN_p < ITR_HIGH) & (energy_alphaN_p > 2) & (energy_alphaN_p < 6))
        idx_alphaN_d  = np.where((radius_alphaN_d < FV) & (itr_alphaN_d > ITR_LOW) & (itr_alphaN_d < ITR_HIGH) & (energy_alphaN_d > 2) & (energy_alphaN_d < 6))
        idx_tl208     = np.where((radius_tl208 < FV) & (itr_tl208 > ITR_LOW) & (itr_tl208 < ITR_HIGH) & (energy_tl208 > 2) & (energy_tl208 < 6))
        idx_tl210     = np.where((radius_tl210 < FV) & (itr_tl210 > ITR_LOW) & (itr_tl210 < ITR_HIGH) & (energy_tl210 > 2) & (energy_tl210 < 6))
        idx_B8_nue    = np.where((radius_B8_nue < FV) & (itr_B8_nue > ITR_LOW) & (itr_B8_nue < ITR_HIGH) & (energy_B8_nue > 2) & (energy_B8_nue < 6))
        idx_B8_numu   = np.where((radius_B8_numu < FV) & (itr_B8_numu > ITR_LOW) & (itr_B8_numu < ITR_HIGH) & (energy_B8_numu > 2) & (energy_B8_numu < 6))
        idx_pa234m    = np.where((radius_pa234m < FV) & (itr_pa234m > ITR_LOW) & (itr_pa234m < ITR_HIGH) & (energy_pa234m > 2) & (energy_pa234m < 6))
    else:
        idx_bi214     = np.where((radius_bi214 < FV)  & (energy_bi214 > 2.5) & (energy_bi214 < 5.0))
        idx_po214     = np.where((radius_po214 < FV) & (energy_po214 > 2.5) & (energy_po214 < 5.0))
        idx_bi212     = np.where((radius_bi212 < FV) &  (energy_bi212 > 2.5) & (energy_bi212 < 5.0))
        idx_po212     = np.where((radius_po212 < FV) &  (energy_po212 > 2.5) & (energy_po212 < 5.0))
        idx_alphaN_p  = np.where((radius_alphaN_p < FV) &  (energy_alphaN_p > 2.5) & (energy_alphaN_p < 5.0))
        idx_alphaN_d  = np.where((radius_alphaN_d < FV) &  (energy_alphaN_d > 2.5) & (energy_alphaN_d < 5.0))
        idx_tl208     = np.where((radius_tl208 < FV) & (energy_tl208 > 2.5) & (energy_tl208 < 5.0))
        idx_tl210     = np.where((radius_tl210 < FV) &  (energy_tl210 > 2.5) & (energy_tl210 < 5.0))
        idx_B8_nue    = np.where((radius_B8_nue < FV) &  (energy_B8_nue > 2.5) & (energy_B8_nue < 5.0))
        idx_B8_numu   = np.where((radius_B8_numu < FV) &  (energy_B8_numu > 2.5) & (energy_B8_numu < 5.0))
        idx_pa234m    = np.where((radius_pa234m < FV) &  (energy_pa234m > 2.5) & (energy_pa234m < 5.0))

    # apply idx
    if dist_type == "energy":
        binning            = np.arange(2.5, 5.0 + bin_width, bin_width)
        energy_bi214     = energy_bi214[idx_bi214]
        energy_po214     = energy_po214[idx_po214]
        energy_bi212     = energy_bi212[idx_bi212]
        energy_po212     = energy_po212[idx_po212]
        energy_alphaN_p  = energy_alphaN_p[idx_alphaN_p]
        energy_alphaN_d  = energy_alphaN_d[idx_alphaN_d]
        energy_tl208     = energy_tl208[idx_tl208]
        energy_tl210     = energy_tl210[idx_tl210]
        energy_B8_nue    = energy_B8_nue[idx_B8_nue]
        energy_B8_numu   = energy_B8_numu[idx_B8_numu]
        energy_pa234m    = energy_pa234m[idx_pa234m]
    elif dist_type == "nhits_cleaned":
        
        binning = np.arange(600, 1500, 20)
        energy_bi214     = nhitsclean_bi214[idx_bi214]
        energy_bi212     = nhitsclean_bi212[idx_bi212]
        energy_alphaN_p  = nhitsclean_alphaN_p[idx_alphaN_p]
        energy_alphaN_d  = nhitsclean_alphaN_d[idx_alphaN_d]
        energy_tl208     = nhitsclean_tl208[idx_tl208]
        energy_tl210     = nhitsclean_tl210[idx_tl210]
        energy_B8_nue    = nhitsclean_B8_nue[idx_B8_nue]
        energy_B8_numu   = nhitsclean_B8_numu[idx_B8_numu]
        energy_pa234m    = nhitsclean_pa234m[idx_pa234m]
    elif dist_type == "nhits_corrected":
        binning = np.arange(600, 1500, 20)
        energy_bi214     = nhitscorrected_bi214[idx_bi214]
        energy_po214     = nhitscorrected_po214[idx_po214]
        energy_bi212     = nhitscorrected_bi212[idx_bi212]
        energy_po212     = nhitscorrected_po212[idx_po212]
        energy_alphaN_p  = nhitscorrected_alphaN_p[idx_alphaN_p]
        energy_alphaN_d  = nhitscorrected_alphaN_d[idx_alphaN_d]
        energy_tl208     = nhitscorrected_tl208[idx_tl208]
        energy_tl210     = nhitscorrected_tl210[idx_tl210]
        energy_B8_nue    = nhitscorrected_B8_nue[idx_B8_nue]
        energy_B8_numu   = nhitscorrected_B8_numu[idx_B8_numu]
        energy_pa234m    = nhitscorrected_pa234m[idx_pa234m]
    else:
        print("Dist type not recognised.")
        return 0
    
    # map the inputted FV floats to the names (strings) of the data TTrees
    name_map  = {3000.0: "3m", 3500.0: "3p5m", 4000.0: "4m", 4500.0: "4p5m", 5000.0: "5m", 5500.0: "5p5m", 6000.0: "6m"}
    name_tree = name_map[float(FV)]
    energy_data = np.load(f"./data_energy_spectrums/{dist_type}_dist_directionality_quiet_4500_newDAG.npy")

    # load JSON file containing dictionary of isotopes and their calculated expected rates
    expected_backgrounds_rates = open("expected_rates_database.json")
    background_rates_db        = json.load(expected_backgrounds_rates)

    # find the expected number of events for each background from the database
    if coincidence_cuts == True:
        # expected_bipo214   = background_rates_db["post_cuts"][name_tree]["BiPo214"] 
        # expected_bipo212   = background_rates_db["post_cuts"][name_tree]["BiPo212"]
        # expected_tl208     = background_rates_db["post_cuts"][name_tree]["Tl208"]
        # expected_tl210     = background_rates_db["post_cuts"][name_tree]["Tl210"]
        # expected_B8_nue    = background_rates_db["post_cuts"][name_tree]["B8_nue"]
        # expected_B8_numu   = background_rates_db["post_cuts"][name_tree]["B8_numu"]
        # expected_alphaN_p  = background_rates_db["post_cuts"][name_tree]["alphaN_prompt"]
        # expected_alphaN_d  = background_rates_db["post_cuts"][name_tree]["alphaN_delayed"]
        # expected_pa234m    = background_rates_db["post_cuts"][name_tree]["Pa234m"]

        expected_bipo214   = background_rates_db["quiet_period"][name_tree]["BiPo214"] 
        expected_bipo212   = background_rates_db["quiet_period"][name_tree]["BiPo212"]
        expected_tl208     = background_rates_db["quiet_period"][name_tree]["Tl208"]
        expected_tl210     = background_rates_db["quiet_period"][name_tree]["Tl210"]
        expected_B8_nue    = background_rates_db["quiet_period"][name_tree]["B8_nue"]
        expected_B8_numu   = background_rates_db["quiet_period"][name_tree]["B8_numu"]
        expected_alphaN_p  = background_rates_db["quiet_period"][name_tree]["alphaN_prompt"]
        expected_alphaN_d  = background_rates_db["quiet_period"][name_tree]["alphaN_delayed"]
        expected_pa234m    = background_rates_db["quiet_period"][name_tree]["Pa234m"]
    if coincidence_cuts == False:
        # expected_bipo214   = background_rates_db["pre_cuts"][name_tree]["BiPo214"] 
        # expected_bipo212   = background_rates_db["pre_cuts"][name_tree]["BiPo212"]
        # expected_tl208     = background_rates_db["pre_cuts"][name_tree]["Tl208"]
        # expected_tl210     = background_rates_db["pre_cuts"][name_tree]["Tl210"]
        # expected_B8_nue    = background_rates_db["pre_cuts"][name_tree]["B8_nue"]
        # expected_B8_numu   = background_rates_db["pre_cuts"][name_tree]["B8_numu"]
        # expected_alphaN_p  = background_rates_db["pre_cuts"][name_tree]["alphaN_prompt"]
        # expected_alphaN_d  = background_rates_db["pre_cuts"][name_tree]["alphaN_delayed"]
        expected_bipo214   = background_rates_db["quiet_period"][name_tree]["BiPo214"] 
        expected_bipo212   = background_rates_db["quiet_period"][name_tree]["BiPo212"]
        expected_tl208     = background_rates_db["quiet_period"][name_tree]["Tl208"]
        expected_tl210     = background_rates_db["quiet_period"][name_tree]["Tl210"]
        expected_B8_nue    = background_rates_db["quiet_period"][name_tree]["B8_nue"]
        expected_B8_numu   = background_rates_db["quiet_period"][name_tree]["B8_numu"]
        expected_alphaN_p  = background_rates_db["quiet_period"][name_tree]["alphaN_prompt"]
        expected_alphaN_d  = background_rates_db["quiet_period"][name_tree]["alphaN_delayed"]
        expected_pa234m    = background_rates_db["quiet_period"][name_tree]["Pa234m"]
    
    # scale the MC to the expected number of events
    print(binning)
    print("Scaled 1")
    pdf_bi214, scaled_bi214, mids = apply_scaling(energy_bi214, expected_bipo214, binning)
    print("Scaled 2")
    pdf_bi212, scaled_bi212, _    = apply_scaling(energy_bi212, expected_bipo212, binning)
    print("Scaled 3")
    pdf_tl208, scaled_tl208, _    = apply_scaling(energy_tl208, expected_tl208, binning)
    print("Scaled 4")
    pdf_tl210, scaled_tl210, _    = apply_scaling(energy_tl210, expected_tl210, binning)
    print("Scaled 5")
    pdf_alphaN_p, scaled_alphaN_p, _ = apply_scaling(energy_alphaN_p, expected_alphaN_p, binning)
    print("Scaled 6")
    pdf_alphaN_d, scaled_alphaN_d, _ = apply_scaling(energy_alphaN_d, expected_alphaN_d, binning)
    print("Scaled 7")
    pdf_B8_nue, scaled_B8_nue, _   = apply_scaling(energy_B8_nue, expected_B8_nue, binning)
    print("Scaled 8")
    pdf_B8_numu, scaled_B8_numu, _  = apply_scaling(energy_B8_numu, expected_B8_numu, binning)
    print("Scaled 9")
    # scaled_pa234m, _   = apply_scaling(energy_pa234m, expected_pa234m, binning)
    print("Scaled 10")

    print("Number Scaled BiPo212: ", np.sum(scaled_bi212))
    print("Number Scaled BiPo214: ", np.sum(scaled_bi214))
    # bin the data and find poisson error on the counts
    counts_data, _ = np.histogram(energy_data, bins = binning)
    error_data = np.sqrt(counts_data)
    print("last bin counts", counts_data[-1])
    # create the plot!
    fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (12, 12))
    # plt.figure(figsize = (12, 6)) 
    axes[0].errorbar(mids, counts_data, yerr = error_data, capsize = 2, marker = "o", markersize = 5, linestyle = "", color = "black", label = f"Data | Total Counts: {round(sum(counts_data), 2)}")
    axes[0].errorbar(mids, scaled_bi214, capsize = 2, marker = "x", markersize = 2, label = rf"$^{{214}}$Bi | {round(sum(scaled_bi214), 2)}")
    axes[0].errorbar(mids, scaled_bi212, capsize = 2, marker = "x", markersize = 2, label = rf"$^{{212}}$Bi | {round(sum(scaled_bi212), 2)}")
    axes[0].errorbar(mids, scaled_tl208, capsize = 2, marker = "x", markersize = 2, label = rf"$^{{208}}$Tl | {round(sum(scaled_tl208), 2)}")
    axes[0].errorbar(mids, scaled_tl210, capsize = 2, marker = "x", markersize = 2, label = rf"$^{{210}}$Tl | {round(sum(scaled_tl210), 2)}")
    # axes[0].errorbar(mids, scaled_alphaN_p, capsize = 2, marker = "x", markersize = 2, label = r"$(\alpha, n)$" + "| Prompt")
    # axes[0].errorbar(mids, scaled_alphaN_d, capsize = 2, marker = "x", markersize = 2, label = r"$(\alpha, n)$" + "| Delayed")
    axes[0].errorbar(mids, scaled_B8_nue, capsize = 2, marker = "x", markersize = 2, label = f"$^8$B " + rf"$\nu _e$ | {round(sum(scaled_B8_nue), 2)}")
    axes[0].errorbar(mids, scaled_B8_numu, capsize = 2, marker = "x", markersize = 2, label = f"$^8$B " + rf"$\nu _\mu$ | {round(sum(scaled_B8_numu), 2)}")
    # axes[0].errorbar(mids, scaled_pa234m, capsize = 2, marker = "x", markersize = 2, label = "Pa234m")
    
    total_model    = scaled_bi214 +  scaled_bi212  + scaled_tl208 + scaled_tl210 + scaled_B8_nue + scaled_B8_numu #+ scaled_alphaN_d + scaled_alphaN_p #+ scaled_pa234m
    
    # find absolute error on each of the predicted model counts in whole ROI
    stat_err_norms = np.array([0.569116091, 5.29325162, 42.16642339, 0.01375316126])  # ignoring error on the B8
    
    # find the error in each bin by multiplying each absolute error in the norms by the normalised PDF counts for each bin
    err_bi214 = stat_err_norms[0] * pdf_bi214
    err_bi212 = stat_err_norms[1] * pdf_bi212
    err_tl208 = stat_err_norms[2] * pdf_tl208
    err_tl210 = stat_err_norms[3] * pdf_tl210
    print(err_bi214.shape, err_bi212.shape, err_tl208.shape, err_tl210.shape)

    # quadrature sum the errors
    error_model    = np.sqrt(err_bi214 **2 + err_bi212 **2 + err_tl208 **2 + err_tl210 **2)
    print("Total model last bin: ", total_model[-1])
    print("Model error las bin: ", error_model[-1])
    # WORK OUT THE Chi2 / NDF for the model!
    counts_data_chi2  = counts_data 
    counts_model_chi2 = total_model
    bins_chi2         = mids
    error_model       = np.sqrt(total_model)
    print(error_model)
    print(len(error_model), len(counts_data), len(counts_model_chi2), len(mids))
    chi2 = np.sum(((counts_data_chi2 - counts_model_chi2)**2/counts_model_chi2))
    print(chi2)
    NDF  = len(bins_chi2) - 6 # num PDFs (ignoring Pa234m as it's outside the ROI zone)
    axes[0].plot([], [], linestyle = "", label = r"$\frac{\chi^2}{NDF} = $ " + f"{round(chi2/NDF, 2)}")
    axes[0].errorbar(mids, total_model, capsize = 2, marker = "x", markersize = 2, color = "red", label = f"Model: Total Counts {round(sum(total_model), 2)}")
    axes[0].legend(fontsize = 10, frameon = False)
    # axes[0].set_title(f"Background Model | FV: {name_tree} | Coincidence Cuts: {coincidence_cuts} | ITR Cut: {itr_cuts}", fontsize = 20)
    
    axes[0].set_ylabel(f"Counts per {bin_width} MeV Bin", fontsize = 20, fontproperties = prop_font)
    # axes[0].set_yscale("log")
    if dist_type == "energy":
        axes[0].set_xlabel("Reconstructed Energy (MeV)", fontsize = 20,fontproperties = prop_font)
        axes[0].set_xlim((2.5,5.0))
        axes[0].set_ylim((0, 60))
    if dist_type == "nhits_corrected":
        axes[0].set_xlabel("Nhits Corrected", fontsize = 20, fontproperties = prop_font)
        axes[0].set_xlim((700,1500))
        axes[0].set_ylim((0, 100))
    if dist_type == "nhits_cleaned":
        axes[0].set_xlabel("Nhits Clean", fontsize = 20, fontproperties = prop_font)
        axes[0].set_xlim((600,1500))
        axes[0].set_ylim((0, 100))
    
    for label in axes[0].get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in axes[0].get_yticklabels():
        label.set_fontproperties(prop_font)
    idx = np.where(error_data != 0)
    # draw the Data - Model Residual
    check = (1/NDF) * np.sum(((counts_data - total_model)**2 / error_model**2))
    print(check)
    # axes[1].bar(mids,  ((counts_data - total_model)**2 / error_model**2), width = 0.05, color = "none", edgecolor = "black")
    axes[1].bar(mids,  (counts_data - total_model) / error_data, width = 0.05, color = "none", edgecolor = "black")
    axes[1].axhline(y = 0, color = "black")
    # axes[1].title(f"Background Model Residual | FV: {name_tree} | Coincidence Cuts: {coincidence_cuts} | ITR Cut: {itr_cuts}")
    axes[1].set_ylabel(r"Residuals $\left(\frac{\text{Data} - \text{Model}}{\sigma_{D}^2}\right)$", fontsize = 20, fontproperties = prop_font)
    
    if dist_type == "energy":
        axes[1].set_xlabel("Reconstructed Energy (MeV)", fontsize = 20, fontproperties = prop_font)
        axes[1].set_xlim((2.5,5.0))
        # axes[1].set_ylim((0, 3))
    if dist_type == "nhits_corrected":
        axes[1].set_xlabel("Nhits Corrected", fontsize = 20, fontproperties = prop_font)
        axes[1].set_xlim((700,1500))
        axes[1].set_ylim((-4, 4))
    if dist_type == "nhits_cleaned":
        axes[1].set_xlabel("Nhits Clean", fontsize = 20, fontproperties = prop_font)
        axes[1].set_xlim((700,1500))
        axes[1].set_ylim((-3, 3))
    for label in axes[1].get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in axes[1].get_yticklabels():
        label.set_fontproperties(prop_font)
    # axes[1].axhline(y = 0, color = "red")
    axes[1].set_ylim((-5, 5))
    # plt.show()
    # plt.savefig(f"residual_{name_tree}_quiet.png")
    plt.savefig(f"../../plots/background_models/{dist_type}_{name_tree}_quiet_68runsAdded_newDAG_2.pdf")
    # plt.savefig("./BI214_NHITS_CLEANED.pdf")
def create_data_spectrum_with_runlist():
    """
    Create a spectrum from the data by opening the individual run by run files.
    """

    # set up
    livetime_path = "/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence" 
    data_path     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/processed_dataset"
    run_list      = np.loadtxt("/data/snoplus3/hunt-stokes/multisite_clean/data_studies/runlists/quiet_period.txt", dtype = int)
    counter       = 0 

    coincidence = True # apply coincidence cuts to data
    fv_cut      = 4500 # fv cut
    itr_low     = 0.18 # itr cuts
    itr_high    = 0.3
    
    # track missing / negative livetime runs in the run list
    bad_runs     = []
    missing_runs = []

    # outputs to save
    dist_energy          = []
    dist_nhits_clean     = []
    dist_nhits_corrected = []
    total_livetime       = 0
    
    for irun in run_list:

        # try and open the corresponding root file
        try:
            file     = ROOT.TFile.Open(f"{data_path}/{irun}.root")
            livetime = np.loadtxt(f"{livetime_path}/livetime/{irun}.txt")
            if livetime < 0:
               bad_runs.append(irun)
               print(irun, livetime)
               continue
            if coincidence == True:
                # ntuple = file.Get("clean_4p5m")
                ntuple = file.Get("2p5_5p0")
                print(ntuple.GetEntries())
            if coincidence == False:
                ntuple = file.Get("full_4p5m")
                print(ntuple.GetEntries())
        except:
            missing_runs.append(irun)
            continue
        
        total_livetime += livetime

        for ientry in ntuple:
            
            # apply ITR cut
            # itr = ientry.ITR
            # if itr < itr_low or itr > itr_high:
            #     continue

            energy          = ientry.energy

            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > fv_cut:
                continue
            # nhits_corrected = ientry.nhitsCorrected
            nhits_clean     = ientry.nhitsCleaned
            dist_energy.append(energy)
            dist_nhits_clean.append(nhits_clean)
            # dist_nhits_corrected.append(nhits_corrected)

        counter += 1 
        print(f"Completed {counter} / {len(run_list)} ({(100 * counter) / len(run_list)} %)")
    
    print("Num events: ", len(dist_energy))
    print(f"\nTotal Livetime: {total_livetime / (60 * 60)} hrs | {total_livetime / (60 * 60 * 24)} days.")
    for i in missing_runs:
        print(i)
    print(f"Missing {len(missing_runs)} runs.")
    print(f"Found {len(bad_runs)} bad runs:")
    print(bad_runs)

    
    # np.save(f"./data_energy_spectrums/energy_dist_directionality_quiet_{fv_cut}.npy", dist_energy)
    # np.save(f"./data_energy_spectrums/nhits_corrected_dist_directionality_quiet_{fv_cut}.npy", dist_nhits_corrected)
    # np.save(f"./data_energy_spectrums/nhits_cleaned_dist_directionality_quiet_{fv_cut}.npy", dist_nhits_clean)

    np.save(f"./data_energy_spectrums/energy_dist_directionality_quiet_{fv_cut}_newDAG.npy", dist_energy)
    # np.save(f"./data_energy_spectrums/nhits_corrected_dist_directionality_quiet_{fv_cut}_newDAG.npy", dist_nhits_corrected)
    np.save(f"./data_energy_spectrums/nhits_cleaned_dist_directionality_quiet_{fv_cut}_newDAG.npy", dist_nhits_clean)


def compare_energy_dists():

    nhits = np.load("./nhit_cleaned_dist_gold_4500_counts.npy")
    binning = np.arange(0, 1500, 10)
    counts, edges = np.histogram(nhits, bins = binning)
    mids = edges[:-1] + np.diff(edges)[0]/2
    error = np.sqrt(counts)

    plt.figure(figsize = (12,6))
    plt.errorbar(mids, counts, yerr = error, capsize=2, linestyle = "", color = "black", marker = "o")
    plt.legend()
    plt.xlabel("Nhits Cleaned")
    plt.ylabel("Counts per ")
    plt.xlim((600, 1500))
    plt.ylim((0, 100))
    
    plt.savefig("nhit_dist.pdf")

def count_bump_events():
    fvs         = [3000, 3500, 4000, 4500, 5000, 5500]
    num_bump_fv = []
    bump_low    = 4.2
    bump_high   = 4.4

    for i in range(len(fvs)):
        # load the energy spectrum associated with this fv
        energy = np.load(f"energy_dist_directionality_{fvs[i]}_counts.npy")

        # count how many fall into the bump window
        bump_events = energy[(energy >= bump_low) & (energy <= bump_high)]

        num_bump_fv.append(len(bump_events))

    print(num_bump_fv)

def check_uniformity():
    """
    Function checks the uniformity of data events in the 2.5 --> 3.0 MeV range.
    This is to check whether the excess BiPo214 events I'm seeing in the data come
    from the internal ropes. An antisymmetric distribution in 4.5 m is expected,
    but no asymmetry in 3.3 m FV.
    
    This will also check the asymmetry (or otherwise) of the Tl208 events (3.5 -> 4.0 MeV).
    """

    # set up
    data_path    = "/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence" 
    run_list     = np.loadtxt("../../runlists/directionality_list.txt", dtype = int)
    counter      = 0 

    coincidence  = True                     # apply coincidence cuts to data
    fv_cuts      = [4500, 3300]             # the FV cuts to apply to the data
    energy_cuts  = [[2.5, 3.0], [3.5, 4.0]] # energy cuts to check the BiPo and the Tl208
    
    events_bipo_smallFV_x = []
    events_bipo_smallFV_y = []
    events_bipo_smallFV_z = []
    events_bipo_smallFV_p = []
    events_bipo_smallFV_itr = []
    events_bipo_bigFV_x   = []
    events_bipo_bigFV_y   = []
    events_bipo_bigFV_z   = []
    events_bipo_bigFV_p   = []
    events_bipo_bigFV_itr = []
    events_tl_smallFV_x   = []
    events_tl_smallFV_y   = []
    events_tl_smallFV_z   = []
    events_tl_smallFV_p   = []
    events_tl_smallFV_itr = []
    events_tl_bigFV_x     = []
    events_tl_bigFV_y     = []
    events_tl_bigFV_z     = []
    events_tl_bigFV_p     = []
    events_tl_bigFV_itr   = []

    # save the negative livetime fucked up runs
    bad_runs     = []
    missing_runs = []
    total_livetime = 0
    for irun in run_list:

        # try and open the corresponding root file
        try:
            file     = ROOT.TFile.Open(f"{data_path}/output_{irun}.root")
            livetime = np.loadtxt(f"{data_path}/livetime/{irun}.txt")
            if livetime < 0:
               bad_runs.append(irun)
               print(irun, livetime)
               continue
            if coincidence == True:
                ntuple = file.Get("clean_6m")
            if coincidence == False:
                ntuple = file.Get("full_6m")
        except:
            missing_runs.append(irun)
            continue
        
        total_livetime += livetime

        for ientry in ntuple:
            
            # apply ITR cut
            itr    = ientry.ITR
            energy = ientry.energy
            x      = ientry.x
            y      = ientry.y
            z      = ientry.z
            r      = np.sqrt(x**2 + y**2 + z**2)
            rho    = (x**2 + y**2) / (6000**2)
            if r < fv_cuts[0]:

                # we are inside the 4.5 m FV! 
                if energy >= energy_cuts[0][0] and energy <= energy_cuts[0][1]:
                    # we are inside the BiPo energy cuts
                    events_bipo_bigFV_x.append(x)
                    events_bipo_bigFV_y.append(y)
                    events_bipo_bigFV_z.append(z)
                    events_bipo_bigFV_p.append(rho)
                    events_bipo_bigFV_itr.append(itr)
                if energy >= energy_cuts[1][0] and energy <= energy_cuts[1][1]:
                    # we are inside the Tl energy cuts
                    events_tl_bigFV_x.append(x)
                    events_tl_bigFV_y.append(y)
                    events_tl_bigFV_z.append(z)
                    events_tl_bigFV_p.append(rho)
                    events_tl_bigFV_itr.append(itr)
                if r < fv_cuts[1]:
                    
                    # we are inside the 3.3 m FV!
                    if energy >= energy_cuts[0][0] and energy <= energy_cuts[0][1]:
                        # we are inside the BiPo energy cuts
                        events_bipo_smallFV_x.append(x)
                        events_bipo_smallFV_y.append(y)
                        events_bipo_smallFV_z.append(z)
                        events_bipo_smallFV_p.append(rho)
                        events_bipo_smallFV_itr.append(itr)
                    if energy >= energy_cuts[1][0] and energy <= energy_cuts[1][1]:
                        # we are inside the Tl energy cuts
                        events_tl_smallFV_x.append(x)
                        events_tl_smallFV_y.append(y)
                        events_tl_smallFV_z.append(z)
                        events_tl_smallFV_p.append(rho)
                        events_tl_smallFV_itr.append(itr)

        counter += 1 
        print(f"Completed {counter} / {len(run_list)} ({(100 * counter) / len(run_list)} %)")
    print(missing_runs)
    print(bad_runs)

    font_s = 20
    # create the 4 output plots (2 x 2)
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 15))

    axes[0,0].hist2d(events_bipo_bigFV_x, events_bipo_bigFV_y, bins = [np.arange(-6000, 6100, 200), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[0,0].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[0,0].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[0,0].set_title(r"2.5 $\rightarrow$ 3.0 MeV | 4.5 m FV | $N_{events}$ = " + f"{len(events_bipo_bigFV_p)}", fontsize = font_s)

    axes[0,1].hist2d(events_bipo_bigFV_p, events_bipo_bigFV_z, bins = [np.arange(0, 1.05, 0.05), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[0,1].set_xlabel(r"$\left(\frac{\rho}{\rho _{AV}}\right)^2$", fontsize = font_s)
    axes[0,1].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[0,1].set_title(r"2.5 $\rightarrow$ 3.0 MeV | 4.5 m FV | $N_{events}$ = " + f"{len(events_bipo_bigFV_p)}", fontsize = font_s)

    axes[1,0].hist2d(events_bipo_smallFV_x, events_bipo_smallFV_y, bins = [np.arange(-6000, 6200, 200), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[1,0].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[1,0].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[1,0].set_title(r"2.5 $\rightarrow$ 3.0 MeV | 3.3 m FV | $N_{events}$ = " + f"{len(events_bipo_smallFV_p)}", fontsize = font_s)

    axes[1,1].hist2d(events_bipo_smallFV_p, events_bipo_smallFV_z, bins = [np.arange(0, 1.05, 0.05), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[1,1].set_xlabel(r"$\left(\frac{\rho}{\rho _{AV}}\right)^2$", fontsize = font_s)
    axes[1,1].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[1,1].set_title(r"2.5 $\rightarrow$ 3.0 MeV | 3.3 m FV | $N_{events}$ = " + f"{len(events_bipo_smallFV_p)}", fontsize = font_s)

    fig.tight_layout()
    plt.savefig("../../plots/background_models/non_uniformity_test_bipo.png")
    plt.close()

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 15))

    axes[0,0].hist2d(events_tl_bigFV_x, events_tl_bigFV_y, bins = [np.arange(-6000, 6200, 200), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[0,0].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[0,0].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[0,0].set_title(r"3.5 $\rightarrow$ 4.0 MeV | 4.5 m FV | $N_{events}$ = " + f"{len(events_tl_bigFV_p)}", fontsize = font_s)
    
    axes[0,1].hist2d(events_tl_bigFV_p, events_tl_bigFV_z, bins = [np.arange(0, 1.05, 0.05), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[0,1].set_xlabel(r"$\left(\frac{\rho}{\rho _{AV}}\right)^2$", fontsize = font_s)
    axes[0,1].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[0,1].set_title(r"3.5 $\rightarrow$ 4.0 MeV | 4.5 m FV | $N_{events}$ = " + f"{len(events_tl_bigFV_p)}", fontsize = font_s)

    axes[1,0].hist2d(events_tl_smallFV_x, events_tl_smallFV_y, bins = [np.arange(-6000, 6200, 200), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[1,0].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[1,0].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[1,0].set_title(r"3.5 $\rightarrow$ 4.0 MeV | 3.3 m FV | $N_{events}$ = " + f"{len(events_tl_smallFV_p)}", fontsize = font_s)

    axes[1,1].hist2d(events_tl_smallFV_p, events_tl_smallFV_z, bins = [np.arange(0, 1.05, 0.05), np.arange(-6000, 6200, 200)], cmin = 1e-4)
    axes[1,1].set_xlabel(r"$\left(\frac{\rho}{\rho _{AV}}\right)^2$", fontsize = font_s)
    axes[1,1].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[1,1].set_title(r"3.5 $\rightarrow$ 4.0 MeV | 3.3 m FV | $N_{events}$ = " + f"{len(events_tl_smallFV_p)}", fontsize = font_s)

    fig.tight_layout()
    plt.savefig("../../plots/background_models/non_uniformity_test_tl208.png")
    plt.close()
    
def count_bipo_in_runlist():
    """
    Function counts the number of bipo events tagged in each run of the run list.
    Returns the total tagged in each FV.
    """

    isotope       = 212
    if isotope == 214:
        folder        = f"/data/snoplus3/hunt-stokes/clean_multisite/bipo_{isotope}_cleanSelec" 
    if isotope == 212:
        folder        = f"/data/snoplus3/hunt-stokes/clean_multisite/bipo_{isotope}_cleanSelec" 
    run_list      = np.loadtxt("../../runlists/quiet_period.txt", dtype = int)
    missing       = []           # keep track of missing runs
    number_per_fv = np.zeros(7)  # record number of bipo tagged in each FV
    fv_names      = ["6m", "5p5m", "5m", "4p5m", "4m", "3p5m", "3m"]
    for irun in run_list:
            
        # attempt to open the corresponding bipo tagged run file
        try:
            file   = ROOT.TFile.Open(f"{folder}/output_{isotope}_{irun}.root")
            ntuple = file.Get("6m")

        except:
            missing.append(irun)
            continue

        for ifv in range(len(fv_names)):
            ntuple = file.Get(fv_names[ifv])

            num_entries = ntuple.GetEntries()
            number_per_fv[ifv] += num_entries

    print(f"BiPo{isotope} tagged in each FV: {number_per_fv}")
    for i in missing:
        print(i)
    print(f"Missing {len(missing)} runs.")

print("hi")
# check_uniformity()
# create_data_spectrum_with_runlist()
create_background_model(4500, True, False, "energy", 0.05)
# create_background_model(4500, True, True, "nhits_corrected", 0.05)

# count_bipo_in_runlist()
