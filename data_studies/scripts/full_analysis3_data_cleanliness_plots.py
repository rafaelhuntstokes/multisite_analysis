import ROOT
import numpy as np
import matplotlib
import matplotlib.font_manager as fm
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import rat 
from ROOT import RAT
import glob

"""
USE env_poster in nuPhysPoster!
"""
path = '/data/snoplus3/hunt-stokes/nuPhysPoster/scripts/Times_New_Roman_Normal.ttf'
prop_font = fm.FontProperties(fname=path, size = 15)

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
Given removal of retriggers and orphans from data set, plot the data cleanliness variables
to identify any remaining bad events in the dataset.
"""

def load_event_info():
    # load every event in 6 m into numpy arrays
    # file = ROOT.TFile.Open("../extracted_data/full_analysis3/reprocessed_dataset_7.0.15/total.root")
    file = ROOT.TFile.Open("../extracted_data/full_analysis3/processed_dataset/total.root")
    tree = file.Get("2p5_5p0")

    x           = []
    y           = []
    z           = []
    r           = []
    energy      = []
    nhits_clean = []
    HC_ratio    = []
    posFOM      = []
    posFOM_hits = []
    itr         = []
    dlogL       = []

    for ientry in tree:
        x.append(ientry.x)
        y.append(ientry.y)
        z.append(ientry.z)
        r.append(np.sqrt(ientry.x**2 + ientry.y**2 + ientry.z**2))
        energy.append(ientry.energy)
        nhits_clean.append(ientry.nhitsCleaned)
        HC_ratio.append(ientry.HC_ratio)
        posFOM.append(ientry.posFOM)
        posFOM_hits.append(ientry.posFOM_hits)
        dlogL.append(ientry.dlogL)
        itr.append(ientry.itr)

    np.save("./full_analysis3_6m_x.npy", x)
    np.save("./full_analysis3_6m_y.npy", y)
    np.save("./full_analysis3_6m_z.npy", z)
    np.save("./full_analysis3_6m_r.npy", r)
    np.save("./full_analysis3_6m_energy.npy", energy)
    np.save("./full_analysis3_6m_nhits_clean.npy", nhits_clean)
    np.save("./full_analysis3_6m_HC_ratio.npy", HC_ratio)
    np.save("./full_analysis3_6m_posFOM.npy", posFOM)
    np.save("./full_analysis3_6m_posFOM_hits.npy", posFOM_hits)
    np.save("./full_analysis3_6m_itr.npy", itr)
    np.save("./full_analysis3_6m_dlogL.npy", dlogL)

    # np.save("./full_analysis3_6m_x_reproc.npy", x)
    # np.save("./full_analysis3_6m_y_reproc.npy", y)
    # np.save("./full_analysis3_6m_z_reproc.npy", z)
    # np.save("./full_analysis3_6m_r_reproc.npy", r)
    # np.save("./full_analysis3_6m_energy_reproc.npy", energy)
    # np.save("./full_analysis3_6m_nhits_clean_reproc.npy", nhits_clean)
    # np.save("./full_analysis3_6m_HC_ratio_reproc.npy", HC_ratio)
    # np.save("./full_analysis3_6m_posFOM_reproc.npy", posFOM)
    # np.save("./full_analysis3_6m_posFOM_hits_reproc.npy", posFOM_hits)
    # np.save("./full_analysis3_6m_itr_reproc.npy", itr)
    # np.save("./full_analysis3_6m_dlogL_reproc.npy", dlogL)
def create_1D_plots(FV_CUT):
    """
    Function loads the variables of interest and plots each quantity in 1 D, 
    for a given FV.
    """

    ### create the plots ##
    plt.figure()
    plt.hist(energy, bins = energy_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("Reconstructed Energy [MeV]", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/energy_{FV_CUT}_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(itr, bins = itr_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("ITR", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/itr_{FV_CUT}_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(hc_ratio, bins = hc_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("HC Ratio", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.35, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/hc_ratio_{FV_CUT}_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(posFOM/posFOM_hits, bins = posFOM_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("scintFitter posFOM/posFOM_hits", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/posFOM_{FV_CUT}_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(dlogL, bins = dlogL_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.15, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/dlogL_{FV_CUT}_reproc.pdf")
    plt.close()

def create_1D_plots_postCuts(FV_CUT):
    """
    Function loads the variables of interest and plots each quantity in 1 D, 
    for a given FV.
    """

    ### create the plots ##
    plt.figure()
    plt.hist(energy, bins = energy_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("Reconstructed Energy [MeV]", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/energy_{FV_CUT}_postCuts_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(itr, bins = itr_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("ITR", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/itr_{FV_CUT}_postCuts_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(hc_ratio, bins = hc_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("HC Ratio", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.35, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/hc_ratio_{FV_CUT}_postCuts_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(posFOM/posFOM_hits, bins = posFOM_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel("scintFitter posFOM/posFOM_hits", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.55, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/posFOM_{FV_CUT}_postCuts_reproc.pdf")
    plt.close()

    plt.figure()
    plt.hist(dlogL, bins = dlogL_bins, histtype="step", color = "black", linewidth = 2)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font)
    plt.ylabel("Counts", fontproperties = prop_font)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax.text(0.15, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/dlogL_{FV_CUT}_postCuts_reproc.pdf")
    plt.close()
def create_2D_plots(quant1, quant2, bins1, bins2, name1, name2):
    """
    Function creates a 2D histogram from quantity 1 and quantity2 with their
    defined binning.
    """

    plt.figure()
    plt.hist2d(quant1, quant2, bins = [bins1, bins2], cmin = 1e-4)
    plt.xlabel(f"{name1}", fontproperties = prop_font)
    plt.ylabel(f"{name2}", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    # ax.text(0.15, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    plt.title(f"FV {FV_CUT/1000:.1f} m")
    plt.savefig(f"{save_dir}/{name1}_{name2}_{FV_CUT}_reproc.pdf")
    plt.close()

def create_2D_plots_postCuts(quant1, quant2, bins1, bins2, name1, name2):
    """
    Function creates a 2D histogram from quantity 1 and quantity2 with their
    defined binning.
    """

    plt.figure()
    plt.hist2d(quant1, quant2, bins = [bins1, bins2], cmin = 1e-4)
    plt.xlabel(f"{name1}", fontproperties = prop_font)
    plt.ylabel(f"{name2}", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    # ax.text(0.15, 0.65, f"FV {FV_CUT/1000:.1f} m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    plt.title(f"FV {FV_CUT/1000:.1f} m")
    plt.savefig(f"{save_dir}/{name1}_{name2}_{FV_CUT}_postCuts_reproc.pdf")
    plt.close() 

def event_cuts(ITR_CUT, posFOM_CUT):
    """
    Function checks the events cut out by a posFOM and ITR cut. Are the same
    events cut out from both these cuts? How many events does the ITR and posFOM
    remove individually? How many combined?
    """

    # normalise posFOM by nhits
    posFOM_norm  = posFOM / posFOM_hits
    itr_id_rm    = []
    posFOM_id_rm = []
    both_id_rm   = []
    only_itr     = []
    only_posFOM  = []
    # loop over each event and find which cuts remove which events
    print(f"Total number of events: {itr.size}")
    for iev in range(itr.size):
        itr_flag    = False
        posFOM_flag = False
        both_flag   = False
        if itr[iev] < ITR_CUT:
            itr_flag = True
            itr_id_rm.append(iev)
    
        if posFOM_norm[iev] < posFOM_CUT:
            posFOM_flag = True
            posFOM_id_rm.append(iev)
        if posFOM_norm[iev] < posFOM_CUT and itr[iev] < ITR_CUT:
            both_flag = True
            both_id_rm.append(iev)

        if itr_flag == True and posFOM_flag == False:
            only_itr.append(iev)
        if itr_flag == False and posFOM_flag == True:
            only_posFOM.append(iev)
    print(f"ITR removed: {len(itr_id_rm)} ({(len(itr_id_rm)/len(itr))*100}) %")
    print(f"posFOM removed: {len(posFOM_id_rm)} ({(len(posFOM_id_rm)/len(posFOM_norm)) * 100}) % ")
    
    print(f"Only ITR: {len(only_itr)}")
    print(f"Only posFOM: {len(only_posFOM)}")
    print(f"Both Removed: {len(both_id_rm)} ({(len(both_id_rm)/len(posFOM_norm)) * 100}) %")

def compare_processing():
    """
    Function compares the distributions for 7.0.8 and 7.0.15 processed datasets.
    """

    # load the distributions
    x_708           = np.load("./full_analysis3_6m_x.npy")
    y_708           = np.load("./full_analysis3_6m_y.npy")
    z_708           = np.load("./full_analysis3_6m_z.npy")
    r_708           = np.load("./full_analysis3_6m_r.npy")
    energy_708      = np.load("./full_analysis3_6m_energy.npy")
    nhits_clean_708 = np.load("./full_analysis3_6m_nhits_clean.npy")
    hc_ratio_708    = np.load("./full_analysis3_6m_HC_ratio.npy")
    posFOM_708      = np.load("./full_analysis3_6m_posFOM.npy")
    posFOM_hits_708 = np.load("./full_analysis3_6m_posFOM_hits.npy")
    itr_708         = np.load("./full_analysis3_6m_itr.npy")
    dlogL_708       = np.load("./full_analysis3_6m_dlogL.npy")

    idx_keep = np.where(r_708 < FV_CUT)
    x_708           = x_708[idx_keep]
    y_708           = y_708[idx_keep]
    z_708           = z_708[idx_keep]
    r_708           = r_708[idx_keep]
    energy_708      = energy_708[idx_keep]
    nhits_clean_708 = nhits_clean_708[idx_keep]
    hc_ratio_708    = hc_ratio_708[idx_keep]
    posFOM_708      = posFOM_708[idx_keep]
    posFOM_hits_708 = posFOM_hits_708[idx_keep]
    itr_708         = itr_708[idx_keep]
    dlogL_708       = dlogL_708[idx_keep]
    rho2_708        = (x_708**2 + y_708**2) / 6000**2

    x_7015           = np.load("./full_analysis3_6m_x_reproc.npy")
    y_7015           = np.load("./full_analysis3_6m_y_reproc.npy")
    z_7015           = np.load("./full_analysis3_6m_z_reproc.npy")
    r_7015           = np.load("./full_analysis3_6m_r_reproc.npy")
    energy_7015      = np.load("./full_analysis3_6m_energy_reproc.npy")
    nhits_clean_7015 = np.load("./full_analysis3_6m_nhits_clean_reproc.npy")
    hc_ratio_7015    = np.load("./full_analysis3_6m_HC_ratio_reproc.npy")
    posFOM_7015      = np.load("./full_analysis3_6m_posFOM_reproc.npy")
    posFOM_hits_7015 = np.load("./full_analysis3_6m_posFOM_hits_reproc.npy")
    itr_7015         = np.load("./full_analysis3_6m_itr_reproc.npy")
    dlogL_7015       = np.load("./full_analysis3_6m_dlogL_reproc.npy")

    idx_keep = np.where(r_7015 < FV_CUT)
    x_7015           = x_7015[idx_keep]
    y_7015           = y_7015[idx_keep]
    z_7015           = z_7015[idx_keep]
    r_7015           = r_7015[idx_keep]
    energy_7015      = energy_7015[idx_keep]
    nhits_clean_7015 = nhits_clean_7015[idx_keep]
    hc_ratio_7015    = hc_ratio_7015[idx_keep]
    posFOM_7015      = posFOM_7015[idx_keep]
    posFOM_hits_7015 = posFOM_hits_7015[idx_keep]
    itr_7015         = itr_7015[idx_keep]
    dlogL_7015       = dlogL_7015[idx_keep]
    rho2_7015        = (x_7015**2 + y_7015**2) / 6000**2
    
    # ITR 
    plt.figure()
    plt.hist(itr_708, bins = itr_bins, density=True, histtype = "step", color = "black", label = f"7.0.8 | Num EVs: {len(itr_708)}")
    plt.hist(itr_7015, bins = itr_bins, density=True, histtype = "step", color = "red", label = f"7.0.15 | Num EVs: {len(itr_7015)}")
    plt.legend(frameon=False, prop = prop_font)
    plt.xlabel("ITR", fontproperties = prop_font)
    plt.ylabel("Normalised Counts", fontproperties = prop_font)
    plt.title(f"FV: {FV_CUT/1000:.1f} m", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.savefig(f"../plots/full_analysis3_data_quality/ITR_708_7015_{FV_CUT}.pdf")
    plt.close()

    # dlogL
    plt.figure()
    plt.hist(dlogL_708, bins = dlogL_bins, density=True, histtype = "step", color = "black", label = f"7.0.8 | Num EVs: {len(dlogL_708)}")
    plt.hist(dlogL_7015, bins = dlogL_bins, density=True, histtype = "step", color = "red", label = f"7.0.15 | Num EVs: {len(dlogL_7015)}")
    plt.legend(frameon=False, prop = prop_font)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font)
    plt.ylabel("Normalised Counts", fontproperties = prop_font)
    plt.title(f"FV: {FV_CUT/1000:.1f} m", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.savefig(f"../plots/full_analysis3_data_quality/dlogL_708_7015_{FV_CUT}.pdf")
    plt.close()

    # posFOM
    posFOM_708  /= posFOM_hits_708 
    posFOM_7015 /= posFOM_hits_7015
    plt.figure()
    plt.hist(posFOM_708, bins = posFOM_bins, density=True, histtype = "step", color = "black", label = f"7.0.8 | Num EVs: {len(dlogL_708)}")
    plt.hist(posFOM_7015, bins = posFOM_bins, density=True, histtype = "step", color = "red", label = f"7.0.15 | Num EVs: {len(dlogL_7015)}")
    plt.legend(frameon=False, prop = prop_font)
    plt.xlabel(r"posFOM / Nhits", fontproperties = prop_font)
    plt.ylabel("Normalised Counts", fontproperties = prop_font)
    plt.title(f"FV: {FV_CUT/1000:.1f} m", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.savefig(f"../plots/full_analysis3_data_quality/posFOM_708_7015_{FV_CUT}.pdf")
    plt.close()

    # x - y
    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    axes[0].hist2d(x_708, y_708, bins = [coordinate_bins, coordinate_bins], density = False, cmin = 1e-4)
    axes[0].set_title(f"7.0.8 | FV {FV_CUT/1000:.1f} m | Num EVs: {len(x_708)}", fontproperties = prop_font)
    axes[0].set_xlabel("Reconstructed X [mm]", fontproperties = prop_font)
    axes[0].set_ylabel("Reconstructed Y [mm]", fontproperties = prop_font)
    axes[0].set_aspect('equal')

    axes[1].hist2d(x_7015, y_7015, bins = [coordinate_bins, coordinate_bins], density = False, cmin = 1e-4)
    axes[1].set_title(f"7.0.15 | FV {FV_CUT/1000:.1f} m | Num EVs: {len(x_7015)}", fontproperties = prop_font)
    axes[1].set_xlabel("Reconstructed X [mm]", fontproperties = prop_font)
    axes[1].set_ylabel("Reconstructed Y [mm]", fontproperties = prop_font)
    axes[1].set_aspect('equal')
    fig.tight_layout()
    plt.savefig(f"../plots/full_analysis3_data_quality/X_Y_708_7015_{FV_CUT}.pdf")
    plt.close()

    # rho2 - z
    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    axes[0].hist2d(rho2_708, z_708, bins = [rho2_bins, coordinate_bins], density = False, cmin = 1e-4)
    axes[0].set_title(f"7.0.8 | FV {FV_CUT/1000:.1f} m | Num EVs: {len(x_708)}", fontproperties = prop_font)
    axes[0].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontproperties = prop_font)
    axes[0].set_ylabel("Reconstructed Z [mm]", fontproperties = prop_font)
    # axes[0].set_aspect('equal')

    axes[1].hist2d(rho2_7015, z_7015, bins = [rho2_bins, coordinate_bins], density = False, cmin = 1e-4)
    axes[1].set_title(f"7.0.15 | FV {FV_CUT/1000:.1f} m | Num EVs: {len(x_7015)}", fontproperties = prop_font)
    axes[1].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontproperties = prop_font)
    axes[1].set_ylabel("Reconstructed Z [mm]", fontproperties = prop_font)
    # axes[1].set_aspect('equal')
    fig.tight_layout()
    plt.savefig(f"../plots/full_analysis3_data_quality/rho2_Z_708_7015_{FV_CUT}.pdf")
    plt.close()

    # energy
    plt.figure()
    plt.hist(energy_708, bins = energy_bins, density=True, histtype = "step", color = "black", label = f"7.0.8 | Num EVs: {len(dlogL_708)}")
    plt.hist(energy_7015, bins = energy_bins, density=True, histtype = "step", color = "red", label = f"7.0.15 | Num EVs: {len(dlogL_7015)}")
    plt.legend(frameon=False, prop = prop_font)
    plt.xlabel(r"Reconstructed Energy [MeV]", fontproperties = prop_font)
    plt.ylabel("Normalised Counts", fontproperties = prop_font)
    plt.title(f"FV: {FV_CUT/1000:.1f} m", fontproperties = prop_font)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.savefig(f"../plots/full_analysis3_data_quality/energy_708_7015_{FV_CUT}.pdf")
    plt.close()

def reproc_in_detail():
    """
    Function loads run by run the 7.0.8 and 7.0.15 data and compares, GTID by GTID,
    the position and energy reconstructions results.
    """

    regenerate_matching = True
    energy_old      = []
    x_old           = []
    y_old           = []
    z_old           = []

    energy_new      = []
    x_new           = []
    y_new           = []
    z_new           = []

    energy_nomatch = []
    x_nomatch      = []
    y_nomatch      = []
    z_nomatch      = []
    gtid_nomatch   = []
    runnum_nomatch = []

    runlist     = np.loadtxt("../runlists/quiet_period.txt", dtype = int)
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3" 
    
    bad_runs = []
    total_evs = 0
    if regenerate_matching == True:
        for irun in runlist:
            # if irun in bad_runs:
            #     continue
            print(irun)

            # load the old and new files
            try:
                old_file   = ROOT.TFile.Open(f"{working_dir}/processed_dataset/{irun}.root")
                old_ntuple = old_file.Get("2p5_5p0")
                total_evs += old_ntuple.GetEntries()
                # print(f"{irun}, {old_ntuple.GetEntries()}")
            except:
                bad_runs.append(irun)
                continue
            new_file   = ROOT.TFile.Open(f"{working_dir}/reprocessed_dataset_7.0.15/{irun}.root")
            new_ntuple = new_file.Get("2p5_5p0")

            for ientry in old_ntuple:
                matched = False
                old_gtid = ientry.gtid

                for jentry in new_ntuple:

                    new_gtid = jentry.gtid

                    # find the matching event in the reprocessed dataset
                    if old_gtid == new_gtid:
                        energy_old.append(ientry.energy)
                        x_old.append(ientry.x)
                        y_old.append(ientry.y)
                        z_old.append(ientry.z)

                        energy_new.append(jentry.energy)
                        x_new.append(jentry.x)
                        y_new.append(jentry.y)
                        z_new.append(jentry.z)
                        matched = True
                    
                    if matched == True:
                        break

                # if we get here there's no matching in the reproc dataset
                if matched == False:
                    energy_nomatch.append(ientry.energy)
                    x_nomatch.append(ientry.x)
                    y_nomatch.append(ientry.y)
                    z_nomatch.append(ientry.z)
                    gtid_nomatch.append(ientry.gtid)
                    runnum_nomatch.append(irun)

            print(f"Found {len(x_old)} matching events.")
            print(f"Found {total_evs} in old set events.")
            print(f"Found {len(x_nomatch)} missing events in 7.0.15.")
        for irun in bad_runs:
            print(irun)
        print(f"Found {len(bad_runs)} bad runs.")

        # save these arrays
        np.save("./matched_energy_708.npy", energy_old)
        np.save("./matched_x_708.npy", x_old)
        np.save("./matched_y_708.npy", y_old)
        np.save("./matched_z_708.npy", z_old)
        np.save("./matched_energy_7015.npy", energy_new)
        np.save("./matched_x_7015.npy", x_new)
        np.save("./matched_y_7015.npy", y_new)
        np.save("./matched_z_7015.npy", z_new)
        np.save("./nomatched_energy_708.npy", energy_nomatch)
        np.save("./nomatched_x_708.npy", x_nomatch)
        np.save("./nomatched_y_708.npy", y_nomatch)
        np.save("./nomatched_z_708.npy", z_nomatch)
        np.save("./nomatched_gtid_708.npy", gtid_nomatch)
        np.save("./nomatched_runnum_708.npy", runnum_nomatch)

    energy_old = np.load("./matched_energy_708.npy")
    x_old      = np.load("./matched_x_708.npy")
    y_old      = np.load("./matched_y_708.npy")
    z_old      = np.load("./matched_z_708.npy")
    energy_new = np.load("./matched_energy_7015.npy")
    x_new      = np.load("./matched_x_7015.npy")
    y_new      = np.load("./matched_y_7015.npy")
    z_new      = np.load("./matched_z_7015.npy")

    no_match_x = np.load("./nomatched_x_708.npy")
    no_match_y = np.load("./nomatched_y_708.npy")
    no_match_z = np.load("./nomatched_z_708.npy")
    no_match_e = np.load("./nomatched_energy_708.npy")
    no_match_gtid = np.load("./nomatched_gtid_708.npy") 
    no_match_run = np.load("./nomatched_runnum_708.npy")
    

    r_nomatch = np.sqrt(no_match_x**2 + no_match_y**2 + no_match_z**2)
    idx_nomatch = np.where(r_nomatch < 4500)
    print(idx_nomatch)
    print("Number of no match events inside FV: ", len(idx_nomatch[0]))
    print(f"r: {r_nomatch[idx_nomatch]} | energy: {no_match_e[idx_nomatch]}")
    print(f"run num: {no_match_run[idx_nomatch]} | gtid: {no_match_gtid[idx_nomatch]}")
    r_old = np.sqrt(x_old**2 + y_old**2 + z_old**2)
    r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
    print(x_old[r_old < 4500].size)

    idx = np.where((r_old < 4500) & (r_new < 4500))
    print("Number of events found in both processings ROI + 6m which reconstruct in 4.5m: ", len(idx[0]))

    idx = np.where((r_old < 4500) & (r_new >= 4500))
    print("Number of events found in 708 which reconstruct outside in 7015: ", len(idx[0]))

    idx = np.where((r_old >= 4500) & (r_new < 4500))
    print("Number of events found outside in 708 which reconstruct inside in 7015: ", len(idx[0]))
    # interested in seeing where the events I lose from the FV go to
    
    
    x_pairs = np.array([x_old[r_old <= 4500], x_new[r_old <= 4500]])
    y_pairs = np.array([y_old[r_old <= 4500], y_new[r_old <= 4500]])
    import matplotlib.patches as patches
    # x - y line plot thing
    plt.scatter(x_old[r_old <= 5000], y_old[r_old <= 5000], color = "black", s = 1, label = "7.0.8")
    plt.scatter(x_new[r_old <= 5000], y_new[r_old <= 5000], color = "red", s = 1, label = "7.0.15")
    plt.plot(x_pairs, y_pairs, linestyle = "dashed", color = "black", linewidth = 1)    

    # add the FV circle to uide the eye
    circle = patches.Circle((0, 0), 4500, edgecolor="black", facecolor='none')
    ax = plt.gca()
    # Add the circle to the axes
    ax.add_patch(circle)
    
    # Set the aspect ratio to be equal
    ax.set_aspect('equal', adjustable='box')
    
    # Set limits to make sure circle is fully visible
    # ax.set_xlim(center[0] - radius - 1, center[0] + radius + 1)
    # ax.set_ylim(center[1] - radius - 1, center[1] + radius + 1)
    plt.legend(frameon = False)
    # plt.title(f"Matched Pairs: Old Processing Evs: {len(x_old[r_old < 4500])} | New Processing Evs: {len(x_new[r_new < 4500])}")
    plt.title("7.0.8 vs 7.0.15 Positions")
    plt.savefig("../plots/full_analysis3_data_quality/reproc_x_y_diffs.pdf")
    plt.close()

def inspect_specific_no_match_event():
    """
    Function looks at the specific events that are no longer present in the 6 m and 2.5 to 5.0 MeV ROI following
    reprocessing.
    """

    run_numbers = [307582, 308406]
    gtids       = [8619400, 14105664]
    run_numbers = [308406]
    gtids       = [14105664]

    for irun in range(len(run_numbers)):

        # we are going to open each reprocessed ratds file with dsreader, loop through it and extract the energy and position of the event
        fname = f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/reprocessed_ratds_7.0.15/{run_numbers[irun]}_4.root"
        # fname = f"/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/extracted_ratds/{run_numbers[irun]}*.root"
        # fnames = glob.glob(fname)
        # fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/condor/single_308406_reproc.root"
        print(fname)
        # for i in fnames:
        for ientry, _ in rat.dsreader(fname):

            # light path calculator and point3D stuff loaded after ratds constructor
            PMTCalStatus = RAT.DU.Utility.Get().GetPMTCalStatus()
            light_path = rat.utility().GetLightPathCalculator()
            group_velocity = rat.utility().GetGroupVelocity()
            pmt_info = rat.utility().GetPMTInfo()
            psup_system_id = RAT.DU.Point3D.GetSystemId("innerPMT")
            av_system_id = RAT.DU.Point3D.GetSystemId("av")
            
            #### RECONSTRUCTION INFORMATION EXTRACTED ####
            reconEvent = ientry.GetEV(0)

            # obtain the gtid of this event and see if it matches the target
            gtid = reconEvent.GetGTID()
            print(gtid)
            if gtid == gtids[irun]:
                # extract the position and energy and fit valid etc.
                # print("Found reprocessed event! SUBRUN FILE: ", i)
                fit_name = reconEvent.GetDefaultFitName()
                if not reconEvent.FitResultExists(fit_name):
                    continue

                vertex = reconEvent.GetFitResult(fit_name).GetVertex(0)
                reconPosition  = vertex.GetPosition() # returns in PSUP coordinates
                reconEnergy    = vertex.GetEnergy()        
                reconEventTime = vertex.GetTime()
                
                # apply AV offset to position
                event_point = RAT.DU.Point3D(psup_system_id, reconPosition)
                event_point.SetCoordinateSystem(av_system_id)
                radius = event_point.Mag()

                print(f"Radius: {radius} mm | Energy: {reconEnergy} MeV")
                break

def bipo_data_cleanliness_plots():
    """
    Function plots the energy, dR, rho2-Z and X-Y distributions of tagged BiPo212
    tagged for the Tl208 systematic constraint.
    """

    # load the file
    file   = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/systematic_tagged212/total.root")
    ntuple = file.Get("tag_info")

    num          = ntuple.GetEntries()
    bi_positions = np.zeros((num, 3))
    po_positions = np.zeros((num, 3))
    dR           = np.zeros(num)
    bi_energy    = np.zeros(num)
    po_energy    = np.zeros(num)

    count = 0
    for ientry in ntuple:
        bi_energy[count] = ientry.energy_bi
        po_energy[count] = ientry.energy_po
        dR[count] = ientry.dR
        bi_positions[count, 0] = ientry.x_bi
        bi_positions[count, 1] = ientry.y_bi
        bi_positions[count, 2] = ientry.z_bi
        po_positions[count, 0] = ientry.x_po
        po_positions[count, 1] = ientry.y_po
        po_positions[count, 2] = ientry.z_po
        count += 1
    energy_binning     = np.arange(0, 4.0, 0.1)
    coordinate_binning = np.arange(-6000, 6050, 50)

    rho2_binning       = np.arange(0, 1.2, 0.2)
    dR_binning         = np.arange(0, 1020, 20) 

    idx = np.where(bi_energy > 2.1)

    ## energy ##
    plt.figure()
    plt.hist(bi_energy, bins = energy_binning, histtype = "step", label = "Bi212")
    plt.hist(po_energy, bins = energy_binning, histtype = "step", label = "Po212")
    plt.xlabel("Reconstructed Energy [MeV]")
    plt.legend()
    plt.savefig("../plots/full_analysis3_data_quality/bipo212_energy.pdf")
    plt.close()

    ## dR ##
    plt.figure()
    plt.hist(dR[idx], bins = dR_binning, histtype = "step")
    plt.xlabel("dR")
    plt.legend()
    plt.savefig("../plots/full_analysis3_data_quality/bipo212_dR.pdf")
    plt.close()

    # X-Y
    x_coords = np.array([bi_positions[:,0], po_positions[:,0]])
    y_coords = np.array([bi_positions[:,1], po_positions[:,1]])
    plt.figure()
    # plt.plot(x_coords, y_coords, linestyle = "--", color = "black")
    plt.scatter(bi_positions[idx][:,0], bi_positions[idx][:,1], color = "blue", s = 2)
    plt.scatter(po_positions[idx][:,0], po_positions[idx][:,1], color = "red", s = 2)
    plt.xlabel("Reconstructed X [mm]")
    plt.ylabel("Reconstructed Y [mm]")
    plt.xlim((-6000, 6000))
    plt.ylim((-6000, 6000))
    plt.savefig("../plots/full_analysis3_data_quality/bipo212_X_Y.pdf")
    plt.close()


    ## work out the number of BiPo tagged in 4.5 m, 4.5 +- 0.35 m systematic
    # using just the Bi214 position
    r_bi    = np.sqrt(np.sum(bi_positions[:,0:3]**2, axis = 1))
    idx_fv  = np.where(r_bi < 5000)
    idx_pFV = np.where(r_bi < 4500 + 35)
    idx_nFV = np.where(r_bi < 4500 - 35)

    print(f"in FV: {len(idx_fv[0])}\n+ve: {len(idx_pFV[0])}\n-ve: {len(idx_nFV[0])}")  
    # rho2 z
    # print(len(idx[0]))
    # bi_positions = bi_positions[idx,:]
    # rho2 = (bi_positions[:,0]**2 + bi_positions[:,1]**2) / (6000**2)
    # plt.figure()
    # plt.hist2d(rho2, bi_positions[:,2], bins = [rho2_binning, coordinate_binning], cmin = 1e-4)
    # plt.xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$")
    # plt.ylabel("Z [mm]")
    # plt.savefig("../plots/full_analysis3_data_quality/bipo212_rho2_Z.pdf")



bipo_data_cleanliness_plots()
# reproc_in_detail()                
# inspect_specific_no_match_event()


FV_CUT     = 4500
ITR_CUT    = 0.22
posFOM_CUT = 13.7
save_dir   = "../plots/full_analysis3_data_quality"

## load up the quantities ## 
# x           = np.load("./full_analysis3_6m_x.npy")
# y           = np.load("./full_analysis3_6m_y.npy")
# z           = np.load("./full_analysis3_6m_z.npy")
# r           = np.load("./full_analysis3_6m_r.npy")
# energy      = np.load("./full_analysis3_6m_energy.npy")
# nhits_clean = np.load("./full_analysis3_6m_nhits_clean.npy")
# hc_ratio    = np.load("./full_analysis3_6m_HC_ratio.npy")
# posFOM      = np.load("./full_analysis3_6m_posFOM.npy")
# posFOM_hits = np.load("./full_analysis3_6m_posFOM_hits.npy")
# itr         = np.load("./full_analysis3_6m_itr.npy")
# dlogL       = np.load("./full_analysis3_6m_dlogL.npy")
load_event_info()
x           = np.load("./full_analysis3_6m_x_reproc.npy")
y           = np.load("./full_analysis3_6m_y_reproc.npy")
z           = np.load("./full_analysis3_6m_z_reproc.npy")
r           = np.load("./full_analysis3_6m_r_reproc.npy")
energy      = np.load("./full_analysis3_6m_energy_reproc.npy")
nhits_clean = np.load("./full_analysis3_6m_nhits_clean_reproc.npy")
hc_ratio    = np.load("./full_analysis3_6m_HC_ratio_reproc.npy")
posFOM      = np.load("./full_analysis3_6m_posFOM_reproc.npy")
posFOM_hits = np.load("./full_analysis3_6m_posFOM_hits_reproc.npy")
itr         = np.load("./full_analysis3_6m_itr_reproc.npy")
dlogL       = np.load("./full_analysis3_6m_dlogL_reproc.npy")

# apply FV cut - return idx of events to keep in each array
idx_keep = np.where(r < FV_CUT)

x           = x[idx_keep]
y           = y[idx_keep]
z           = z[idx_keep]
r           = r[idx_keep]
energy      = energy[idx_keep]
nhits_clean = nhits_clean[idx_keep]
hc_ratio    = hc_ratio[idx_keep]
posFOM      = posFOM[idx_keep]
posFOM_hits = posFOM_hits[idx_keep]
itr         = itr[idx_keep]
dlogL       = dlogL[idx_keep]

# rho2        = (x**2 + y**2) / 6000**2

# # apply the cuts
# idx_keep = np.where((itr > 0.2) & (itr < 0.3) & ((posFOM/posFOM_hits) > 13.7))
# x           = x[idx_keep]
# y           = y[idx_keep]
# z           = z[idx_keep]
# r           = r[idx_keep]
# energy      = energy[idx_keep]
# nhits_clean = nhits_clean[idx_keep]
# hc_ratio    = hc_ratio[idx_keep]
# posFOM      = posFOM[idx_keep]
# posFOM_hits = posFOM_hits[idx_keep]
# itr         = itr[idx_keep]
# dlogL       = dlogL[idx_keep]

# rho2        = (x**2 + y**2) / 6000**2

# define the binning
if FV_CUT == 4500:
    factor = 5
else:
    factor = 1

energy_bins = np.arange(2.5, 5.0, 0.01 * factor)
itr_bins    = np.arange(0.1, 0.5, 0.005 * factor/3)
hc_bins     = np.arange(0.8, 1, 0.005)
posFOM_bins = np.arange(13, 15, 0.01 * factor)
dlogL_bins  = np.arange(-1.24, -1.185, 0.0002*factor)
coordinate_bins = np.arange(-6000, 6050, 50*factor)
rho2_bins = np.arange(0, 1, 0.01*factor)
radial_bins = np.arange(0, 6050, 50*factor)

# compare_processing()
# event_cuts(ITR_CUT, posFOM_CUT)

# create_1D_plots(FV_CUT)
# create_2D_plots(x, y, coordinate_bins, coordinate_bins, "X", "Y")
# create_2D_plots(rho2, z, rho2_bins, coordinate_bins, "rho2", "Z")
# create_2D_plots(x, itr, coordinate_bins, itr_bins, "X", "ITR")
# create_2D_plots(y, itr, coordinate_bins, itr_bins, "Y", "ITR")
# create_2D_plots(z, itr, coordinate_bins, itr_bins, "Z", "ITR")
# create_2D_plots(r, itr, radial_bins, itr_bins, "R", "ITR")

# create_2D_plots(x, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "X", "posFOM_div_posFOM_hits")
# create_2D_plots(y, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Y", "posFOM_div_posFOM_hits")
# create_2D_plots(z, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Z", "posFOM_div_posFOM_hits")
# create_2D_plots(r, posFOM/posFOM_hits, radial_bins, posFOM_bins, "R", "posFOM_div_posFOM_hits")

# create_2D_plots(x, hc_ratio, coordinate_bins, hc_bins, "X", "HC_ratio")
# create_2D_plots(y, hc_ratio, coordinate_bins, hc_bins, "Y", "HC_ratio")
# create_2D_plots(z, hc_ratio, coordinate_bins, hc_bins, "Z", "HC_ratio")
# create_2D_plots(r, hc_ratio, radial_bins, hc_bins, "R", "HC_ratio")

# create_2D_plots(itr, hc_ratio, itr_bins, hc_bins, "ITR", "HC_ratio")
# create_2D_plots(itr, posFOM/posFOM_hits, itr_bins, posFOM_bins, "ITR", "posFOM_div_posFOM_hits")
# create_2D_plots(hc_ratio, posFOM/posFOM_hits, hc_bins, posFOM_bins, "HC_ratio", "posFOM_div_posFOM_hits")


# create_1D_plots_postCuts(FV_CUT)
# create_2D_plots_postCuts(x, y, coordinate_bins, coordinate_bins, "X", "Y")
# create_2D_plots_postCuts(rho2, z, rho2_bins, coordinate_bins, "rho2", "Z")
# create_2D_plots_postCuts(x, itr, coordinate_bins, itr_bins, "X", "ITR")
# create_2D_plots_postCuts(y, itr, coordinate_bins, itr_bins, "Y", "ITR")
# create_2D_plots_postCuts(z, itr, coordinate_bins, itr_bins, "Z", "ITR")
# create_2D_plots_postCuts(r, itr, radial_bins, itr_bins, "R", "ITR")

# create_2D_plots_postCuts(x, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "X", "posFOM_div_posFOM_hits")
# create_2D_plots_postCuts(y, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Y", "posFOM_div_posFOM_hits")
# create_2D_plots_postCuts(z, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Z", "posFOM_div_posFOM_hits")
# create_2D_plots_postCuts(r, posFOM/posFOM_hits, radial_bins, posFOM_bins, "R", "posFOM_div_posFOM_hits")

# create_2D_plots_postCuts(x, hc_ratio, coordinate_bins, hc_bins, "X", "HC_ratio")
# create_2D_plots_postCuts(y, hc_ratio, coordinate_bins, hc_bins, "Y", "HC_ratio")
# create_2D_plots_postCuts(z, hc_ratio, coordinate_bins, hc_bins, "Z", "HC_ratio")
# create_2D_plots_postCuts(r, hc_ratio, radial_bins, hc_bins, "R", "HC_ratio")

# create_2D_plots_postCuts(itr, hc_ratio, itr_bins, hc_bins, "ITR", "HC_ratio")
# create_2D_plots_postCuts(itr, posFOM/posFOM_hits, itr_bins, posFOM_bins, "ITR", "posFOM_div_posFOM_hits")
# create_2D_plots_postCuts(hc_ratio, posFOM/posFOM_hits, hc_bins, posFOM_bins, "HC_ratio", "posFOM_div_posFOM_hits")