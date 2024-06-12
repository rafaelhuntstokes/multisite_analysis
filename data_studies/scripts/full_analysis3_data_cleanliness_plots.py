import ROOT
import numpy as np
import matplotlib
import matplotlib.font_manager as fm
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    plt.savefig(f"{save_dir}/energy_{FV_CUT}.pdf")
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
    plt.savefig(f"{save_dir}/itr_{FV_CUT}.pdf")
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
    plt.savefig(f"{save_dir}/hc_ratio_{FV_CUT}.pdf")
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
    plt.savefig(f"{save_dir}/posFOM_{FV_CUT}.pdf")
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
    plt.savefig(f"{save_dir}/dlogL_{FV_CUT}.pdf")
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
    plt.savefig(f"{save_dir}/{name1}_{name2}_{FV_CUT}.pdf")
    plt.close() 

FV_CUT   = 4500
save_dir = "../plots/full_analysis3_data_quality"

## load up the quantities ## 
x           = np.load("./full_analysis3_6m_x.npy")
y           = np.load("./full_analysis3_6m_y.npy")
z           = np.load("./full_analysis3_6m_z.npy")
r           = np.load("./full_analysis3_6m_r.npy")
energy      = np.load("./full_analysis3_6m_energy.npy")
nhits_clean = np.load("./full_analysis3_6m_nhits_clean.npy")
hc_ratio    = np.load("./full_analysis3_6m_HC_ratio.npy")
posFOM      = np.load("./full_analysis3_6m_posFOM.npy")
posFOM_hits = np.load("./full_analysis3_6m_posFOM_hits.npy")
itr         = np.load("./full_analysis3_6m_itr.npy")
dlogL       = np.load("./full_analysis3_6m_dlogL.npy")

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

rho2        = (x**2 + y**2) / 6000**2

# define the binning
if FV_CUT == 4500:
    factor = 5
else:
    factor = 1

energy_bins = np.arange(2.5, 5.0, 0.01)
itr_bins    = np.arange(0.1, 0.5, 0.005)
hc_bins     = np.arange(0.8, 1, 0.005)
posFOM_bins = np.arange(13, 15, 0.01)
dlogL_bins  = np.arange(-1.24, -1.185, 0.0002*factor)
coordinate_bins = np.arange(-6000, 6050, 50*factor)
rho2_bins = np.arange(0, 1, 0.01*factor)
radial_bins = np.arange(0, 6050, 50*factor)

# load_event_info()
create_1D_plots(FV_CUT)
create_2D_plots(x, y, coordinate_bins, coordinate_bins, "X", "Y")
create_2D_plots(rho2, z, rho2_bins, coordinate_bins, "rho2", "Z")
create_2D_plots(x, itr, coordinate_bins, itr_bins, "X", "ITR")
create_2D_plots(y, itr, coordinate_bins, itr_bins, "Y", "ITR")
create_2D_plots(z, itr, coordinate_bins, itr_bins, "Z", "ITR")
create_2D_plots(r, itr, radial_bins, itr_bins, "R", "ITR")

create_2D_plots(x, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "X", "posFOM_div_posFOM_hits")
create_2D_plots(y, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Y", "posFOM_div_posFOM_hits")
create_2D_plots(z, posFOM/posFOM_hits, coordinate_bins, posFOM_bins, "Z", "posFOM_div_posFOM_hits")
create_2D_plots(r, posFOM/posFOM_hits, radial_bins, posFOM_bins, "R", "posFOM_div_posFOM_hits")

create_2D_plots(x, hc_ratio, coordinate_bins, hc_bins, "X", "HC_ratio")
create_2D_plots(y, hc_ratio, coordinate_bins, hc_bins, "Y", "HC_ratio")
create_2D_plots(z, hc_ratio, coordinate_bins, hc_bins, "Z", "HC_ratio")
create_2D_plots(r, hc_ratio, radial_bins, hc_bins, "R", "HC_ratio")

create_2D_plots(itr, hc_ratio, itr_bins, hc_bins, "ITR", "HC_ratio")
create_2D_plots(itr, posFOM/posFOM_hits, itr_bins, posFOM_bins, "ITR", "posFOM_div_posFOM_hits")
create_2D_plots(hc_ratio, posFOM/posFOM_hits, hc_bins, posFOM_bins, "HC_ratio", "posFOM_div_posFOM_hits")
