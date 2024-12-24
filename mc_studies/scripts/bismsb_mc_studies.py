import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
from scipy import stats 
"""
Script looks at the new optics table MC and checks the resolution / bias of the reconstruction.

This is to be compared with the bias / resolution of the 2.2 g/l PPO full-fill phase MC.
"""
path = '/data/snoplus3/hunt-stokes/nuPhysPoster/scripts/Times_New_Roman_Normal.ttf'
prop_font = fm.FontProperties(fname=path, size = 24)

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

def func(x, a, x0, sigma): 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 

def load_datasets(file_name, FV_CUT):
    """
    Function takes in a set of ntuples and extracts all the event reconstructed
    and true position and energy. Returns .npy arrays.
    """
    
    file   = ROOT.TFile.Open(file_name)
    ntuple = file.Get("output")

    dX = []
    dY = []
    dZ = []
    dE = []

    # some parameters of interest to compare between optics models 
    E_rec = []
    nhits = [] 
    for ientry in ntuple:
        
        if ientry.fitValid == False or ientry.scintFit == False:
            continue
        x_rec = ientry.posx
        x_tru = ientry.mcPosx

        y_rec = ientry.posy
        y_tru = ientry.mcPosy

        z_rec = ientry.posz
        z_tru = ientry.mcPosz

        r_tru = np.sqrt(x_tru**2 + y_tru**2 + z_tru**2)
        if r_tru > FV_CUT:
            continue
        e_rec = ientry.energy
        e_tru = ientry.mcke1

        if e_tru < E_LOW or e_tru > E_HIGH:
            continue

        # save the differences
        dX.append(x_tru - x_rec)
        dY.append(y_tru - y_rec) 
        dZ.append(z_tru - z_rec)
        dE.append(e_tru - e_rec)
        E_rec.append(e_rec)
        nhits.append(ientry.correctedNhits)

    return dX, dY, dZ, dE, E_rec, nhits

def res_plot_maker(bis, ppo, binning, unit, fit_range, quantity):

    plt.figure()



    plt.hist(bis, bins = binning, density = True, linewidth = 2, histtype = "step", color = "red", label = f"Bismsb Optics\nBias = {np.mean(bis):.2f} {unit}\nRes. = {np.std(bis):.2f} {unit}")
    plt.hist(ppo, bins = binning, density = True, linewidth = 2, histtype = "step", color = "black", label = f"PPO Optics\nBias = {np.mean(ppo):.2f} {unit}\nRes. = {np.std(ppo):.2f} {unit}")

    fit   = stats.norm.pdf(fit_range, np.mean(bis), np.std(bis))
    plt.plot(fit_range, fit, color = "red", linestyle = "dashed")
    fit   = stats.norm.pdf(fit_range, np.mean(ppo), np.std(ppo))
    plt.plot(fit_range, fit, color = "black", linestyle = "dashed")

    plt.legend(fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 15),frameon=False)
    plt.xlabel(rf"${quantity}_{{True}} - {quantity}_{{Reco.}}$", fontproperties = prop_font)
    ax = plt.gca()

    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    ax.text(0.05, 0.70, r"$^8$B - $\nu_e$ MC" + f"\n\nFV {FV_CUT/1000:.1f} m\n\n" + rf"{E_LOW} $\leq$ $E_{{reco}}$ $\leq$ {E_HIGH} MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    plt.tight_layout()
    # plt.yscale("log")
    plt.savefig(f"../plots/bismsb_optics/{quantity}_{FV_CUT}_{E_LOW}_{E_HIGH}.pdf")
    plt.close()


# arrays to save the distributions
ppo_E, ppo_nhits, ppo_dX, ppo_dY, ppo_dZ, ppo_dE = [], [], [], [], [], []
bis_E, bis_nhits, bis_dX, bis_dY, bis_dZ, bis_dE = [], [], [], [], [], []
FV_CUT = 4500
E_LOW  = 2.5
E_HIGH = 3.5
# load the 2.2 PPO full-fill data
runlist_ppo = np.loadtxt("../../data_studies/runlists/simulated_list.txt", dtype = int)

for irun in runlist_ppo:

    file = f"/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationB8_solar_nue_{irun}.ntuple.root"

    try:
        dX, dY, dZ, dE, E_rec, nhits = load_datasets(file, FV_CUT)
    except:
        continue
    
    ppo_E += E_rec
    ppo_nhits += nhits
    ppo_dX += dX
    ppo_dY += dY
    ppo_dZ += dZ
    ppo_dE += dE

# load the 2.2 PPO + bismsb data
runlist_bis = np.loadtxt("../runlists/bismsb_mc_runlist.txt", dtype = int)

for irun in runlist_bis:

    file = f"/data/snoplus2/miniPROD_bismsb_newOptics_oldFitter7015/ntuples/B8_solar_nue_{irun}_bismsb.ntuple.root"

    # try:
    dX, dY, dZ, dE, E_rec, nhits = load_datasets(file, FV_CUT)
    # except:
        # continue
    
    bis_E += E_rec
    bis_nhits += nhits
    bis_dX += dX
    bis_dY += dY
    bis_dZ += dZ
    bis_dE += dE

coord_bins = np.arange(-500, 500, 15)
coord_mids = coord_bins[:-1] + np.diff(coord_bins)[0]/2


# now we can make the plots

res_plot_maker(bis_dX, ppo_dX, coord_bins, "mm", np.linspace(-500, 500, 1000), "X")
res_plot_maker(bis_dY, ppo_dY, coord_bins, "mm", np.linspace(-500, 500, 1000), "Y")
res_plot_maker(bis_dZ, ppo_dZ, coord_bins, "mm", np.linspace(-500, 500, 1000), "Z")

energy_bins = np.arange(-2.0, 2.0, 0.1)
res_plot_maker(bis_dE, ppo_dE, energy_bins, "MeV", np.linspace(-2, 2, 1000), "Energy")

plt.figure()
plt.hist(bis_nhits, bins = 100,  density = True, linewidth = 2, histtype = "step", color = "red", label = f"Bismsb Optics")
plt.hist(ppo_nhits, bins = 100, density = True, linewidth = 2, histtype = "step", color = "black", label = f"PPO Optics")
plt.legend()
plt.xlabel("Corrected Nhits")
plt.savefig("../plots/bismsb_optics/nhits_true.pdf")
plt.close()

plt.figure()
plt.hist(bis_E, bins = 100,  density = True, linewidth = 2, histtype = "step", color = "red", label = f"Bismsb Optics")
plt.hist(ppo_E, bins = 100, density = True, linewidth = 2, histtype = "step", color = "black", label = f"PPO Optics")
plt.legend()
plt.xlabel("Reconstructed Energy")
plt.savefig("../plots/bismsb_optics/energy_true.pdf")