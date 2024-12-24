import numpy as np 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT
import matplotlib.font_manager as fm

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
Creates a plot of the energy dependence of the multisite classifier.
"""
def energy_dependence_of_classifier_plot():
    # load the root files containing the distributions
    tl208_dists = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/full_analysis3_Tl208/total.root")
    b8_dists = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/full_analysis3_B8_solar_nue/total.root")

    # load the multisite PDFs
    ranges = ["2p5_3p0", "3p0_3p5", "3p5_4p0", "4p0_4p5", "4p5_5p0"]
    colours = ["orange", "green", "blue", "red", "black"]
    fig, axes = plt.subplots(nrows = 1, ncols = 5, figsize = (13, 3))
    #binning = np.linspace(-2.12, -2.06, 100)#np.linspace(-0.015, 0.006, 100)
    binning = 100#np.arange(-1.375, -1.325, 0.0005)
    for i in range(len(ranges)):
        # load each distribution in turn
        b8 = b8_dists.Get(f"{ranges[i]}")
        tl208 = tl208_dists.Get(f"{ranges[i]}")
        dist_b8 = []
        dist_tl = []
        for ievent in b8:
            dlogl = ievent.dlogL
            dist_b8.append(dlogl)
        for ievent in tl208:
            dlogl = ievent.dlogL
            dist_tl.append(dlogl)

        # shift each disto so the mean is centred on 0 for display
        mean_tl = np.mean(dist_tl)
        mean_b8 = np.mean(dist_b8)

        print(f"{ranges[i]} - min: {np.min(dist_b8)}, {np.min(dist_tl)} | max: {np.max(dist_b8)}, {np.max(dist_tl)}")

        axes[i].hist(dist_b8, bins = binning, linewidth = 2, color = "black", histtype = "step", density = True)
        axes[i].hist(dist_tl, bins = binning, linewidth = 2, color = "red", linestyle = "dotted", histtype= "step", density = True)
        axes[i].set_xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font)
        axes[i].legend(prop=fm.FontProperties(fname=path2, size = 14), handles = [plt.plot([], [], color = "black", linewidth = 2)[0], plt.plot([], [], color = "red", linestyle = "dotted", linewidth =2)[0]], labels = [r"$^8$B", r"$^{208}$Tl"], frameon = False, loc = "upper left")
        for label in axes[i].get_xticklabels():
            label.set_fontproperties(prop_font)

        for label in axes[i].get_yticklabels():
            label.set_fontproperties(prop_font)

    plt.tight_layout()
    plt.savefig("multisite_e_dependence2.pdf")


energy_dependence_of_classifier_plot()