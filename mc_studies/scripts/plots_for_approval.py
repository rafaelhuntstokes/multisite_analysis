import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import seaborn as sns

"""
USE env_poster in nuPhysPoster!
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

def create_time_residual_pdfs():
    """
    Create a plot of the time residuals for B8 and Tl208 MC between 2.5 and 
    5.0 MeV.
    """
    ## create plot ##
    # cmap = plt.get_cmap("inferno")
    # col1 = cmap(0.1)
    # col2 = cmap(0.4)
    # col3 = cmap(0.7)
    # col4 = cmap(0.0)

    map = "colorblind"
    cmap = sns.color_palette(map, 3)
    col1 = cmap[1]
    col2 = cmap[2]
    # col3 = cmap[2]
    # col4 = cmap[3]
    # col5 = cmap[4]

    pdf_B8_file    = ROOT.TFile.Open("../run_by_run_pdf/full_analysis2_B8_solar_nue/total.root")
    pdf_Tl208_file = ROOT.TFile.Open("../run_by_run_pdf/full_analysis2_Tl208/total.root")

    # the time residuals are already binned --> so making a bar graph
    multisite_pdf_B8    = pdf_B8_file.Get("multi_2.5_5.0")
    multisite_pdf_Tl208 = pdf_Tl208_file.Get("multi_2.5_5.0")

    multi_B8_bins      = [multisite_pdf_B8.GetBinLowEdge(i) for i in range(1, multisite_pdf_B8.GetNbinsX() + 2)]
    multi_B8_counts    = [multisite_pdf_B8.GetBinContent(i) / multisite_pdf_B8.Integral() for i in range(1, multisite_pdf_B8.GetNbinsX()+1)]
    multi_Tl208_counts = [multisite_pdf_Tl208.GetBinContent(i) / multisite_pdf_Tl208.Integral() for i in range(1, multisite_pdf_Tl208.GetNbinsX()+1)]


    fig = plt.figure()
    b8 = plt.step(multi_B8_bins[:-1], multi_B8_counts,  where = 'post', linewidth = 2, color = col1, label = r"$^8B$")
    tl208 = plt.step(multi_B8_bins[:-1], multi_Tl208_counts,  where = 'post', linewidth = 2, linestyle = "dashed", alpha = 1, color = col2, label = r"$^{208}Tl$")
    plt.legend(handles = [plt.plot([], [], color = col1)[0], plt.plot([], [], linestyle = "dashed", color = col2)[0]], labels = [r"${^8}$B", r"$^{208}$Tl"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    plt.xlabel("Time Residual [ns]", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Normalised Counts", fontproperties = prop_font, fontsize = 26)
    # plt.title(r"Time Residuals for $^8$B and $^{208}$Tl Events", fontproperties = prop_font, fontsize = 26)
    plt.xlim((-5, 100))
    plt.ylim((0, 0.05))
    ax = plt.gca()
    ax.xaxis.set_minor_locator(MultipleLocator(5))

    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    ax.text(0.55, 0.35, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/time_residuals_{map}_coloured_2p5_5p0.pdf")
    plt.close()

def create_dlogL_plot():
    """
    Create a plot of the dlogL multisite discriminant for 2.5 to 5.0 MeV between
    B8 and Tl208.
    
    Add an 'inset' plot showing the ROC curve.
    """
    
    # cmap = plt.get_cmap("inferno")
    # col1 = cmap(0.1)
    # col2 = cmap(0.4)
    # col3 = cmap(0.7)
    # col4 = cmap(0.0)

    map = "colorblind"
    cmap = sns.color_palette(map, 3)
    col1 = cmap[1]
    col2 = cmap[2]
    # col3 = cmap[2]
    # col4 = cmap[3]
    # col5 = cmap[4]

    # find the multisite PDFs for this energy
    pdfs      = np.load("multisite_pdf_array_2p5_5p0.npy")

    # extract the two of interest
    b8_pdf    = pdfs[0,:]
    tl208_pdf = pdfs[1,:]
    print(np.sum(b8_pdf), np.sum(tl208_pdf))

    binning_full = np.arange(-1.375, -1.325, 0.0005)
    binning_3p5_5p0 = np.arange(-1.309, -1.26, 0.0005)
    fig = plt.figure()
    # plt.step(binning_3p5_5p0, b8_pdf.tolist() + [0], where = "post", linewidth = 2, color = col1)
    # plt.step(binning_3p5_5p0, tl208_pdf.tolist() + [0], where = "post", linestyle = "dashed", linewidth = 2, color = col2)
    plt.step(binning_full, b8_pdf.tolist() + [0], where = "post", linewidth = 2, color = col1)
    plt.step(binning_full, tl208_pdf.tolist() + [0], where = "post", linestyle = "dashed", linewidth = 2, color = col2)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Normalised Counts", fontproperties = prop_font, fontsize = 26)
    # plt.title(r"Multisite Discriminants for $^8$B and $^{208}$Tl", fontproperties = prop_font, fontsize = 26)
    plt.ylim((0, 0.12))

    plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], color = col2, linewidth = 2, linestyle = "dashed")[0]], labels = [r"${^8}$B", r"$^{208}$Tl"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    ax.text(0.125, 0.6, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)

    #### create ROC curve from the PDFs ####
    signal_acceptance     = []
    background_acceptance = []

    for ibin in range(len(binning_3p5_5p0) - 1):
        signal_accepted = np.sum(b8_pdf[ibin:]) 
        background_accepted = np.sum(tl208_pdf[ibin:]) 
        signal_acceptance.append(signal_accepted)
        background_acceptance.append(background_accepted)
    
    # create the ROC curve inset graph
    factor = 1.5
    left, bottom, width, height = [0.25, 0.28, 0.2*factor, 0.2*factor]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(signal_acceptance, background_acceptance, color = "black", linewidth = 2)
    ax2.set_xlim((0,1))
    ax2.set_ylim((0,1))
    ax2.set_xlabel(r"$^8$B Acceptance", fontproperties = prop_font, fontsize = 15)
    ax2.set_ylabel(r"$^{208}$Tl Acceptance", fontproperties = prop_font, fontsize = 15)
    for label in ax2.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax2.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)

    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/multisite_pdfs_{map}_coloured_2p5_5p0.pdf")
    plt.close()

def create_asimov_dataset_graphs():
    """
    Function creates the asimov dataset graphs in terms of energy and multisite
    PDFs. The function creates a combined subplot style, and individual plots.
    """
   
    # cmap = plt.get_cmap("inferno")
    # col1 = cmap(0.0)
    # col2 = cmap(0.2)
    # col3 = cmap(0.4)
    # col4 = cmap(0.7)
    # col5 = cmap(0.58)
    map = "colorblind"
    cmap = sns.color_palette(map, 8)
    col1 = cmap[1]
    col2 = cmap[2]
    col3 = cmap[7]
    col4 = cmap[3]
    col5 = cmap[4]

    # load the pdfs in terms of energy and multisite discriminant
    # pdfs_multisite = np.load("multisite_pdf_array_3p5_5p0.npy")
    # pdfs_energy    = np.load("energy_pdf_array_3p5_5p0.npy")
    pdfs_multisite = np.load("multisite_pdf_array_2p5_5p0.npy")
    pdfs_energy    = np.load("energy_pdf_array_2p5_5p0.npy")

    # define the normalisations of each - from background model and livetime
    normalisations = np.array([66.3, 468.9, 0.9081, 59.227, 38.528])
    weights        = np.array([0.551, 0.634, 0.620, 0, 0])
    weights        = np.array([1, 1, 1, 1, 1])
    normalisations = normalisations * weights

    # define everything to the solar signal, which has norm 1
    normalisations = normalisations / normalisations[0]

    # multiply each PDF by the respective normalisations 
    scaled_energy    = normalisations[:, None] * pdfs_energy
    scaled_multisite = normalisations[:, None] * pdfs_multisite

    energy_bins            = np.arange(2.5, 5.05, 0.05)
    multisite_bins         = np.arange(-1.375, -1.325, 0.0005)
    multisite_bins_3p5_5p0 = np.arange(-1.309, -1.26, 0.0005)
    energy_bins_3p5_5p0    = np.arange(2.5, 5.05, 0.05)
    # create the individual plots
    
    ### energy plot ###
    fig = plt.figure()
    # plt.step(energy_bins_3p5_5p0, scaled_energy[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    # plt.step(energy_bins_3p5_5p0, scaled_energy[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # # plt.step(energy_bins, scaled_energy[2, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # # plt.step(energy_bins_3p5_5p0, scaled_energy[3, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # # plt.step(energy_bins_3p5_5p0, scaled_energy[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10,2], linewidth = 2)
    # plt.step(energy_bins_3p5_5p0, np.sum(scaled_energy, axis= 0).tolist() + [0], color = col3, linestyle = "dotted", where = 'post', linewidth = 2)

    plt.step(energy_bins, scaled_energy[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    plt.step(energy_bins, scaled_energy[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(energy_bins, scaled_energy[2, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    plt.step(energy_bins, scaled_energy[3, :].tolist() + [0], color = col5, where = 'post', dashes = [5, 2], linewidth = 2)
    plt.step(energy_bins, scaled_energy[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10,2], linewidth = 2)
    plt.step(energy_bins, np.sum(scaled_energy, axis= 0).tolist() + [0], color = col3, linestyle = "dotted", where = 'post', linewidth = 2)


    plt.xlabel("Reconstructed Energy [MeV]", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Relative Counts", fontproperties = prop_font, fontsize = 26)
    # plt.title("Asimov Dataset", fontproperties= prop_font, fontsize = 26)
    # plt.xlim((3.5, 5.0))
    plt.xlim((2.5, 5.0))
    # plt.ylim((0, 0.8))
    plt.ylim((0, 0.4))
    
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], color = col2, linewidth = 2, linestyle = "dashdot")[0],  \
                          plt.plot([], [], color = col5, linewidth = 2, dashes = [5, 2])[0],\
                          plt.plot([], [], color = col4, linewidth = 2, dashes = [10, 2])[0], plt.plot([], [], color = col3, linewidth = 2)[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "BiPo212", "BiPo214", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    # plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], color = col2, linewidth = 2, linestyle = "dashdot")[0],  \
    #                       plt.plot([], [], color = col3, linewidth = 2, linestyle = "dotted")[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)

    # ax.text(0.66, 0.42, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    ax.text(0.06, 0.72, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/asimov_energy_seaborn_{map}_coloured_2p5_5p0.pdf")
    plt.close()

    ### multisite plot ###
    fig = plt.figure()
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(multisite_bins, scaled_multisite[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[3, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10, 2], linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, np.sum(scaled_multisite, axis = 0).tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth =2)

    plt.step(multisite_bins, scaled_multisite[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    plt.step(multisite_bins, scaled_multisite[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(multisite_bins, scaled_multisite[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    plt.step(multisite_bins, scaled_multisite[3, :].tolist() + [0], color = col5, where = 'post', dashes = [5, 2], linewidth = 2)
    plt.step(multisite_bins, scaled_multisite[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10, 2], linewidth = 2)
    plt.step(multisite_bins, np.sum(scaled_multisite, axis = 0).tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth =2)

    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Relative Counts", fontproperties = prop_font, fontsize = 26)
    # plt.title("Asimov Dataset", fontproperties= prop_font, fontsize = 26)
    plt.xlim((-1.355, -1.335))
    # plt.xlim((-1.29, -1.265))
    plt.ylim((0, 0.85))
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], linestyle = "dashdot", linewidth = 2, color = col2)[0],  \
                          plt.plot([], [], color = col5, dashes = [5, 2], linewidth = 2)[0], plt.plot([], [], linewidth = 2, dashes = [10, 2], color = col4)[0], plt.plot([], [], linewidth = 2, linestyle = "dotted", color = col3)[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "BiPo212", "BiPo214", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)

    ax = plt.gca()
    # ax.text(0.05, 0.55, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"3.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    ax.text(0.06, 0.25, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/asimov_multisite_searborn_{map}_coloured_2p5_5p0.pdf")
    plt.close()

    ### combined subplot style ###
    # fig, axes = plt.subplots(nrows = 1, ncols = 2)

    # axes[0].step(energy_bins, scaled_energy[0, :].tolist() + [0], color = "black", where = 'post', linewidth = 2)
    # axes[0].step(energy_bins, scaled_energy[1, :].tolist() + [0], color = "red", where = 'post', linewidth = 2)
    # axes[0].step(energy_bins, scaled_energy[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    # axes[0].step(energy_bins, scaled_energy[3, :].tolist() + [0], color = "orange", where = 'post', linewidth = 2)
    # axes[0].step(energy_bins, scaled_energy[4, :].tolist() + [0], color = "blue", where = 'post', linewidth = 2)

    # axes[0].set_xlabel("Reconstructed Energy [MeV]", fontproperties = prop_font, fontsize = 13)
    # axes[0].set_ylabel("Relative Counts", fontproperties = prop_font, fontsize = 13)
    # # axes[0].set_title("Asimov Dataset", fontproperties= prop_font, fontsize = 13)
    # axes[0].set_xlim((2.5, 5.0))
    # axes[0].set_ylim((0, 0.4))
    
    # for label in axes[0].get_xticklabels():
    #     label.set_fontproperties(prop_font)

    # for label in axes[0].get_yticklabels():
    #     label.set_fontproperties(prop_font)
    # axes[0].tick_params(labelsize = 12)

    # axes[0].legend(handles = [plt.plot([], [], color = "black")[0], plt.plot([], [], color = "red")[0],  \
    #                       plt.plot([], [], color = "orange")[0],\
    #                       plt.plot([], [], color = "blue")[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "BiPo212", "BiPo214"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 10),frameon=False)

    # axes[0].text(0.05, 0.83, "SNO+ Preliminary\n\nMC Only\n\nFV 4 m", transform=axes[0].transAxes, fontproperties = prop_font, fontsize = 10)
    
    # axes[1].step(multisite_bins, scaled_multisite[0, :].tolist() + [0], color = "black", where = 'post', linewidth = 2)
    # axes[1].step(multisite_bins, scaled_multisite[1, :].tolist() + [0], color = "red", where = 'post', linewidth = 2)
    # axes[1].step(multisite_bins, scaled_multisite[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    # axes[1].step(multisite_bins, scaled_multisite[3, :].tolist() + [0], color = "orange", where = 'post', linewidth = 2)
    # axes[1].step(multisite_bins, scaled_multisite[4, :].tolist() + [0], color = "blue", where = 'post', linewidth = 2)

    # axes[1].set_xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 13)
    # axes[1].set_ylabel("Relative Counts", fontproperties = prop_font, fontsize = 13)
    # # axes[1].set_title("Asimov Dataset", fontproperties= prop_font, fontsize = 13)
    # axes[1].set_xlim((-1.355, -1.335))
    # axes[1].set_ylim((0, 0.8))

    # axes[1].xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    # for label in axes[1].get_xticklabels():
    #     label.set_fontproperties(prop_font)

    # for label in axes[1].get_yticklabels():
    #     label.set_fontproperties(prop_font)
    # axes[1].tick_params(labelsize = 12)
    # axes[1].legend(handles = [plt.plot([], [], color = "black")[0], plt.plot([], [], color = "red")[0],  \
    #                       plt.plot([], [], color = "orange")[0],\
    #                       plt.plot([], [], color = "blue")[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "BiPo212", "BiPo214"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 10),frameon=False)

    # axes[1].text(0.05, 0.83, "SNO+ Preliminary\n\nMC Only\n\nFV 4 m", transform=axes[1].transAxes, fontproperties = prop_font, fontsize = 10)

    # fig.tight_layout()
    # fig.subplots_adjust(wspace=0.30)
    # plt.savefig("../plots/plots_for_approval/asimov_combined.pdf")
    # plt.close()

def calculate_uncertainty(profile_ll):
        """
        Function returns the 1 sigma confidence intervals on the minimum of profile LL
        fit.
        """
        
        # find idx of point on every log-likelihood curve closes in value to 1
        distance_to_interval = np.abs(profile_ll - 1) # absolute value of the difference
        closest_to_interval  = np.argmin(distance_to_interval, axis = 1) # find first intercept with 1 sigma level

        # create a masked array and remove the first intercept
        masked_1sig_full = np.ma.masked_array(distance_to_interval, mask = False)

        # mask out the closest interval in each row
        masked_1sig_full.mask[0, closest_to_interval[0]] = True
        masked_1sig_full.mask[1, closest_to_interval[1]] = True
        masked_1sig_full.mask[2, closest_to_interval[2]] = True

        # find second intercept
        second_closest = np.argmin(masked_1sig_full, axis = 1)

        # find minimum LL values as variables to save typing
        combined_min   = signal_hypothesis[np.argmin(profile_ll[0,:])]
        combined_upper = abs(combined_min - signal_hypothesis[second_closest[0]]) 
        combined_lower = abs(combined_min - signal_hypothesis[closest_to_interval[0]])
        multi_min      = signal_hypothesis[np.argmin(profile_ll[1,:])]
        multi_upper    = abs(multi_min - signal_hypothesis[second_closest[1]]) 
        multi_lower    = abs(multi_min - signal_hypothesis[closest_to_interval[1]])
        energy_min     = signal_hypothesis[np.argmin(profile_ll[2,:])]
        energy_upper   = abs(energy_min - signal_hypothesis[second_closest[2]]) 
        energy_lower   = abs(energy_min - signal_hypothesis[closest_to_interval[2]])

        return [combined_min, combined_upper, combined_lower], [multi_min, multi_upper, multi_lower], [energy_min, energy_upper, energy_lower]

def create_profile_LL_asimov_plot():
    """
    Create a plot of the result of the profile LL for the energy, multisite and
    combined fit based on the Asimov results.
    """
    
    # load the profile LL curves
    profile_ll = np.load("profile_ll_asimov_3p5_5p0.npy")
    # profile_ll = np.load("profile_ll_asimov.npy")
    combined_error, multisite_error, energy_error = calculate_uncertainty(profile_ll)
    signal_hypothesis = np.arange(0, 800, 1)
    
    ## create plot ##
    # cmap = plt.get_cmap("inferno")
    # col1 = cmap(0.1)
    # col2 = cmap(0.4)
    # col3 = cmap(0.7)
    # col4 = cmap(0.0)

    map = "colorblind"
    cmap = sns.color_palette(map, 8)
    col1 = cmap[1]
    col2 = cmap[2]
    col3 = cmap[7]
    col4 = cmap[0]
    # col5 = cmap[4]
    plt.figure()
    plt.plot(signal_hypothesis, profile_ll[2, :], linewidth = 2.5, linestyle = "--", color = col1)
    plt.plot(signal_hypothesis, profile_ll[1, :], linewidth = 2.5, color = col2)
    plt.plot(signal_hypothesis, profile_ll[0, :], linewidth = 2.5, linestyle = "dashdot", color = col3)
    
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    plt.xlabel(r"$^8$B Hypothesis", fontproperties = prop_font, fontsize = 26)
    plt.ylabel(r"$-2log(\mathcal{L})$", fontproperties = prop_font, fontsize = 26)

    plt.xlim((0, 80))
    # plt.xlim((0, 150))
    plt.ylim((0, 10))
    plt.axhline(1, color = col4, linestyle = "dotted", linewidth = 2.5, alpha = 1)
    # plt.axvline(66.3, color = "orange", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)

    plt.text(0.32, 0.28, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)
    # plt.legend(handles = [plt.plot([], [], color = "black", linewidth = 2.5)[0], plt.plot([], [], color = "green", linewidth = 2.5)[0], plt.plot([], [], color = "blue", linewidth = 2.5)[0], plt.plot([], [], color = "red", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)[0], plt.plot([], [], color = "orange", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)[0]], labels = [rf"Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$", rf"Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$", rf"Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$", r"1 $\sigma$ Frequentist", "True Signal Counts"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    plt.legend(handles = [plt.plot([], [], color = col1, linestyle = "--", linewidth = 2.5)[0], plt.plot([], [], color = col2, linewidth = 2.5)[0], plt.plot([], [], color = col3, linestyle = "dashdot", linewidth = 2.5)[0], plt.plot([], [], color = col4, linestyle = "dotted", linewidth = 2.5, alpha = 1)[0]], labels = [rf"Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$", rf"Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$", rf"Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$", r"1 $\sigma$ Frequentist"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)

    plt.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/profile_ll_asimov_{map}_coloured_3p5_5p0.pdf")
    plt.close()

def create_data_vs_mc_plot():
    """
    Function loads the data and MC dlog(L) information for B8 events above 5 MeV.
    Compares the dlog(L) distributions.
    """

    # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    data_solar = ROOT.TFile.Open(working_dir + "/extracted_data/above_5MeV/solar_data_discriminants/total.root")
    mc_solar   = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/B8_solar_nue/total.root")
    data_bi214 = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/discriminants_7.0.8_2.5_3.125_PDF/total.root")
    # data_bi214 = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.0_PDF/total.root")
    mc_bi214   = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/Bi214/total.root")

    # grab the above 5MeV TTrees
    solar_data_tree = data_solar.Get("5p0_plus")
    solar_mc_tree   = mc_solar.Get("5p0_plus")
    # bi214_data_tree = data_bi214.Get("1p25_3p0")
    # bi214_data_tree = data_bi214.Get("1p25_3p0")
    bi214_data_tree = data_bi214.Get("bi214")
    bi214_mc_tree   = mc_bi214.Get("1p25_3p0")
    
    FV_CUT            = 4500
    ITR_LOW           = 0.0 # not applied
    ITR_HIGH          = 1.0  # not applied
    BI214_E_CUT_LOW   = 1.25
    BI214_E_CUT_HIGH  = 3.0
    SOLAR_E_CUT_LOW   = 6.0
    SOLAR_E_CUT_HIGH  = 15.0
    solar_data_dlogL  = []
    solar_mc_dlogL    = []
    bi214_data_dlogL  = []
    bi214_mc_dlogL    = []
    
    for ientry in solar_data_tree:
        if ientry.energy >= SOLAR_E_CUT_LOW and ientry.energy < SOLAR_E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            solar_data_dlogL.append(ientry.itr)
    for ientry in solar_mc_tree:
        if ientry.energy >= SOLAR_E_CUT_LOW and ientry.energy < SOLAR_E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            solar_mc_dlogL.append(ientry.itr)

    for ientry in bi214_data_tree:
        if ientry.energy >= BI214_E_CUT_LOW and ientry.energy < BI214_E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
                continue
            bi214_data_dlogL.append(ientry.dlogL)
    for ientry in bi214_mc_tree:
        if ientry.energy >= BI214_E_CUT_LOW and ientry.energy < BI214_E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
                continue
            bi214_mc_dlogL.append(ientry.dlogL)       
    print("Num data events: ", len(bi214_data_dlogL))
    #### PLOTS ####
    # binning_solar = np.linspace(1.096, 1.12, 20)
    binning_solar = np.linspace(np.min(solar_data_dlogL), np.max(solar_data_dlogL), 25)
    solar_mids    = binning_solar[:-1] + np.diff(binning_solar)[0] / 2
    # binning_bi214 = np.linspace(-0.125, -0.11, 20)
    binning_bi214 = np.linspace(-0.125, -0.11, 20)
    # binning_bi214 = np.arange(0.15, 0.30, 0.005)
    # binning_bi214 = np.arange(-1.3531751965639947, -1.2049409232845083, 0.0005)
    print(np.min(bi214_mc_dlogL), np.max(bi214_mc_dlogL), np.min(bi214_data_dlogL), np.max(bi214_data_dlogL))
    # binning_bi214 = np.linspace(np.min(bi214_data_dlogL), np.max(bi214_data_dlogL), 25)
    bi214_mids    = binning_bi214[:-1] + np.diff(binning_bi214)[0] / 2

    # normalise the counts
    solar_counts_data, _ = np.histogram(solar_data_dlogL, bins = binning_solar)
    solar_counts_mc, _   = np.histogram(solar_mc_dlogL,  bins = binning_solar)
    solar_scale          = np.sum(solar_counts_data) / np.sum(solar_counts_mc)
    solar_mc_err         = np.sqrt(solar_counts_mc) * solar_scale
    solar_data_err       = np.sqrt(solar_counts_data)

    bi214_counts_data, _ = np.histogram(bi214_data_dlogL, bins = binning_bi214)
    bi214_counts_mc, _   = np.histogram(bi214_mc_dlogL,  bins = binning_bi214)
    bi214_scale          = np.sum(bi214_counts_data) / np.sum(bi214_counts_mc)
    print(bi214_scale)
    bi214_mc_err         = np.sqrt(bi214_counts_mc) * bi214_scale
    bi214_data_err       = np.sqrt(bi214_counts_data)

    # individual solar plot #
    fig = plt.figure()
    
    solar_counts_mc, _, _  = plt.hist(solar_mc_dlogL, bins = binning_solar, weights = np.full_like(solar_mc_dlogL, solar_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    solar_err_data = plt.errorbar(solar_mids, solar_counts_data, yerr = solar_data_err, linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
    # plt.xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.xlabel("Multisite", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Counts", fontproperties = prop_font, fontsize = 26)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(handles = [solar_err_data, plt.plot([], [], linewidth = 2, alpha = 0.7, color = "red")[0]], labels = ["Data", "MC"],  fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    plt.text(0.05, 0.68, "SNO+ Preliminary\n\nFV 4.5 m\n\nSolar Candidates\n\nE " + r"$\geq$ 6 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 15)
    fig.tight_layout()
    plt.savefig("../plots/plots_for_approval/solar_data_vs_mc.pdf")
    plt.close()

    # individual bi214 plot #
    fig = plt.figure()
    
    bi214_counts_mc, _, _  = plt.hist(bi214_mc_dlogL, bins = binning_bi214, weights = np.full_like(bi214_mc_dlogL, bi214_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # bi214_counts_mc, _, _  = plt.hist(bi214_mc_dlogL, bins = binning_bi214, density = True, alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    bi214_err_data = plt.errorbar(bi214_mids, bi214_counts_data, yerr = bi214_data_err, linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
    # plt.xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Counts", fontproperties = prop_font, fontsize = 26)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(handles = [solar_err_data, plt.plot([], [], linewidth = 2, alpha = 0.7, color = "red")[0]], labels = ["Data", "MC"],  fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    plt.text(0.05, 0.58, "SNO+ Preliminary\n\nFV 4.5 m\n\nBi214 Candidates\n\n" + rf"{BI214_E_CUT_LOW}$\leq$" " E " + rf"$\leq$ {BI214_E_CUT_HIGH} MeV" + "\n\n" + rf"{ITR_LOW} $\leq$ ITR $\leq$ {ITR_HIGH}", transform=ax.transAxes, fontproperties = prop_font, fontsize = 15)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/bi214_data_vs_mc_reprocessed_{BI214_E_CUT_LOW}_{BI214_E_CUT_HIGH}.pdf")
    plt.close()

    # combined plot #
    # fig, axes = plt.subplots(nrows = 1, ncols = 2)

    # solar_counts_mc, _, _  = axes[0].hist(solar_mc_dlogL, bins = binning_solar, weights = np.full_like(solar_mc_dlogL, solar_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    # solar_err_data = axes[0].errorbar(solar_mids, solar_counts_data, yerr = solar_data_err, linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
    # axes[0].set_xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 13)
    # axes[0].set_ylabel("Counts", fontproperties = prop_font, fontsize = 13)

    # for label in axes[0].get_xticklabels():
    #     label.set_fontproperties(prop_font)

    # for label in axes[0].get_yticklabels():
    #     label.set_fontproperties(prop_font)
    # axes[0].tick_params(labelsize = 12)

    # axes[0].legend(handles = [solar_err_data, plt.plot([], [], linewidth = 2, alpha = 0.7, color = "red")[0]], labels = ["Data", "MC"],  fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 10),frameon=False)
    # axes[0].text(0.04, 0.78, "SNO+ Preliminary\n\nFV 5.25 m\n\nSolar Candidates\n\nE " + r"$\geq$ 6 MeV", transform=axes[0].transAxes, fontproperties = prop_font, fontsize = 10)

    
    # bi214_counts_mc, _, _  = axes[1].hist(bi214_mc_dlogL, bins = binning_bi214, weights = np.full_like(bi214_mc_dlogL, bi214_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    # bi214_err_data = axes[1].errorbar(bi214_mids, bi214_counts_data, yerr = bi214_data_err, linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
    # axes[1].set_xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 13)
    # axes[1].set_ylabel("Counts", fontproperties = prop_font, fontsize = 13)

    # for label in axes[1].get_xticklabels():
    #     label.set_fontproperties(prop_font)

    # for label in axes[1].get_yticklabels():
    #     label.set_fontproperties(prop_font)
    # axes[1].tick_params(labelsize = 12)
    # axes[1].legend(handles = [solar_err_data, plt.plot([], [], linewidth = 2, alpha = 0.7, color = "red")[0]], labels = ["Data", "MC"],  fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 10),frameon=False)
    # axes[1].text(0.04, 0.78, "SNO+ Preliminary\n\nFV 5.25 m\n\nBi214 Candidates\n\n" + r"1.25$\leq$" " E " + r"$\leq$ 3.00 MeV", transform=axes[1].transAxes, fontproperties = prop_font, fontsize = 10)
    
    # fig.tight_layout()
    # fig.subplots_adjust(wspace = 0.3)
    # plt.savefig("../plots/plots_for_approval/bi214_solar_data_vs_mc.pdf")
    # plt.close()
    
def impact_of_reprocessing():
    """
    Create plots showing the difference in reconstruction following reprocessing
    from 7.0.8 --> 7.0.15. The 7.0.8 events have their multisite and recon
    variables plotted, without applying any energy, itr, or position cuts.

    This allows us to directly compare the changes after reprocessing.

    This will show changes in the data in terms of:

    recon X
    recon Y
    recon Z
    recon R
    recon E
    ITR
    Multisite
    """

     # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    data_bi214_708  = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/discriminants_7.0.8_2.5_3.125_PDF/total.root")
    data_bi214_7015 = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/total.root")

    bi214_data_708_tree  = data_bi214_708.Get("bi214")
    bi214_data_7015_tree = data_bi214_7015.Get("bi214")

    x_7015      = []
    x_708       = []
    y_7015      = []
    y_708       = []
    z_7015      = []
    z_708       = []
    r_7015      = []
    r_708       = []

    energy_708  = []
    energy_7015 = []

    itr_7015    = []
    itr_708     = []

    multi_7015  = []
    multi_708   = []


    for ientry in bi214_data_708_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        r = np.sqrt(x**2 + y **2 + z**2)
        
        itr = ientry.itr

        multisite = ientry.dlogL

        energy = ientry.energy

        x_708.append(x)
        y_708.append(y) 
        z_708.append(z) 
        r_708.append(r) 
        itr_708.append(itr) 
        multi_708.append(multisite) 
        energy_708.append(energy)
    for ientry in bi214_data_7015_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        r = np.sqrt(x**2 + y **2 + z**2)
        
        itr = ientry.itr

        multisite = ientry.dlogL

        energy = ientry.energy

        x_7015.append(x)
        y_7015.append(y) 
        z_7015.append(z) 
        r_7015.append(r) 
        itr_7015.append(itr) 
        multi_7015.append(multisite) 
        energy_7015.append(energy)
    
    x_7015      = np.array(x_7015)
    y_7015      = np.array(y_7015)
    z_7015      = np.array(z_7015)
    r_7015      = np.array(r_7015)
    energy_7015 = np.array(energy_7015)
    itr_7015    = np.array(itr_7015)
    multi_7015  = np.array(multi_7015)
    x_708       = np.array(x_708)
    y_708       = np.array(y_708)
    z_708       = np.array(z_708)
    r_708       = np.array(r_708)
    energy_708  = np.array(energy_708)
    itr_708     = np.array(itr_708)
    multi_708   = np.array(multi_708)

    #### position reconstruction differences ####
    plt.figure()
    plt.hist(x_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(x_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Reconstructed X (mm)")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_x_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(y_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(y_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Reconstructed Y (mm)")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_y_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(z_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(z_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Reconstructed Z (mm)")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_z_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(r_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(r_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Reconstructed R (mm)")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_r_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(energy_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(energy_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Reconstructed Energy (MeV)")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_energy_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(itr_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(itr_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("ITR")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_itr_recon_diffs.pdf")
    plt.close()

    plt.figure()
    plt.hist(multi_708, bins = 50, density = True, histtype = "step", label = "7.0.8")
    plt.hist(multi_7015, bins = 50, density = True, histtype = "step", label = "7.0.15")
    plt.xlabel("Multisite")
    plt.legend()
    plt.savefig("../plots/plots_for_approval/reprocessing/bi214_multi_recon_diffs.pdf")
    plt.close()

signal_hypothesis = np.arange(0, 800, 1)
plt.rcParams['xtick.major.pad'] = '6' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '6'
# create_time_residual_pdfs()
# create_dlogL_plot()
plt.rcParams['xtick.major.pad'] = '12' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '12'
# create_asimov_dataset_graphs()
plt.rcParams['xtick.major.pad'] = '6' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '6'
# create_profile_LL_asimov_plot()
# create_data_vs_mc_plot()
impact_of_reprocessing()