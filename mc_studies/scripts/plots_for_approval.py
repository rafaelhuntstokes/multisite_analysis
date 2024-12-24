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
from scipy.optimize import curve_fit

"""
USE env_poster in nuPhysPoster!
"""
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
    ax = fig.add_axes([0.13, 0.13, 0.8, 0.8])
    b8 = ax.step(multi_B8_bins[:-1], multi_B8_counts,  where = 'post', linewidth = 2, color = col1, label = r"$^8B$")
    tl208 = ax.step(multi_B8_bins[:-1], multi_Tl208_counts,  where = 'post', linewidth = 2, linestyle = "dashed", alpha = 1, color = col2, label = r"$^{208}Tl$")
    plt.legend(handles = [plt.plot([], [], color = col1)[0], plt.plot([], [], linestyle = "dashed", color = col2)[0]], labels = [r"${^8}$B $\nu _e$ - e$^-$ES", r"$^{208}$Tl - $\beta \gamma$"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 18),frameon=False)
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
    ax = fig.add_axes([0.13, 0.13, 0.8, 0.8])
    # plt.step(binning_3p5_5p0, b8_pdf.tolist() + [0], where = "post", linewidth = 2, color = col1)
    # plt.step(binning_3p5_5p0, tl208_pdf.tolist() + [0], where = "post", linestyle = "dashed", linewidth = 2, color = col2)
    ax.step(binning_full, b8_pdf.tolist() + [0], where = "post", linewidth = 2, color = col1)
    ax.step(binning_full, tl208_pdf.tolist() + [0], where = "post", linestyle = "dashed", linewidth = 2, color = col2)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Normalised Counts", fontproperties = prop_font, fontsize = 26)
    # plt.title(r"Multisite Discriminants for $^8$B and $^{208}$Tl", fontproperties = prop_font, fontsize = 26)
    plt.ylim((0, 0.139))

    plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], color = col2, linewidth = 2, linestyle = "dashed")[0]], labels = [r"${^8}$B $\nu _e$ - e$^-$ES", r"$^{208}$Tl - $\beta \gamma$"], loc = "upper right", fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 18),frameon=False)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    ax.text(0.06, 0.6, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)

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
    left, bottom, width, height = [0.23, 0.26, 0.2*factor, 0.2*factor]
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
    col3 = "black"#cmap[7]
    col4 = cmap[3]
    col5 = cmap[4]

    # load the pdfs in terms of energy and multisite discriminant
    # pdfs_multisite = np.load("multisite_pdf_array_3p5_5p0.npy")
    # pdfs_energy    = np.load("energy_pdf_array_3p5_5p0.npy")
    pdfs_multisite = np.load("multisite_pdf_array_2p5_5p0.npy")
    pdfs_energy    = np.load("energy_pdf_array_2p5_5p0.npy")

    # define the normalisations of each - from background model and livetime
    # normalisations = np.array([66.3, 468.9, 0.9081, 59.227, 38.528]) # 145.7 days livetime (what I have in dataset)
    normalisations = np.array([230.0, 1174.61, 2.27, 148.38, 96.52]) # 1 yr scaled livetime
    weights        = np.array([0.551, 0.634, 0.620, 0, 0])
    weights        = np.array([1, 1, 1, 1, 1])
    normalisations = normalisations * weights

    # define everything to the solar signal, which has norm 1
    normalisations = normalisations / normalisations[0]

    # multiply each PDF by the respective normalisations 
    scaled_energy    = normalisations[:, None] * pdfs_energy
    scaled_multisite = normalisations[:, None] * pdfs_multisite

    energy_bins            = np.arange(2.5, 5.05, 0.1)
    multisite_bins         = np.linspace(-1.355, -1.335, 50)#np.arange(-1.375, -1.325, 0.0005)
    multisite_bins_3p5_5p0 = np.arange(-1.309, -1.26, 0.0005)
    energy_bins_3p5_5p0    = np.arange(2.5, 5.05, 0.05)
    # create the individual plots
    
    ### energy plot ###
    fig = plt.figure()
    ax = fig.add_axes([0.13, 0.13, 0.8, 0.8])
    # plt.step(energy_bins_3p5_5p0, scaled_energy[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    # plt.step(energy_bins_3p5_5p0, scaled_energy[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # # plt.step(energy_bins, scaled_energy[2, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # # plt.step(energy_bins_3p5_5p0, scaled_energy[3, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # # plt.step(energy_bins_3p5_5p0, scaled_energy[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10,2], linewidth = 2)
    # plt.step(energy_bins_3p5_5p0, np.sum(scaled_energy, axis= 0).tolist() + [0], color = col3, linestyle = "dotted", where = 'post', linewidth = 2)

    ax.step(energy_bins, scaled_energy[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    ax.step(energy_bins, scaled_energy[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(energy_bins, scaled_energy[2, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    ax.step(energy_bins, scaled_energy[3, :].tolist() + [0], color = col5, where = 'post', dashes = [5, 2], linewidth = 2)
    ax.step(energy_bins, scaled_energy[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10,2], linewidth = 2)
    ax.step(energy_bins, np.sum(scaled_energy, axis= 0).tolist() + [0], color = col3, linestyle = "dotted", where = 'post', linewidth = 2)


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
                          plt.plot([], [], color = col4, linewidth = 2, dashes = [10, 2])[0], plt.plot([], [], color = col3, linestyle = "dotted", linewidth = 2)[0]], labels = [r"${^8}$B $\nu _e$ - e$^-$ES", r"$^{208}$Tl - $\beta \gamma$", "BiPo212 (pileup)", "$^{214}$Bi", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path2, size = 16),frameon=False)
    # plt.legend(handles = [plt.plot([], [], color = col1, linewidth = 2)[0], plt.plot([], [], color = col2, linewidth = 2, linestyle = "dashdot")[0],  \
    #                       plt.plot([], [], color = col3, linewidth = 2, linestyle = "dotted")[0]], labels = [r"${^8}$B", r"$^{208}$Tl", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)

    # ax.text(0.66, 0.42, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    ax.text(0.06, 0.72, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/asimov_energy_seaborn_{map}_coloured_2p5_5p0_2.pdf")
    plt.close()

    ### multisite plot ###
    fig = plt.figure()
    ax = fig.add_axes([0.13, 0.13, 0.8, 0.8])
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(multisite_bins, scaled_multisite[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[3, :].tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, scaled_multisite[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10, 2], linewidth = 2)
    # plt.step(multisite_bins_3p5_5p0, np.sum(scaled_multisite, axis = 0).tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth =2)

    ax.step(multisite_bins, scaled_multisite[0, :].tolist() + [0], color = col1, where = 'post', linewidth = 2)
    ax.step(multisite_bins, scaled_multisite[1, :].tolist() + [0], color = col2, where = 'post', linestyle = "dashdot", linewidth = 2)
    # plt.step(multisite_bins, scaled_multisite[2, :].tolist() + [0], color = "green", where = 'post', linewidth = 2)
    ax.step(multisite_bins, scaled_multisite[3, :].tolist() + [0], color = col5, where = 'post', dashes = [5, 2], linewidth = 2)
    ax.step(multisite_bins, scaled_multisite[4, :].tolist() + [0], color = col4, where = 'post', dashes = [10, 2], linewidth = 2)
    ax.step(multisite_bins, np.sum(scaled_multisite, axis = 0).tolist() + [0], color = col3, where = 'post', linestyle = "dotted", linewidth =2)

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
                          plt.plot([], [], color = col5, dashes = [5, 2], linewidth = 2)[0], plt.plot([], [], linewidth = 2, dashes = [10, 2], color = col4)[0], plt.plot([], [], linewidth = 2, linestyle = "dotted", color = col3)[0]], labels = [r"${^8}$B $\nu _e$ - e$^-$ES", r"$^{208}$Tl - $\beta \gamma$", "BiPo212 (pileup)", "$^{214}$Bi", "Total Model"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path2, size = 16),frameon=False)

    ax = plt.gca()
    # ax.text(0.05, 0.55, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"3.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 20)
    ax.text(0.06, 0.25, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 18)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/asimov_multisite_searborn_{map}_coloured_2p5_5p0_2.pdf")
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
    # plt.savefig("../plots/plots_for_approval/asimov_combined.pdf")bounds
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
    # profile_ll = np.load("profile_ll_asimov_3p5_5p0.npy")
    # profile_ll = np.load("profile_ll_asimov.npy")
    profile_ll = np.load('profile_ll_asimov_1yr_fullROI_2.npy')
    combined_error, multisite_error, energy_error = calculate_uncertainty(profile_ll)
    signal_hypothesis = np.arange(0, 400, 1)
    
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
    fig = plt.figure()
    ax = fig.add_axes([0.13, 0.13, 0.8, 0.8])
    ax.plot(signal_hypothesis, profile_ll[2, :], linewidth = 2.5, linestyle = "--", color = col1)
    ax.plot(signal_hypothesis, profile_ll[1, :], linewidth = 2.5, linestyle = "dashdot", color = col2)
    ax.plot(signal_hypothesis, profile_ll[0, :], linewidth = 2.5, color = "black")
    
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    plt.xlabel(r"Expected $^8$B $\nu _e$ - e$^-$ ES  Interactions / Year", fontproperties = prop_font, fontsize = 26)
    plt.ylabel(r"$-2log(\mathcal{L})$", fontproperties = prop_font, fontsize = 26)

    # plt.xlim((0, 80))
    plt.xlim((0, 460))
    plt.ylim((0, 10))
    plt.axhline(1, color = "black", linestyle = "dotted", linewidth = 2, alpha = 1)
    # plt.axvline(66.3, color = "orange", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)

    plt.text(0.03, 0.15, "SNO+ Preliminary\n\nMC Only\n\nFV 4.5 m\n\n" + r"2.5 $\leq$ E $\leq$ 5.0 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 15)
    # plt.legend(handles = [plt.plot([], [], color = "black", linewidth = 2.5)[0], plt.plot([], [], color = "green", linewidth = 2.5)[0], plt.plot([], [], color = "blue", linewidth = 2.5)[0], plt.plot([], [], color = "red", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)[0], plt.plot([], [], color = "orange", linestyle = "dashed", linewidth = 2.5, alpha = 0.5)[0]], labels = [rf"Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$", rf"Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$", rf"Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$", r"1 $\sigma$ Frequentist", "True Signal Counts"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    plt.legend(handles = [plt.plot([], [], color = col1, linestyle = "--", linewidth = 2.5)[0], plt.plot([], [], color = col2, linestyle = "dashdot", linewidth = 2.5)[0], plt.plot([], [], color = "black", linewidth = 2.5)[0], plt.plot([], [], color = "black", linestyle = "dotted", linewidth = 2, alpha = 1)[0]], labels = [rf"Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$", rf"Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$", rf"Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$", r"1 $\sigma$ Frequentist"], fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path2, size = 16),frameon=False)

    plt.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/profile_ll_asimov_{map}_coloured_2.pdf")
    plt.close()

def create_data_vs_mc_plot():
    """
    Function loads the data and MC dlog(L) information for B8 events above 5 MeV.
    Compares the dlog(L) distributions.
    """

    # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    data_solar = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/reprocessed_solar_multisite_7.0.15_highE_PDF/total.root")
    mc_solar   = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/solar_mc_multisite_highE_PDF/total.root")
    data_bi214 = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/reprocessed_bi214_multisite_7.0.15/total.root")
    # data_bi214 = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.0_PDF/total.root")
    mc_bi214   = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/bi214_mc_multisite/total.root")

    # grab the above 5MeV TTrees
    solar_data_tree = data_solar.Get("above_5MeV")
    solar_mc_tree   = mc_solar.Get("above_5MeV")
    # bi214_data_tree = data_bi214.Get("1p25_3p0")
    # bi214_data_tree = data_bi214.Get("1p25_3p0")
    bi214_data_tree = data_bi214.Get("bi214")
    bi214_mc_tree   = mc_bi214.Get("bi214")
    
    FV_CUT            = 4500
    ITR_LOW           = 0.22
    ITR_HIGH          = 0.3 
    BI214_E_CUT_LOW   = 1.25
    BI214_E_CUT_HIGH  = 3.0
    SOLAR_E_CUT_LOW   = 6.0
    SOLAR_E_CUT_HIGH  = 15.0
    SOLAR_ITR_LOW     = 0.21
    SOLAR_ITR_HIGH    = 0.3
    SOLAR_HC_RATIO_LOW = 0.7
    SOLAR_FV_CUT       = 4500
    
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
            if r > SOLAR_FV_CUT:
                continue
            if ientry.itr < SOLAR_ITR_LOW or ientry.itr > SOLAR_ITR_HIGH:
                continue
            if ientry.dlogL < 1.105:
                print(f"dlogL: {ientry.dlogL}\nGTID: {ientry.gtid}")
            if ientry.HC_ratio < SOLAR_HC_RATIO_LOW:
                continue
            solar_data_dlogL.append(ientry.dlogL)
    print("Num solar events: ", len(solar_data_dlogL))
    for ientry in solar_mc_tree:
        if ientry.energy >= SOLAR_E_CUT_LOW and ientry.energy < SOLAR_E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > SOLAR_FV_CUT:
                continue
            if ientry.itr < SOLAR_ITR_LOW or ientry.itr > SOLAR_ITR_HIGH:
                continue

            if ientry.HC_ratio < SOLAR_HC_RATIO_LOW:
                continue
            solar_mc_dlogL.append(ientry.dlogL)

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
    print("Num Bi214 events: ", len(bi214_data_dlogL))
    #### PLOTS ####
    binning_solar = np.linspace(1.096, 1.12, 20)
    # binning_solar = np.linspace(-0.121,-0.111, 12)
    solar_mids    = binning_solar[:-1] + np.diff(binning_solar)[0] / 2
    # binning_bi214 = np.linspace(-0.125, -0.11, 20)
    binning_bi214 = np.linspace(-0.125, -0.11, 20)
    # binning_bi214 = np.arange(0.15, 0.30, 0.005)
    # binning_bi214 = np.arange(-1.3531751965639947, -1.2049409232845083, 0.0005)
    print(np.min(solar_mc_dlogL), np.max(solar_mc_dlogL), np.min(solar_data_dlogL), np.max(solar_data_dlogL))
    # binning_bi214 = np.linspace(np.min(bi214_data_dlogL), np.max(bi214_data_dlogL), 25)
    bi214_mids    = binning_bi214[:-1] + np.diff(binning_bi214)[0] / 2

    # normalise the counts
    solar_counts_data, _ = np.histogram(solar_data_dlogL, bins = binning_solar)
    # identify zero data counts to remove
    id_nonzero_solar = np.nonzero(solar_counts_data)
    solar_counts_mc, _   = np.histogram(solar_mc_dlogL,  bins = binning_solar)
    solar_scale          = np.sum(solar_counts_data) / np.sum(solar_counts_mc)
    solar_mc_err         = np.sqrt(solar_counts_mc) * solar_scale
    solar_data_err       = np.sqrt(solar_counts_data)

    bi214_counts_data, _ = np.histogram(bi214_data_dlogL, bins = binning_bi214)
    id_nonzero_bi = np.nonzero(bi214_counts_data)
    bi214_counts_mc, _   = np.histogram(bi214_mc_dlogL,  bins = binning_bi214)
    bi214_scale          = np.sum(bi214_counts_data) / np.sum(bi214_counts_mc)
    print(bi214_scale)
    bi214_mc_err         = np.sqrt(bi214_counts_mc) * bi214_scale
    bi214_data_err       = np.sqrt(bi214_counts_data)

    # individual solar plot #
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.12, 0.8, 0.8])
    solar_counts_mc, _, _  = ax.hist(solar_mc_dlogL, bins = binning_solar, weights = np.full_like(solar_mc_dlogL, solar_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    solar_err_data = ax.errorbar(solar_mids[id_nonzero_solar], solar_counts_data[id_nonzero_solar], yerr = solar_data_err[id_nonzero_solar], linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
    # plt.xlabel("Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.xlabel(r"$\Delta log(\mathcal{L})$ Multisite Discriminant", fontproperties = prop_font, fontsize = 26)
    plt.ylabel("Counts", fontproperties = prop_font, fontsize = 26)

    plt.yticks(np.arange(0, 41, 2), fontsize = 20)
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.xticks(fontsize = 20)
    
    plt.legend(handles = [solar_err_data, plt.plot([], [], linewidth = 2, alpha = 0.7, color = "red")[0]], labels = ["Data", "MC"],  fancybox=False, numpoints=1,prop=fm.FontProperties(fname=path, size = 16),frameon=False)
    # plt.text(0.05, 0.48, f"SNO+ Preliminary\n\n{np.sum(solar_counts_data)} Solar Candidates\n\nFV 4.5 m\n\nE " + r"$\geq$ 6 MeV" + "\n\n" + rf"{SOLAR_ITR_LOW} $\leq$ ITR $\leq$ {SOLAR_ITR_HIGH}" + "\n\nHC Ratio " + rf"< {SOLAR_HC_RATIO_LOW}", transform=ax.transAxes, fontproperties = prop_font, fontsize = 16)
    plt.text(0.05, 0.48, f"SNO+ Preliminary\n\nSolar Candidates ({np.sum(solar_counts_data)})\n\nFV 4.5 m\n\nE " + r"$\geq$ 6 MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 16)
    fig.tight_layout()
    plt.savefig("../plots/plots_for_approval/solar_data_vs_mc_NEW_highE.pdf")
    plt.close()

    # individual bi214 plot #
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.12, 0.8, 0.8])
    bi214_counts_mc, _, _  = ax.hist(bi214_mc_dlogL, bins = binning_bi214, weights = np.full_like(bi214_mc_dlogL, bi214_scale), alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # bi214_counts_mc, _, _  = plt.hist(bi214_mc_dlogL, bins = binning_bi214, density = True, alpha = 0.7, linewidth = 2, histtype = "step", color = "red")
    # solar_err_mc   = plt.errorbar(solar_mids, solar_counts_mc, yerr = solar_mc_err, linewidth = 2, color = "red", capsize = 1.5)
    bi214_err_data = ax.errorbar(bi214_mids[id_nonzero_bi], bi214_counts_data[id_nonzero_bi], yerr = bi214_data_err[id_nonzero_bi], linestyle = "", marker = "o", color = "black", capsize = 1.5)
    
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
    plt.text(0.05, 0.52, f"SNO+ Preliminary\n\n" + r"$^{214}Bi$ Candidates " + f"({np.sum(bi214_counts_data)})\n\nFV {FV_CUT / 1000:.1f} m\n\n" + rf"{BI214_E_CUT_LOW}$\leq$" " E " + rf"$\leq$ {BI214_E_CUT_HIGH} MeV", transform=ax.transAxes, fontproperties = prop_font, fontsize = 16)
    fig.tight_layout()
    plt.savefig(f"../plots/plots_for_approval/bi214_data_vs_mc_reprocessed_{BI214_E_CUT_LOW}_{BI214_E_CUT_HIGH}_NEW.pdf")
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

def data_quality():
    """
    Make some 2D histos to find best cut values to apply to Bi214 and Solar
    data.
    """

    # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    data_solar = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/reprocessed_solar_multisite_7.0.15_highE_PDF/total.root")
    mc_solar   = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/solar_mc_multisite_highE_PDF/total.root")
    data_bi214 = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/reprocessed_bi214_multisite_7.0.15/total.root")
    mc_bi214   = ROOT.TFile.Open(working_dir + "/extracted_data/better_tagging/bi214_mc_multisite/total.root")

    # grab the above 5MeV TTrees
    solar_data_tree = data_solar.Get("above_5MeV")
    solar_mc_tree   = mc_solar.Get("above_5MeV")
    bi214_data_tree = data_bi214.Get("bi214")
    bi214_mc_tree   = mc_bi214.Get("bi214")

    # apply FV and energy cuts
    FV_CUT       = 4500
    SOLAR_E_LOW  = 6
    SOLAR_E_HIGH = 15
    BI_E_LOW     = 1.25
    BI_E_HIGH    = 3.00

    # keep track of DC variables
    solar_mc_HC       = []
    solar_mc_posFOM   = []
    solar_mc_itr      = []
    solar_data_HC     = []
    solar_data_posFOM = []
    solar_data_itr    = []

    bi214_mc_HC       = []
    bi214_mc_posFOM   = []
    bi214_mc_itr      = []
    bi214_data_HC     = []
    bi214_data_posFOM = []
    bi214_data_itr    = []

    ## solar ##
    for ientry in solar_data_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        e = ientry.energy

        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue
        if e > SOLAR_E_HIGH or e < SOLAR_E_LOW:
            continue

        hc_ratio = ientry.HC_ratio
        itr      = ientry.itr
        posFOM   = ientry.posFOM / ientry.posFOM_hits

        solar_data_HC.append(hc_ratio)
        solar_data_itr.append(itr)
        solar_data_posFOM.append(posFOM)
    for ientry in solar_mc_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        e = ientry.energy

        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue
        if e > SOLAR_E_HIGH or e < SOLAR_E_LOW:
            continue

        hc_ratio = ientry.HC_ratio
        itr      = ientry.itr
        posFOM   = ientry.posFOM / ientry.posFOM_hits

        solar_mc_HC.append(hc_ratio)
        solar_mc_itr.append(itr)
        solar_mc_posFOM.append(posFOM)

    ## make some comparison plots ## 
    fig, axes = plt.subplots(nrows = 3, ncols = 3)
    bins_itr = np.linspace(0.1, 0.3, 20)
    bins_HC  = np.linspace(0.5, 1, 40)
    bins_posFOM = np.linspace(10, 15, 40)
    axes[0,0].hist(solar_data_itr, bins = bins_itr, density = True, histtype = "step", label = "Data")
    axes[0,0].hist(solar_mc_itr, bins = bins_itr, density = True, histtype = "step", label = "MC")
    axes[0,0].set_xlabel("ITR")
    axes[0,0].legend(frameon = False)

    axes[0,1].hist(solar_data_HC, bins = bins_HC, density = True, histtype = "step", label = "Data")
    axes[0,1].hist(solar_mc_HC, bins = bins_HC, density = True, histtype = "step", label = "MC")
    axes[0,1].set_xlabel("HC Ratio")
    axes[0,1].legend(frameon = False)

    axes[0,2].hist(solar_data_posFOM, bins = bins_posFOM, density = True, histtype = "step", label = "Data")
    axes[0,2].hist(solar_mc_posFOM, bins = bins_posFOM, density = True, histtype = "step", label = "MC")
    axes[0,2].set_xlabel("posFOM")
    axes[0,2].legend(frameon = False)

    axes[1,0].hist2d(solar_data_itr, solar_data_HC, bins = [bins_itr, bins_HC], density = True, cmin = 1e-4)
    axes[1,0].set_xlabel("ITR")
    axes[1,0].set_ylabel("HC Ratio")
    axes[1,0].set_title("Data")
    
    axes[1,1].hist2d(solar_data_itr, solar_data_posFOM, bins = [bins_itr, bins_posFOM], density = True, cmin = 1e-4)
    axes[1,1].set_xlabel("ITR")
    axes[1,1].set_ylabel("posFOM")
    axes[1,1].set_title("Data")

    axes[1,2].hist2d(solar_data_posFOM, solar_data_HC, bins = [bins_posFOM, bins_HC], density = True, cmin = 1e-4)
    axes[1,2].set_xlabel("posFOM")
    axes[1,2].set_ylabel("HC Ratio")
    axes[1,2].set_title("Data")

    axes[2,0].hist2d(solar_mc_itr, solar_mc_HC, bins = [bins_itr, bins_HC], density = True, cmin = 1e-4)
    axes[2,0].set_xlabel("ITR")
    axes[2,0].set_ylabel("HC Ratio")
    axes[2,0].set_title("MC")
    
    axes[2,1].hist2d(solar_mc_itr, solar_mc_posFOM, bins = [bins_itr, bins_posFOM], density = True, cmin = 1e-4)
    axes[2,1].set_xlabel("ITR")
    axes[2,1].set_ylabel("posFOM")
    axes[2,1].set_title("MC")

    axes[2,2].hist2d(solar_mc_posFOM, solar_mc_HC, bins = [bins_posFOM, bins_HC], density = True, cmin = 1e-4)
    axes[2,2].set_xlabel("posFOM")
    axes[2,2].set_ylabel("HC Ratio")
    axes[2,2].set_title("MC")

    fig.tight_layout()
    plt.savefig("../plots/plots_for_approval/solar_data_quality.pdf")
    plt.close()

    ## Bi214 ##
    for ientry in bi214_data_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        e = ientry.energy

        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue
        if e > BI_E_HIGH or e < BI_E_LOW:
            continue

        hc_ratio = ientry.HC_ratio
        itr      = ientry.itr
        posFOM   = ientry.posFOM / ientry.posFOM_hits

        bi214_data_HC.append(hc_ratio)
        bi214_data_itr.append(itr)
        bi214_data_posFOM.append(posFOM)
    for ientry in bi214_mc_tree:
        x = ientry.x
        y = ientry.y
        z = ientry.z
        e = ientry.energy

        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue
        if e > BI_E_HIGH or e < BI_E_LOW:
            continue

        hc_ratio = ientry.HC_ratio
        itr      = ientry.itr
        posFOM   = ientry.posFOM / ientry.posFOM_hits

        bi214_mc_HC.append(hc_ratio)
        bi214_mc_itr.append(itr)
        bi214_mc_posFOM.append(posFOM)

    ## make some comparison plots ## 
    fig, axes = plt.subplots(nrows = 3, ncols = 3)
    bins_itr = np.linspace(0.1, 0.3, 20)
    bins_HC  = np.linspace(0.5, 1, 40)
    bins_posFOM = np.linspace(10, 15, 40)
    axes[0,0].hist(bi214_data_itr, bins = bins_itr, density = True, histtype = "step", label = "Data")
    axes[0,0].hist(bi214_mc_itr, bins = bins_itr, density = True, histtype = "step", label = "MC")
    axes[0,0].set_xlabel("ITR")
    axes[0,0].legend(frameon = False)

    axes[0,1].hist(bi214_data_HC, bins = bins_HC, density = True, histtype = "step", label = "Data")
    axes[0,1].hist(bi214_mc_HC, bins = bins_HC, density = True, histtype = "step", label = "MC")
    axes[0,1].set_xlabel("HC Ratio")
    axes[0,1].legend(frameon = False)

    axes[0,2].hist(bi214_data_posFOM, bins = bins_posFOM, density = True, histtype = "step", label = "Data")
    axes[0,2].hist(bi214_mc_posFOM, bins = bins_posFOM, density = True, histtype = "step", label = "MC")
    axes[0,2].set_xlabel("posFOM")
    axes[0,2].legend(frameon = False)

    axes[1,0].hist2d(bi214_data_itr, bi214_data_HC, bins = [bins_itr, bins_HC], density = True, cmin = 1e-4)
    axes[1,0].set_xlabel("ITR")
    axes[1,0].set_ylabel("HC Ratio")
    axes[1,0].set_title("Data")
    
    axes[1,1].hist2d(bi214_data_itr, bi214_data_posFOM, bins = [bins_itr, bins_posFOM], density = True, cmin = 1e-4)
    axes[1,1].set_xlabel("ITR")
    axes[1,1].set_ylabel("posFOM")
    axes[1,1].set_title("Data")

    axes[1,2].hist2d(bi214_data_posFOM, bi214_data_HC, bins = [bins_posFOM, bins_HC], density = True, cmin = 1e-4)
    axes[1,2].set_xlabel("posFOM")
    axes[1,2].set_ylabel("HC Ratio")
    axes[1,2].set_title("Data")

    axes[2,0].hist2d(bi214_mc_itr, bi214_mc_HC, bins = [bins_itr, bins_HC], density = True, cmin = 1e-4)
    axes[2,0].set_xlabel("ITR")
    axes[2,0].set_ylabel("HC Ratio")
    axes[2,0].set_title("MC")
    
    axes[2,1].hist2d(bi214_mc_itr, bi214_mc_posFOM, bins = [bins_itr, bins_posFOM], density = True, cmin = 1e-4)
    axes[2,1].set_xlabel("ITR")
    axes[2,1].set_ylabel("posFOM")
    axes[2,1].set_title("MC")

    axes[2,2].hist2d(bi214_mc_posFOM, bi214_mc_HC, bins = [bins_posFOM, bins_HC], density = True, cmin = 1e-4)
    axes[2,2].set_xlabel("posFOM")
    axes[2,2].set_ylabel("HC Ratio")
    axes[2,2].set_title("MC")

    fig.tight_layout()
    plt.savefig("../plots/plots_for_approval/bi214_data_quality.pdf")
    plt.close()

    plt.figure()
    plt.hist(bi214_mc_itr, bins = bins_itr, density = True, histtype = "step", color = "red", linestyle = "dotted")
    plt.hist(bi214_data_itr, bins = bins_itr, density = True, histtype = "step", color = "black")
    plt.legend(handles = [plt.plot([], [], color = "red", linestyle = "dotted")[0], plt.plot([], [], color = "black")[0]], labels = ["MC", "Data"], frameon = False, prop = fm.FontProperties(fname=path2, size = 16))
    plt.xlabel("ITR", fontproperties = prop_font, fontsize = 14)
    plt.ylabel("Counts", fontproperties = prop_font, fontsize = 14)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    plt.savefig("../plots/plots_for_approval/itr_mc_vs_data.pdf")
    plt.close()

def bipo_tagging_data_plots():
    """
    Function makes the normal energy, DR and fitted dT plots to check the quality
    of the BiPo tagging performed.
    """

    data_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/better_tagging/raw_tagging_output"

    tagging_file = ROOT.TFile.Open(f"{data_dir}/total.root")
    tree         = tagging_file.Get("tag_info")

    FV_CUT    = 4500
    energy_bi = []
    energy_po = []
    dT        = []
    dR        = []
    x_bi      = []
    y_bi      = []
    z_bi      = []

    for ientry in tree:
        x = ientry.x_bi
        y = ientry.y_bi
        z = ientry.z_bi

        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue

        dT.append(ientry.dT)
        dR.append(ientry.dR)
        x_bi.append(x)
        y_bi.append(y)
        z_bi.append(z)
        energy_bi.append(ientry.energy_bi)
        energy_po.append(ientry.energy_po)

    x_bi  = np.array(x_bi)
    y_bi  = np.array(y_bi)
    z_bi  = np.array(z_bi)

    rho2  = x_bi**2 + y_bi**2
    rho2  = rho2 / 6000**2

    # fit in mu seconds
    dT = np.array(dT)
    dT = dT / 1000

    binning_pos = np.arange(-6000, 6020, 50)
    binning_dr  = np.arange(0, 1050, 50)
    binning_e   = np.arange(0, 4.1, 0.05)
    # energy plot #
    plt.hist(energy_bi, bins = binning_e, histtype = "step", color = "C0")
    plt.hist(energy_po, bins = binning_e, histtype = "step", color = "C1")
    plt.legend(handles = [plt.plot([], [], color = "C0")[0], plt.plot([], [], color = "C0")[0]], labels = [r"$^{214}$Bi", r"$^{214}$Po"], frameon = False)
    plt.xlabel("Reconstructed Energy (MeV)")
    plt.ylabel("Counts")
    plt.savefig("../plots/plots_for_approval/bipo_tagging_energy.pdf")
    plt.close()

    # dR plot #
    plt.hist(dR, bins = binning_dr, histtype = "step", color = "C0")
    plt.legend(handles = [plt.plot([], [], color = "C0")[0]], labels = [r"Tagged BiPo214"], frameon = False)
    plt.xlabel(r"$\Delta r$ (mm)")
    plt.ylabel("Counts")
    plt.savefig("../plots/plots_for_approval/bipo_tagging_dR.pdf")
    plt.close()

    # dT #
    width = 5
    binning = np.arange(4, 1000+ width, width)
    counts, bins = np.histogram(dT, bins = binning)
    mids = bins[:-1] + np.diff(bins)[0]/2
    err = plt.errorbar(mids, counts, yerr = np.sqrt(counts), marker = "o", linestyle = "", capsize = 2, color = "black")
    plt.xlabel(r"$\Delta t$ $[\mu s]$", fontproperties = prop_font, loc = "right")
    plt.ylabel(f"Counts per {width}" + r" $\mu s$ Bin", fontproperties = prop_font, loc = "top")
    
    # ax = plt.gca()
    # for label in ax.get_xticklabels():
    #     label.set_fontproperties(prop_font)
    # for label in ax.get_yticklabels():
    #     label.set_fontproperties(prop_font)
    # ax.minorticks_on()
    # ax.get_xaxis().set_tick_params(which='both',direction='in', width=1)
    # ax.get_yaxis().set_tick_params(which='both',direction='in', width=1)
    # ax.xaxis.set_ticks_position('both')
    # ax.yaxis.set_ticks_position('both')
    
    # do the fit ... #
    def exponential(x, A, B, const):
        return A*np.exp(-const*x) + B
    
    popt, cov = curve_fit(exponential, mids, counts, bounds = ([0, 0, 0], [np.inf, np.inf, 5e-3]))
    
    fit_err = np.sqrt(np.diag(cov))
    X_FIT = np.arange(4, 1000, 1)
    Y_FIT = exponential(X_FIT, *popt)
    plt.plot(X_FIT, Y_FIT, color = "red", linewidth = 2)
    plt.legend(handles = [plt.plot([], [], color = "white")[0], plt.plot([], [], color = "white")[0], err, plt.plot([], linewidth =2, ls="-", color="red")[0], plt.plot([],[], color = "white")[0]], labels = ["SNO+ Preliminary", r"R $\leq$ 4.5 m", "BiPo214 Data", f"Exponential Fit\n" + r"$f(t) = Ae^{-\lambda t} + B$", f"A: {popt[0]:.0f}" + r" $\pm$" + f" {fit_err[0]:.2f}\nB: {popt[1]:.2f}" + r" $\pm$" + f" {fit_err[1]:.2f}\n" + r"$\lambda$" + f" : {popt[2] * 1000:.2f}" + r"$m s^{-1}$" + r" $\pm$" + f" {fit_err[2]*1000:.2f}" + r" $m s^{-1}$"], fancybox=False, numpoints=1,prop=prop_font,fontsize = 20, frameon=False)
    print(popt)
    plt.ylim((0, 100))
    plt.xlim((0, 1000))
    plt.savefig("../plots/plots_for_approval/bipo_tagging_dt.pdf")
    plt.close()

    # position #
    plt.hist2d(rho2, z_bi, bins = [np.arange(0, 1.1, 0.05), binning_pos], cmin = 1e-4)
    plt.xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$")
    plt.ylabel("Z (mm)")
    plt.savefig("../plots/plots_for_approval/bipo_tagging_rho2_z.pdf")
    plt.close()

    # position #
    plt.hist2d(x_bi, y_bi ,bins = binning_pos, cmin = 1e-4)
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.savefig("../plots/plots_for_approval/bipo_tagging_x_y.pdf")
    plt.close()

# def calculate_uncertainty(profile_ll, signal_hypothesis):
#     """
#     Function returns the 1 sigma confidence intervals on the minimum of profile LL
#     fit.
#     """
    
#     # find idx of point on every log-likelihood curve closes in value to 1
#     distance_to_interval = np.abs(profile_ll - 1) # absolute value of the difference
#     closest_to_interval  = np.argmin(distance_to_interval, axis = 1) # find first intercept with 1 sigma level

#     # create a masked array and remove the first intercept
#     masked_1sig_full = np.ma.masked_array(distance_to_interval, mask = False)

#     # mask out the closest interval in each row
#     masked_1sig_full.mask[0, closest_to_interval[0]] = True
#     masked_1sig_full.mask[1, closest_to_interval[1]] = True
#     masked_1sig_full.mask[2, closest_to_interval[2]] = True

#     # find second intercept
#     second_closest = np.argmin(masked_1sig_full, axis = 1)

#     # find minimum LL values as variables to save typing
#     combined_min   = signal_hypothesis[np.argmin(profile_ll[0,:])]
#     combined_upper = abs(combined_min - signal_hypothesis[second_closest[0]]) 
#     combined_lower = abs(combined_min - signal_hypothesis[closest_to_interval[0]])
#     multi_min      = signal_hypothesis[np.argmin(profile_ll[1,:])]
#     multi_upper    = abs(multi_min - signal_hypothesis[second_closest[1]]) 
#     multi_lower    = abs(multi_min - signal_hypothesis[closest_to_interval[1]])
#     energy_min     = signal_hypothesis[np.argmin(profile_ll[2,:])]
#     energy_upper   = abs(energy_min - signal_hypothesis[second_closest[2]]) 
#     energy_lower   = abs(energy_min - signal_hypothesis[closest_to_interval[2]])

#     return [combined_min, combined_upper, combined_lower], [multi_min, multi_upper, multi_lower], [energy_min, energy_upper, energy_lower]

def constrained_unconstrained_ll():
    """
    profileLL curves with constrained and unconstrained curves.
    """

    constrained       = np.load("./profileLL_constrained.npy")
    err_constrained   = np.load("./error_constrained.npy") 
    unconstrained     = np.load("./profileLL_unconstrained.npy")
    err_unconstrained = np.load("./error_unconstrained.npy")
    print(err_constrained)
    # Energy 2 multisite 1 combined 0
    signal_hypothesis = np.arange(0, 400, 1)
                                                                               # ${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$"
    plt.plot(constrained[2,:], color = "orange", label = rf"Energy Constrained: ${err_constrained[2, 0]}^{{+{err_constrained[2, 1]}}}_{{-{err_constrained[2,2]}}}$")
    plt.plot(unconstrained[2,:], color = "orange", linestyle = "dotted", label = rf"Energy Unconstrained: ${err_unconstrained[2, 0]}^{{+{err_unconstrained[2, 1]}}}_{{-{err_unconstrained[2,2]}}}$")

    plt.plot(constrained[1,:], color = "green", label = rf"Multisite Constrained: ${err_constrained[1, 0]}^{{+{err_constrained[1, 1]}}}_{{-{err_constrained[1,2]}}}$")
    plt.plot(unconstrained[1,:], color = "green", linestyle = "dotted", label = rf"Multisite Unconstrained: ${err_unconstrained[1, 0]}^{{+{err_unconstrained[1, 1]}}}_{{-{err_unconstrained[1,2]}}}$")

    plt.plot(constrained[0,:], color = "black", label = rf"Combined Constrained: ${err_constrained[0, 0]}^{{+{err_constrained[0, 1]}}}_{{-{err_constrained[0,2]}}}$")
    plt.plot(unconstrained[0,:], color = "black", linestyle = "dotted", label = rf"Combined Unconstrained: ${err_unconstrained[0, 0]}^{{+{err_unconstrained[0, 1]}}}_{{-{err_unconstrained[0,2]}}}$")

    

    plt.xlim((0, 50))
    plt.ylim((0, 10))

    plt.axhline(1, color = '#5e5e5e', linestyle = "dashed", label = r"1 $\sigma$")
    plt.axhline(4, color = '#959595', linestyle = "dashed", label = r"2 $\sigma$")
    plt.axhline(9, color = '#cecece', linestyle = "dashed", label = r"3 $\sigma$")
    plt.axvline(30.7, color = "red", label = r"$^8 B $ Predicted Value: 30.7")
    plt.legend(frameon = True, loc = "upper left")
    plt.savefig("../plots/constrained_unconstrained.pdf")

def tl208_profiles():
    
    back_hypothesis = np.arange(1, 400, 1)
    energy_curve = np.load("./tl208_profile_energy_constrained.npy")
    multi_curve  = np.load("./tl208_profile_multi_constrained.npy")
    combi_curve  = np.load("./tl208_profile_combi_constrained.npy")

    profile = np.zeros((3, len(back_hypothesis)))
    profile[0, :] = energy_curve
    profile[1, :] = multi_curve
    profile[2, :] = combi_curve

    energy, multi, combi = calculate_uncertainty(profile, back_hypothesis)
    print(energy)
    print(multi)
    print(combi)

    plt.plot(energy_curve, color = "orange", label = rf"Energy: ${energy[0]}^{{+{energy[1]}}}_{{-{energy[2]}}}$")
    plt.plot(multi_curve, color = "green", label = rf"Multisite: ${multi[0]}^{{+{multi[1]}}}_{{-{multi[2]}}}$")
    plt.plot(combi_curve, color = "black", label = rf"Combined: ${combi[0]}^{{+{combi[1]}}}_{{-{combi[2]}}}$")
    plt.ylim((0, 10))
    plt.xlim((100, 300))
    plt.xlabel("Tl208 Normalisation")
    plt.legend()
    plt.savefig("../plots/tl208_profile_constrained.pdf")
def obtain_pdf(location, fv_cut, multisite_bins, energy_bins, run_list, energy_range, plot_name):
        """
        Extract the energy and multisite PDFs for a given isotope.
        Only use events that fall within event selection cuts (e.g. FV cut).
        
        Inputs:
            location, str, the path to a folder containing all the ntuple files for 
                           specific isotope (1 for each run in run_list)

            fv_cut, float, the fv cut to apply in mm

            multisite_bins, array, float bin edges of multisite PDF

            energy_bins, array, float bin edges of energy PDF

            run_list, array, run number of each MC simulation to extract

            energy_range, str, string giving branch name of MC files to load, based
            on the energy used to generate the multisite / tres PDFs.
        
        Output: multisite_counts, array float, normalised bin counts of multisite pdf
                energy_counts, array float, normalised bin counts of energy pdf
        """

        energies     = []
        multisites   = []
        missing_runs = [] # runs that don't exist / corrupted but included in run list
        for irun in run_list:

            try:
                # open the MC ntuple file
                file   = ROOT.TFile.Open(f"{location}/{irun}.root")
                ntuple = file.Get(energy_range)
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

                # event passed FV so add info to PDFs
                energies.append(ientry.energy)
                multisites.append(ientry.dlogL)
        print(f"Found {len(missing_runs)} missing runs:\n{missing_runs}")

        

        return energies, multisites
def plot_multisite_energy_pdfs():
    """
    Create a nice 2 plot subplot showing normalised multisite and energy PDFs
    for each isotope used in the energy + multisite fits.
    """
    pdf_runlist    = np.loadtxt("../runlists/full_test.txt", dtype = int)
    energy_bins    = np.arange(2.5, 5.05, 0.05)
    multisite_bins = np.arange(-1.36, -1.325, 0.0005)
    energy_string  = "2p5_5p0"
    path           = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
    
    # load the information
    b8_energies, b8_multisites           = obtain_pdf(f"{path}/full_analysis_B8_solar_nue", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    tl208_energies, tl208_multisites     = obtain_pdf(f"{path}/full_analysis_Tl208", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    tl210_energies, tl210_multisites     = obtain_pdf(f"{path}/full_analysis_Tl210", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    bipo214_energies, bipo214_multisites = obtain_pdf(f"{path}/full_analysis_BiPo214", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    bipo212_energies, bipo212_multisites = obtain_pdf(f"{path}/full_analysis_BiPo212", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    
    # create the plots of the pdfs
    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    fig.set_size_inches((8.27 , 11.69 * 0.3 ))
    axes[0].hist(b8_energies,      bins = energy_bins, linewidth = 2, histtype = "step", density = True)
    axes[0].hist(tl208_energies,   bins = energy_bins, linewidth = 2, histtype = "step", density = True)
    axes[0].hist(tl210_energies,   bins = energy_bins, linewidth = 2, histtype = "step", density = True)
    axes[0].hist(bipo214_energies, bins = energy_bins, linewidth = 2, histtype = "step", density = True)
    axes[0].hist(bipo212_energies, bins = energy_bins, linewidth = 2, histtype = "step", density = True)
    axes[0].set_xlabel("Reconstructed Energy (MeV)", fontproperties = prop_font, fontsize = 15)
    axes[0].plot([], [], color = "C0", label = r"$^8$B $\nu _e$" )
    axes[0].plot([], [], color = "C1", label = r"$^{208}$Tl" )
    axes[0].plot([], [], color = "C2", label = r"$^{210}$Tl" )
    axes[0].plot([], [], color = "C3", label = r"BiPo214" )
    axes[0].plot([], [], color = "C4", label = r"BiPo212" )
    axes[0].legend(fontsize = 6, frameon = False, prop=fm.FontProperties(fname=path2, size = 12))
    axes[0].tick_params(axis='both', labelrotation=0, labelsize = 8)
    
    ax = axes[0]
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)

    axes[1].hist(b8_multisites,      bins = multisite_bins, linewidth = 2, histtype = "step", density = True)
    axes[1].hist(tl208_multisites,   bins = multisite_bins, linewidth = 2, histtype = "step", density = True)
    axes[1].hist(tl210_multisites,   bins = multisite_bins, linewidth = 2, histtype = "step", density = True)
    axes[1].hist(bipo214_multisites, bins = multisite_bins, linewidth = 2, histtype = "step", density = True)
    axes[1].hist(bipo212_multisites, bins = multisite_bins, linewidth = 2, histtype = "step", density = True)
    axes[1].set_xlabel(r"$\Delta log(\mathcal{L})$", fontproperties = prop_font, fontsize = 15)
    axes[1].plot([], [], color = "C0", label = r"$^8$B $\nu _e$" )
    axes[1].plot([], [], color = "C1", label = r"$^{208}$Tl" )
    axes[1].plot([], [], color = "C2", label = r"$^{210}$Tl" )
    axes[1].plot([], [], color = "C3", label = r"BiPo214" )
    axes[1].plot([], [], color = "C4", label = r"BiPo212" )
    axes[1].legend(fontsize = 6, frameon = False, prop=fm.FontProperties(fname=path2, size = 12))
    axes[1].tick_params(axis='both', labelrotation=0, labelsize = 8)

    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    fig.tight_layout()
    plt.savefig(f"../plots/pdf_analytics/energy_multisites_{energy_string}.pdf")

    print(len(b8_energies), len(b8_multisites))

def marginalisation_profileLL():

    profile = np.load('profile_ll_asimov_1yr_fullROI_2.npy')

    multisite_ll = profile[0,:]
    combined     = profile[1,:]
    energy       = profile[2,:]

    print(multisite_ll.shape)





signal_hypothesis = np.arange(0, 800, 1)
plt.rcParams['xtick.major.pad'] = '6' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '6'
# create_time_residual_pdfs()
# constrained_unconstrained_ll()
# tl208_profiles()
# plot_multisite_energy_pdfs()
# create_dlogL_plot()
# plt.rcParams['xtick.major.pad'] = '12' ## change me if the axis labels overlap! ## 
# plt.rcParams['ytick.major.pad'] = '12'
# create_asimov_dataset_graphs()
# plt.rcParams['xtick.major.pad'] = '6' ## change me if the axis labels overlap! ## 
# plt.rcParams['ytick.major.pad'] = '6'
# create_profile_LL_asimov_plot()
# create_data_vs_mc_plot()
# impact_of_reprocessing()
# data_quality()
# bipo_tagging_data_plots()
marginalisation_profileLL()