import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def create_roc(likelihood_signal, likelihood_background):

	# bin the data 
	min_vbb = np.amin(likelihood_signal)
	max_vbb = np.amax(likelihood_signal)
	min_cobalt = np.amin(likelihood_background)
	max_cobalt = np.amax(likelihood_background)
	if min_vbb < min_cobalt:
		bin_min = min_vbb
	else: 
		bin_min = min_cobalt
	if max_vbb > max_cobalt:
		bin_max = max_vbb
	else: 
		bin_max = max_cobalt
	print("Max Bin: {}\nMin Bin: {}".format(bin_max, bin_min))
	steps = 100#int((bin_max-bin_min)/0.2)
	print(steps)
	binning = np.linspace(bin_min, bin_max, steps)
	sig_counts, sig_bins = np.histogram(likelihood_signal, density=True, bins = binning)
	back_counts, back_bins = np.histogram(likelihood_background, density=True, bins = binning)

	max_hist = np.argmax(sig_counts)
	print("Number entries in signal hist: {}\nNum entries in background his: {}".format(np.size(likelihood_signal), np.size(likelihood_background)))
	# plt.figure()
	# plt.hist(likelihood_signal, density=True, histtype="step", bins = binning, label = "signal", linewidth=2)
	# plt.hist(likelihood_background, density=True, histtype="step", bins = binning, label = "background", linewidth=2)
	# plt.title("Log-Likelihood Ratio")
	# plt.legend()
	
	signal_acceptance = []
	background_acceptance = []
	
    # make a cut at each bin 
	tot_sig = np.sum(sig_counts)
	tot_back = np.sum(back_counts)
	print("LEN BINNING: ", len(binning[:-1])) 
	for i in range(len(binning[:-1])):
		# calculate the % of bin counts to the LEFT of the cut (hopefully the same as integrating it!)
		signal_accepted = np.sum(sig_counts[i:] * np.diff(sig_bins)[0]) 
		background_accepted = np.sum(back_counts[i:] * np.diff(back_bins)[0]) 
		signal_acceptance.append(signal_accepted)
		background_acceptance.append(background_accepted)

	return signal_acceptance, background_acceptance, binning

def analysis_output_plot():
    """
    Function creates an absolutely massive (5 x 4) plot array showing: 
        1. dlog(L) for each energy bin
        2. ROC for multisite for each energy bin
        3. cos(theta_sun) for each energy bin
        4. Scatter of (cos(theta_sun), dlog(L)) for each energy bin 
    """

    # create the mega grid of plots
    fig, axes = plt.subplots(nrows = 2, ncols = 5, figsize = (30, 14))
    font_s    = 30
    
    bins_dir = np.linspace(-1, 1, 40)
    # load each of the test sets
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies"
    test_B8     = ROOT.TFile.Open(working_dir + "/run_by_run_test/B8_solar_nue/total_combined.root")
    test_Tl208  = ROOT.TFile.Open(working_dir + "/run_by_run_test/Tl208/total_combined.root")

    # get each of the histograms
    histo_names = ["2p5_5p0", "2p5_3p125", "3p125_3p75", "3p75_4p375", "4p375_5p0"]
    plot_titles = [r"$2.5 \leq E < 5.0$ MeV", r"$2.5 \leq E < 3.125$ MeV", r"$3.125 \leq E < 3.75$ MeV", r"$3.75 \leq E < 4.375$ MeV", r"$4.375 \leq E < 5.0$ MeV"]
    plot_col    = 0
    for name in histo_names:
        
        dlogL_B8      = []
        dlogL_Tl208   = []
        cos_sun_B8    = []
        cos_sun_Tl208 = []
        roc_multi     = []

        # fill in the plots for each histogram in turn
        tree_b8    = test_B8.Get(name)
        tree_tl208 = test_Tl208.Get(name)
        for ientry in tree_b8:
            dlogL_B8.append(ientry.dlogL)
            cos_sun_B8.append(ientry.cos_theta_sun)
        for ientry in tree_tl208:
            dlogL_Tl208.append(ientry.dlogL)
            cos_sun_Tl208.append(ientry.cos_theta_sun)
        
        # compute the ROC curve
        signal_roc, background_roc, binning = create_roc(dlogL_B8, dlogL_Tl208)
		
        # fill the plot!
        axes[0, plot_col].hist(dlogL_B8, bins = binning, density = True, histtype = "step", label = r"$^8B$ | Num Evs: " + f"{len(dlogL_B8)}")
        axes[0, plot_col].hist(dlogL_Tl208, bins = binning, density = True, histtype = "step", label = r"$^{208}Tl$ | Num. Evs: " + f"{len(dlogL_Tl208)}")
        axes[0, plot_col].legend(frameon=False, fontsize = 14, loc = "upper left")
        axes[0, plot_col].set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[0, plot_col].set_ylabel("Normalised Counts", fontsize = font_s)
        axes[0, plot_col].set_title(plot_titles[plot_col], fontsize = font_s)
        
        axes[1, plot_col].plot(signal_roc, background_roc, linewidth = 2)
        axes[1, plot_col].set_xlabel("Signal Acceptance", fontsize = font_s)
        axes[1, plot_col].set_ylabel("Background Acceptance", fontsize = font_s)
        axes[1, plot_col].grid("both")
        axes[1, plot_col].set_xlim((0,1))
        axes[1, plot_col].set_ylim((0,1))
        # axes[2, plot_col].hist(cos_sun_B8, bins = bins_dir, density = True, histtype = "step", label = r"$^8B$")
        # axes[2, plot_col].hist(cos_sun_Tl208, bins = bins_dir, density = True, histtype = "step", label = r"$^{208}Tl$")
        # axes[2, plot_col].legend(frameon=False, fontsize = font_s, loc = "upper left")
        # axes[2, plot_col].set_xlabel(r"$cos(\theta_{sun})$", fontsize = font_s)
        # axes[2, plot_col].set_ylabel("Normalised Counts", fontsize = font_s)

        # axes[3, plot_col].scatter(dlogL_B8, cos_sun_B8, alpha = 0.1, s = 1, color = "red", label = r"$^8B$")
        # axes[3, plot_col].scatter(dlogL_Tl208, cos_sun_Tl208, alpha = 0.1, s = 1, color = "blue", label = r"$^{208}Tl$")
        # axes[3, plot_col].legend(frameon = False, fontsize = font_s, loc = "upper left")

        plot_col += 1
    fig.tight_layout()
    plt.savefig(working_dir + "/plots/results/mc_energy_binned_discrimination_5.pdf")


def plot_itr_as_discriminant():
    """
    Noticed that the ITR is correlated with dlog(L). Can you use it is a discriminant?
    Checking the ITR for B8 and Tl208 MC.
    """
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies"
    test_B8     = ROOT.TFile.Open(working_dir + "/run_by_run_test/B8_solar_nue/total.root")
    test_Tl208  = ROOT.TFile.Open(working_dir + "/run_by_run_test/Bi214/total.root")
	
    tree_b8     = test_B8.Get("2p5_5p0")
    tree_tl208  = test_Tl208.Get("2p5_5p0")
    itr_b8      = []
    itr_tl208   = []
    dlogL_b8    = []
    dlogL_tl208 = []
    for ientry in tree_b8:
        dlogL_b8.append(ientry.dlogL)
        itr_b8.append(ientry.itr)
    for ientry in tree_tl208:
        dlogL_tl208.append(ientry.dlogL)
        itr_tl208.append(ientry.itr)

    fig, axes  = plt.subplots(nrows = 2, ncols = 1)
    bins_itr   = np.linspace(0.18, 0.3, 100)
    bins_dlogL = np.linspace(-1.115, -1.095, 100) 
    axes[0].hist(dlogL_b8, bins = bins_dlogL, histtype = "step", label = "B8")
    axes[0].hist(dlogL_tl208, bins = bins_dlogL, histtype = "step", label = "Bi214")
    axes[0].set_xlabel("dLog(L)")
    axes[0].legend()
    
    axes[1].hist(itr_b8, bins = bins_itr, histtype = "step", label = "B8")
    axes[1].hist(itr_tl208, bins = bins_itr, histtype = "step", label = "Bi214")
    axes[1].legend()
    axes[1].set_xlabel("ITR")
    fig.tight_layout()
    plt.savefig("../plots/itr_dlogL.pdf")

plot_itr_as_discriminant()
# analysis_output_plot()
