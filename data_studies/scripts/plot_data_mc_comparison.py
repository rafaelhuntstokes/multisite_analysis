import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from scipy.optimize import curve_fit

def compare_dlogL_data_mc():
    """
    Function loads the data and MC dlog(L) information for B8 events above 5 MeV.
    Compares the dlog(L) distributions.
    """

    # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    # data = ROOT.TFile.Open(working_dir + "/extracted_data/above_5MeV/solar_data_discriminants/total.root")
    # mc   = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/B8_solar_nue/total.root")
    data = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/bismuth214_data_discriminants/total.root")
    mc   = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/Bi214/total.root")

    # grab the above 5MeV TTrees
    # data_5MeV_tree = data.Get("5p0_plus")
    # mc_5MeV_tree   = mc.Get("5p0_plus")
    data_5MeV_tree = data.Get("1p25_3p0")
    mc_5MeV_tree   = mc.Get("1p25_3p0")
    
    FV_CUT      = 5250
    ITR_LOW     = 0.21 
    ITR_HIGH    = 0.3
    E_CUT_LOW   = 1.25
    E_CUT_HIGH  = 3.0
    data_dlogL  = []
    mc_dlogL    = []
    data_itr    = []
    cos_theta   = []
    data_energy = []
    
    for ientry in data_5MeV_tree:
        if ientry.itr > ITR_LOW and ientry.itr < ITR_HIGH and ientry.energy >= E_CUT_LOW and ientry.energy < E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            data_dlogL.append(ientry.dlogL)
            data_itr.append(ientry.itr)
            cos_theta.append(ientry.cos_theta_sun)
            data_energy.append(ientry.energy)
    for ientry in mc_5MeV_tree:
        if ientry.itr > ITR_LOW and ientry.itr < ITR_HIGH and ientry.energy >= E_CUT_LOW and ientry.energy < E_CUT_HIGH:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            mc_dlogL.append(ientry.dlogL)

    # create output plot
    plt.figure()
    # binning = np.linspace(1.096, 1.12, 20)
    binning = np.linspace(-0.125, -0.11, 20)
    binning_mc = np.linspace(1.08, 1.12, 20)
    plt.title(f"{E_CUT_LOW} <= E < {E_CUT_HIGH} MeV")
    plt.xlabel("dlog(L)")
    plt.ylabel("Counts")
    counts_data, bins_data = np.histogram(data_dlogL, bins = binning)
    counts_mc, bins_mc     = np.histogram(mc_dlogL,  bins = binning)
    
    # scale the counts in the MC to match the counts in the data
    scale = np.sum(counts_data) / np.sum(counts_mc)
    error_mc = np.sqrt(counts_mc) * scale
    error_data = np.sqrt(counts_data)

    # plt.hist(data_dlogL, bins = binning, color = "black", histtype = "step", label = r"$^{214}Bi$ Data | Num Events: " + f"{len(data_dlogL)}")
    counts_mc, _, _ = plt.hist(mc_dlogL, color = "red", bins = binning, weights = np.full_like(mc_dlogL, scale), histtype="step", label = r"$^{214}Bi$ MC")
    # counts_mc, _, _ = plt.hist(mc_dlogL, color = "red", bins = binning, weights = np.full_like(mc_dlogL, scale), histtype="step", label = r"$^{8}B$ MC")
    print(counts_data, counts_mc)
    
    
    mids = bins_data[:-1] + np.diff(bins_data)[0]/2
    plt.errorbar(mids, counts_data, linestyle = "", marker = "o", markersize = 2, yerr = error_data, color = "black", capsize = 2, label = r"$^{214}Bi$ Data | Num Events: " + f"{len(data_dlogL)}")
    # plt.errorbar(mids, counts_data, linestyle = "", marker = "o", markersize = 2, yerr = error_data, color = "black", capsize = 2, label = r"$^{8}B$ Data | Num Events: " + f"{len(data_dlogL)}")
    plt.errorbar(mids, counts_mc, linestyle = "", yerr = error_mc, color = "red", capsize=2)

    print(np.sum(counts_mc), np.sum(counts_data))
    plt.legend(fontsize = 8, frameon = False)
    plt.savefig(f"../plots/bi214_dlogL_{E_CUT_LOW}_{E_CUT_HIGH}MeV_FV_{FV_CUT}mm.pdf")
    plt.close()

    # create 5 subplots showing dlog(L) data vs MC as a function of energy
    fig, axes = plt.subplots(nrows = 1, ncols = 4, figsize = (30, 10))

    # binning = np.linspace(1.096, 1.12, 15)
    # E_CUTS     = [6.0, 7.0, 8.0, 10.0, 12.0]
    E_CUTS       = [1.25, 1.75, 2.25, 2.75, 3.25]
    for i in range(len(E_CUTS)-1):

        data_dlogL = []
        mc_dlogL   = []
        
        # define the cuts to apply
        E_LOW  = E_CUTS[i]
        E_HIGH = E_CUTS[i+1]
        print("Calculating dlog(L) agreement for energy range: ", E_LOW, E_HIGH)

        # data
        for ientry in data_5MeV_tree:
            itr    = ientry.itr
            energy = ientry.energy
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            if itr < ITR_LOW or itr > ITR_HIGH:
                continue
            if energy <= E_LOW or energy > E_HIGH:
                continue
            
            data_dlogL.append(ientry.dlogL)
        # mc 
        for ientry in mc_5MeV_tree:
            itr    = ientry.itr
            energy = ientry.energy
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            if itr < ITR_LOW or itr > ITR_HIGH:
                continue
            if energy <= E_LOW or energy > E_HIGH:
                continue
            
            mc_dlogL.append(ientry.dlogL)

        # create the plot #
        
        # obtain mc normalisation factor
        counts_data, bins_data = np.histogram(data_dlogL, bins = binning)
        counts_mc, bins_mc     = np.histogram(mc_dlogL,  bins = binning)
        scale = np.sum(counts_data) / np.sum(counts_mc)

        error_data = np.sqrt(counts_data)
        error_mc   = np.sqrt(counts_mc) * scale
        # axes[i].hist(data_dlogL, linewidth = 3, bins = binning, color = "black", histtype = "step", label = r"$^8B$ Data" + f"\nNum Events:{len(data_dlogL)}")
        axes[i].errorbar(mids, counts_data, linestyle = "", marker = "o", markersize = 2, yerr = error_data, color = "black", capsize = 2, label = r"$^{214}Bi$ Data | Num Events: " + f"{len(data_dlogL)}")
        # axes[i].errorbar(mids, counts_data, linestyle = "", marker = "o", markersize = 2, yerr = error_data, color = "black", capsize = 2, label = r"$^{8}B$ Data | Num Events: " + f"{len(data_dlogL)}")
        counts_mc, _, _  = axes[i].hist(mc_dlogL, linewidth = 3, color = "red", bins = binning, weights = np.full_like(mc_dlogL, scale), histtype="step", label = r"$^{214}Bi$ MC")
        # counts_mc, _, _  = axes[i].hist(mc_dlogL, linewidth = 3, color = "red", bins = binning, weights = np.full_like(mc_dlogL, scale), histtype="step", label = r"$^{8}B$ MC")
        axes[i].errorbar(mids, counts_mc, linestyle = "", yerr = error_mc, color = "red", capsize=2)
        axes[i].set_title(f"{E_LOW} <= E < {E_HIGH} MeV", fontsize = 20)
        axes[i].set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = 20)
        axes[i].set_ylabel("Counts", fontsize = 20)
        axes[i].legend(fontsize = 15, loc = "upper left", frameon = False)

    fig.tight_layout()
    plt.savefig(f"../plots/bi214_dlogL_energy_comparison_FV_{FV_CUT}mm.pdf")

def check_candidate_health():
    """
    Function used to plot the reconstructed position, energy, itr
    and time residuals of each extracted B8 solar candidate in data.
    
    The individual event residuals are saved in a subfolder as .png files.
    These are used to check the events look 'normal' and aren't 
    'neck hotspot trash'.
    """
    
    # load the data and mc
    working_dir = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    # data           = ROOT.TFile.Open(working_dir + "/extracted_data/above_5MeV/solar_data_discriminants/total.root")
    data           = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/bismuth214_data_discriminants/total.root")
    data_5MeV_tree = data.Get("1p25_3p0")
    # data_5MeV_tree = data.Get("5p0_plus")
    
    energy  = []
    itr     = []
    x       = []
    y       = []
    z       = []
    dlogl   = []
    
    FV_CUT      = 5250
    ITR_LOW     = 0.21
    ITR_HIGH    = 0.3
    E_CUT_LOW   = 1.25
    E_CUT_HIGH  = 3.00

    for ientry in data_5MeV_tree:
        if ientry.itr > ITR_LOW and ientry.itr < ITR_HIGH and ientry.energy >= E_CUT_LOW and ientry.energy < E_CUT_HIGH:
            X = ientry.x
            Y = ientry.y
            Z = ientry.z
            r = np.sqrt(X**2 + Y**2 + Z**2)
            if r > FV_CUT:
                continue
            
            energy.append(ientry.energy)
            itr.append(ientry.itr)
            x.append(ientry.x)
            y.append(ientry.y)
            z.append(ientry.z)
            dlogl.append(ientry.dlogL)
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    rho2 = (x**2 + y**2) / 6000**2
    r = np.sqrt(x**2 + y**2 + z**2)
    
    fig, axes = plt.subplots(nrows = 2, ncols = 5, figsize = (30, 10))
    
    energy_bins = np.arange(1.0, np.max(energy) +2, 0.1) 
    itr_bins    = np.arange(0.1, 0.4, 0.01)
    rho2_bins   = np.linspace(0, 1, 50)
    z_bins      = np.linspace(-6000, 6000, 50)
    dlogl_bins  = np.linspace(-0.125, -0.11, 20)
    r_bins      = np.linspace(0, 6100, 50)
    # energy_bins = np.arange(3, np.max(energy) +2, 0.5) 
    # itr_bins    = np.arange(0.1, 0.4, 0.01)
    # rho2_bins   = np.linspace(0, 1, 50)
    # z_bins      = np.linspace(-6000, 6000, 50)
    # dlogl_bins  = np.linspace(1.1, 1.12, 20)
    # r_bins      = np.linspace(0, 6100, 50)
    axes[0,0].hist(energy, bins = energy_bins, histtype = "step")
    axes[0,0].set_xlabel("Reconstructed Energy (MeV")
    
    axes[0,1].hist(itr, bins = itr_bins, histtype = "step")
    axes[0,1].set_xlabel("ITR")
    
    axes[0,2].hist(dlogl, bins = dlogl_bins, histtype = "step")
    axes[0,2].set_xlabel("dlogL Discriminant")
    
    axes[0,3].hist2d(x, y, bins = z_bins, cmin = 1e-4)
    axes[0,3].set_xlabel("Reconstructed X (mm)")
    axes[0,3].set_ylabel("Reconstructed Y (mm)")
    
    axes[0,4].hist2d(rho2, z, bins = [rho2_bins, z_bins], cmin = 1e-4)
    axes[0,4].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$")
    axes[0,4].set_ylabel("Reconstructed Z (mm)")
    
    axes[1,0].hist(r, bins = r_bins, histtype = "step")
    axes[1,0].set_xlabel("Reconstructed R (mm)")
    
    axes[1,1].hist2d(energy, dlogl, bins = [energy_bins, dlogl_bins], cmin = 1e-8)
    axes[1,1].set_xlabel("Reconstructed Energy (MeV)")
    axes[1,1].set_ylabel("dLog(L) Discriminant")
    
    axes[1,2].hist2d(itr, dlogl, bins = [itr_bins, dlogl_bins], cmin = 1e-8)
    axes[1,2].set_xlabel("ITR")
    axes[1,2].set_ylabel("dLog(L) Discriminant")
    
    axes[1,3].hist2d(r, dlogl, bins = [r_bins, dlogl_bins], cmin = 1e-8)
    axes[1,3].set_xlabel("Reconstructed R (mm)")
    axes[1,3].set_ylabel("dlog(L) Discriminant")

    axes[1,4].hist2d(energy, itr, bins = [energy_bins, itr_bins], cmin = 1e-8)
    axes[1,4].set_xlabel("Reconstructed Energy (MeV)")
    axes[1,4].set_ylabel("ITR")
    
    fig.tight_layout()
    plt.subplots_adjust(top=0.93)
    # plt.suptitle(f"B8 Candidates | Num. EVs: {len(dlogl)}", fontsize = 30)
    plt.suptitle(f"Bi214 Candidates | Num. EVs: {len(dlogl)}", fontsize = 30)
    plt.savefig("../plots/bi214_healthcheck_clean.pdf")
    # plt.savefig("../plots/solar_healthcheck_clean.png")

def get_bipo_info():
    """
    Function is designed to check the quality of BiPo214 tags per energy bin used
    in analysis of the data vs mc dlog(L) matching.

    It is thought that the higher energy bin might be contaminated with some
    'bad stuff'.

    So, for each energy bin, how bad is the dT and dR distributions?
    """

    def exponential_fit(x, A, B, tau):
        return A * np.exp(-x/tau) + B
    
    # load the run list
    runlist = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    bipo_dir = "/data/snoplus3/hunt-stokes/clean_bipo/tagged_ntuples/run_by_run214"
    
    # for each run in the list, load the BiPo tagged NTUPLE information

    energy = []
    dT     = []
    dR     = []
    for irun in runlist:

        file   = ROOT.TFile.Open(f"{bipo_dir}/output_{irun}.root")
        ntuple = file.Get("ntuple214")

        for ientry in ntuple:
            x = ientry.x2
            y = ientry.y2
            z = ientry.z2
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > 5250:
                continue

            energy.append(ientry.energy2)
            dT.append(ientry.dT)
            dR.append(ientry.dR)

    # convert to numpy arrays for easy plotting of info in energy bins
    energy = np.array(energy)
    dT     = np.array(dT)
    dR     = np.array(dR)

    # setup the plots
    fig, axes   = plt.subplots(nrows = 2, ncols = 4, figsize = (20, 10))
    energy_bins = [1.25, 1.75, 2.25, 2.75, 3.25]
    binning_r   = np.arange(0, 1500, 50)
    binning_t   = np.arange(3, 1000, 20) 
    for i in range(len(energy_bins) -1):
        E_LOW  = energy_bins[i]
        E_HIGH = energy_bins[i+1]

        # find the idx of dT and dR where energy lies between these bins
        idx_plot = np.where((energy >= E_LOW) & (energy < E_HIGH))[0]
        print(idx_plot)
        dT_plot = dT[idx_plot] / 1e3

        # plot
        counts, bin_edge, _ = axes[0, i].hist(dT_plot, histtype = "step", bins = binning_t)
        mids = bin_edge[:-1] + np.diff(bin_edge)[0] / 2
        lifetime = 164 / np.log(2)

        # do a fit in the dT
        popt, cov, = curve_fit(exponential_fit, mids, counts, p0 = [100, 10, lifetime], bounds = ([0, 0, 150], [np.inf, np.inf, 350]))
        
        error = np.sqrt(np.diag(cov))
        XFIT = np.arange(0, 1000, 1)
        YFIT = exponential_fit(XFIT, popt[0], popt[1], popt[2])
        axes[0, i].plot(XFIT, YFIT)
        axes[0, i].plot([], [], "", label = f"Fit Parameters: A = {round(popt[0], 3)} +- {round(error[0], 3)}\n Lifetime: {round(popt[2], 3)} +- {round(error[2], 3)} mu s\n B = {round(popt[1], 3)} +- {round(error[1], 3)}")
        axes[0, i].set_title(f"{E_LOW} <= E < {E_HIGH} MeV")
        axes[0, i].set_xlabel(r"$\Delta T$ ($\mu s$)")
        axes[0, i].set_xlim((0, 1000))
        axes[0, i].legend(frameon = False)
        axes[1, i].hist(dR[idx_plot], histtype = "step", bins = binning_r)
        axes[1, i].set_xlabel("$\Delta R$ (mm)")
        axes[1, i].set_xlim((0, 1000))

    plt.suptitle(r"$^{214}$Bi Candidates", fontsize = 20)
    plt.savefig("../plots/bipo_dT_dR_vsE.pdf")

get_bipo_info()
# check_candidate_health()
# compare_dlogL_data_mc()
# compare_single_site_data_mc()
