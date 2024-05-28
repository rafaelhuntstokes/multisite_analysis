import numpy as np
import os
import rat
from ROOT import RAT
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d, binned_statistic

# import plotly.express as px
"""
Script is created to investigate the systematic shift in multisite discriminant
observed in Bi214 data vs MC.

This systematic disappears when an ITR cut is imposed on the dataset:

Only accepting events with 0.21 < ITR removes the systematic.

Want to isolate the Bi214 events failing this cut, and work out if they look 
'weird' in time residual etc. in order to filter them out / work out what is
going on with them a bit better.

Will first identify the GTIDs of Bi214 events failing the ITR cut, and create a
runlist containing them.

Then, running through this run list, save the time residual distributions, position,
energy, ITR, multisite etc. and see if there's any patters.

The ITR might show up neck hotspots (flat in time residual), which would be nice
to exclude easily ...
"""

def write_gtids_of_weird_events():
    """
    Function creates a .txt file per run containing the GTIDs of Bi214 events
    failing the ITR < 0.21 cut.
    """

    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    run_list        = np.loadtxt(working_dir + "/runlists/contains_solar_candidates.txt", dtype = int)

    # loop over all the Bi214 events and extract the itr failed events
    for irun in run_list:
        print("Run ", irun)
        file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/{irun}.root")
        tree = file.Get("bi214")

        # track the Bi214 events that have low itr
        failed_gtids = []
        for ientry in tree:

            itr  = ientry.itr
            gtid = ientry.gtid

            if itr < 0.21:
                failed_gtids.append(gtid)
        # write these gtids to a run by run gtid .txt file
        np.savetxt(working_dir + f"/extracted_data/bi214/gtid_lists/low_itr/{irun}.txt", failed_gtids, fmt = "%i")

def plot_time_residuals():
    """
    Function will make a plot of the individual time residuals per event per run
    that have been identified as having low ITR.
    """

    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    run_list        = np.loadtxt(working_dir + "/runlists/contains_solar_candidates.txt", dtype = int)

    for irun in run_list:
        file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_lowITR/{irun}.root")

        branch = file.Get("hit_level")
        ev_count = 1

        event_itr  = []
        event_multi = []
        event_x    = []
        event_y    = []
        event_z    = []
        event_e    = []
        event_rho2 = []
        for ientry in branch:
            residuals = list(ientry.tRes)
            itr       = ientry.itr
            multi     = ientry.dlogL
            x         = ientry.x
            y         = ientry.y
            z         = ientry.z
            rho2      = x**2 + y**2
            e         = ientry.energy

            event_itr.append(itr)
            event_multi.append(multi)
            event_x.append(x)
            event_y.append(y)
            event_z.append(z)
            event_rho2.append(rho2)
            event_e.append(e)

            plt.figure()
            plt.hist(residuals, bins = np.arange(-100, 300, 1), histtype = "step", density = False)
            plt.xlabel("Time Residual (ns)")
            plt.ylabel("Area Normalised Counts")
            plt.ylim((0, 50))
            plt.legend(handles = [plt.plot([], [], linestyle = "")[0]], labels = [f"ITR: {itr:.3f}"])
            plt.savefig(f"../plots/low_itr_bi214/tres/{irun}_{ev_count}.png")
            plt.close()

            ev_count += 1
        
        plt.figure()
        plt.hist2d(event_itr, event_multi, bins = 50)
        plt.xlabel("ITR")
        plt.ylabel("Multisite")
        plt.tight_layout()
        plt.savefig(f"../plots/low_itr_bi214/{irun}_itr_vs_multi.png")
        plt.close()
        break

def plot_itr_vs_multi():
    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    run_list        = np.loadtxt(working_dir + "/runlists/contains_solar_candidates.txt", dtype = int)
    event_itr  = []
    event_multi = []
    event_x    = []
    event_y    = []
    event_z    = []
    event_e    = []
    event_rho2 = []
    event_r    = []
    removed_itr_r = []
    removed_itr_x = []
    removed_itr_y = []
    removed_itr_z = []
    # for irun in run_list:
    #     print(irun)
        
    #     try:
    #         # file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/discriminants_7.0._2.5_3.125_PDF/{irun}.root")
    #         # file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/{irun}.root")
    #         file = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/bismuth214_data_discriminants/total.root")
    #         branch = file.Get("1p25_3p0")
    #     except:
    #         continue
        
    #     for ientry in branch:
    #         itr       = ientry.itr
    #         multi     = ientry.dlogL
    #         x         = ientry.x
    #         y         = ientry.y
    #         z         = ientry.z
    #         rho2      = x**2 + y**2
    #         e         = ientry.energy

    #         event_itr.append(itr)
    #         event_multi.append(multi)
    #         event_x.append(x)
    #         event_y.append(y)
    #         event_z.append(z)
    #         event_rho2.append(rho2)
    #         event_e.append(e)
        
    # file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/discriminants_7.0._2.5_3.125_PDF/{irun}.root")
    file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/total.root")
    # file = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/bismuth214_data_discriminants/total.root")
    branch = file.Get("bi214")
    
    FV_CUT = 6000
    E_LOW  = 1.25
    E_HIGH = 3.0
    ITR_LOW = 0.21
    
    for ientry in branch:
        itr       = ientry.itr
        multi     = ientry.dlogL
        x         = ientry.x
        y         = ientry.y
        z         = ientry.z
        rho2      = np.sqrt(x**2 + y**2)
        e         = ientry.energy
        if e < E_LOW or e > E_HIGH:
            continue
        r = np.sqrt(x**2 + y**2 + z**2)

        if r > FV_CUT:
            continue
        # if z > 5000:
        #     continue
        # if rho2 > 2000:
        #     continue
        if itr < ITR_LOW:
            removed_itr_r.append(r)
            removed_itr_x.append(x)
            removed_itr_y.append(y)
            removed_itr_z.append(z)
            continue
        event_itr.append(itr)
        event_multi.append(multi)
        event_x.append(x)
        event_y.append(y)
        event_z.append(z)
        event_rho2.append(rho2)
        event_e.append(e)
        event_r.append(r)
        
    # plot the MC versions
    mc       = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/Bi214/total.root")
    mc_tree  = mc.Get("1p25_3p0")
    mc_itr   = []
    mc_multi = []
    mc_x     = []
    mc_y     = []
    mc_z     = []
    mc_e     = []
    mc_rho2  = []
    mc_r     = []
    
    for ientry in mc_tree:
        itr       = ientry.itr
        multi     = ientry.dlogL
        x         = ientry.x
        y         = ientry.y
        z         = ientry.z
        r         = np.sqrt(x**2 + y**2 + z**2)
        if r > 5250:
            print(r)
        rho2 = np.sqrt(x**2+y**2)
        if r > FV_CUT:
            continue
        e         = ientry.energy

        if e < E_LOW or e > E_HIGH:
            continue
        if itr < ITR_LOW:
            continue
        # if rho2 > 2000:
        #     continue
        mc_itr.append(itr)
        mc_multi.append(multi)
        mc_x.append(x)
        mc_y.append(y)
        mc_z.append(z)
        mc_rho2.append(rho2)
        mc_e.append(e)
        mc_r.append(r)
    print(len(event_itr))
    # fig, axes = plt.subplots(nrows = 1, ncols = 2, sharex=True, sharey=True)
    # axes[0].hist2d(event_r, event_multi, density = True, bins = 10, cmin = 1e-4)
    # axes[0].set_xlabel("Radius")
    # axes[0].set_ylabel("Multisite")
    # axes[0].set_title("Data")

    # axes[1].hist2d(mc_r, mc_multi, density = True, bins = 10, cmin = 1e-4)
    # axes[1].set_xlabel("Radius")
    # axes[1].set_ylabel("Multisite")
    # axes[1].set_title("MC")
    # fig.tight_layout()
    # plt.savefig(f"../plots/low_itr_bi214/r_vs_multi_7015.png")
    # plt.close()
    
    plt.figure()
    binning = np.arange(-0.130, -0.105, 0.0005)
    binning = np.linspace(-0.125, -0.11, 20)
    mids = binning[:-1] + np.diff(binning)[0] / 2 
    counts_data, bins_data = np.histogram(event_multi, bins = binning)
    counts_mc, bins_mc     = np.histogram(mc_multi,  bins = binning)
    
    # scale the counts in the MC to match the counts in the data
    scale = np.sum(counts_data) / np.sum(counts_mc)
    error_mc = np.sqrt(counts_mc) * scale
    error_data = np.sqrt(counts_data)
    # plt.hist(event_multi, bins = binning, histtype = "step", density = True, color = "black", label = "data")
    # plt.hist(mc_multi, bins = binning, histtype = "step", density = True, color = "red", label = "MC")
    plt.errorbar(mids, counts_data, yerr = error_data, linestyle = "", color = "black", capsize = 2, label = "data")
    plt.hist(mc_multi, color = "red", bins = binning, weights = np.full_like(mc_multi, scale), histtype="step", label = r"$^{214}Bi$ MC")
    plt.xlabel("Multisite")
    plt.title(f"FV: {FV_CUT} mm | {E_LOW} < E < {E_HIGH} MeV")
    plt.legend()
    plt.savefig(f"../plots/low_itr_bi214/multisite_comparisons/data_vs_mc_multi_{FV_CUT}mm_{E_LOW}_{E_HIGH}MeV.png")
    plt.close()

    plt.figure()
    plt.hist(mc_multi, bins = binning, density = True, histtype = "step", color = "red")
    plt.hist(event_multi, bins = binning, density = True, histtype = "step", color = "black")
    plt.savefig("../plots/low_itr_bi214/zcut_multisite.pdf")
    plt.close()
    
    plt.figure()
    plt.hist2d(event_rho2, event_z, bins = 50, density = True, cmin = 1e-9)
    plt.xlabel("rho")
    plt.ylabel("Z")
    plt.savefig("../plots/low_itr_bi214/rho_vs_z.png")
    plt.close()

    plt.figure()
    plt.hist2d(event_x, event_y, bins = 50, density = True, cmin = 1e-9)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig("../plots/low_itr_bi214/x_vs_z.png")
    plt.close()
    # plt.figure()
    # plt.hist2d(event_itr, event_e, bins = 50, cmin = 1e-4)
    # plt.xlabel("ITR")
    # plt.ylabel("Energy")
    # plt.tight_layout()
    # plt.savefig(f"../plots/low_itr_bi214/{irun}_itr_vs_energy_708.png")
    # plt.close()

    # plt.figure()
    bins_rho = np.linspace(0, 6500, 11)
    bins_z   = np.linspace(-6500, 6500, 11)
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (12, 12))
    mean_itr, x_edge, y_edge, _ = binned_statistic_2d(event_z, event_rho2, event_itr, bins = [bins_z, bins_rho])
    print(mean_itr.shape)
    img = axes[0,0].imshow(mean_itr, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = 0.15, vmax = 0.35)
    plt.colorbar(img, ax = axes[0,0])
    axes[0,0].set_xlabel("Rho (mm)")
    axes[0,0].set_ylabel("Z (mm)")
    axes[0,0].set_title("ITR")

    mean_multi, x_edge, y_edge, _ = binned_statistic_2d(event_z, event_rho2, event_multi, bins = [bins_z, bins_rho])
    print(mean_itr.shape)
    img = axes[0,1].imshow(mean_multi, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = -0.120, vmax = -0.110)
    plt.colorbar(img, ax = axes[0,1])
    axes[0,1].set_xlabel("Rho (mm)")
    axes[0,1].set_ylabel("Z (mm)")
    axes[0,1].set_title("Multisite")

    mean_itr, x_edge, y_edge, _ = binned_statistic_2d(event_y, event_x, event_itr, bins = [bins_z, bins_z])
    print(mean_itr.shape)
    img = axes[1,0].imshow(mean_itr, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = 0.15, vmax = 0.35)
    plt.colorbar(img, ax = axes[1,0])
    axes[1,0].set_xlabel("X (mm)")
    axes[1,0].set_ylabel("Y (mm)")

    mean_multi, x_edge, y_edge, _ = binned_statistic_2d(event_y, event_x, event_multi, bins = [bins_z, bins_z])
    print(mean_itr.shape)
    img = axes[1,1].imshow(mean_multi, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = -0.120, vmax = -0.110)
    plt.colorbar(img, ax = axes[1,1])
    axes[1,1].set_xlabel("X (mm)")
    axes[1,1].set_ylabel("Y (mm)")


    fig.tight_layout()
    plt.savefig(f"../plots/low_itr_bi214/mean_itr_rho2_z_{ITR_LOW}_cut.pdf")
    plt.close()

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (12, 12))
    mean_itr, x_edge, y_edge, _ = binned_statistic_2d(mc_z, mc_rho2, mc_itr, bins = [bins_z, bins_rho])
    print(mean_itr.shape)
    img = axes[0,0].imshow(mean_itr, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = 0.15, vmax = 0.35)
    plt.colorbar(img, ax = axes[0,0])
    axes[0,0].set_xlabel("Rho (mm)")
    axes[0,0].set_ylabel("Z (mm)")
    axes[0,0].set_title("ITR")

    mean_multi, x_edge, y_edge, _ = binned_statistic_2d(mc_z, mc_rho2, mc_multi, bins = [bins_z, bins_rho])
    print(mean_itr.shape)
    img = axes[0,1].imshow(mean_multi, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = -0.120, vmax = -0.110)
    plt.colorbar(img, ax = axes[0,1])
    axes[0,1].set_xlabel("Rho (mm)")
    axes[0,1].set_ylabel("Z (mm)")
    axes[0,1].set_title("Multisite")

    mean_itr, x_edge, y_edge, _ = binned_statistic_2d(mc_y, mc_x, mc_itr, bins = [bins_z, bins_z])
    print(mean_itr.shape)
    img = axes[1,0].imshow(mean_itr, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = 0.15, vmax = 0.35)
    plt.colorbar(img, ax = axes[1,0])
    axes[1,0].set_xlabel("X (mm)")
    axes[1,0].set_ylabel("Y (mm)")

    mean_multi, x_edge, y_edge, _ = binned_statistic_2d(mc_y, mc_x, mc_multi, bins = [bins_z, bins_z])
    print(mean_itr.shape)
    img = axes[1,1].imshow(mean_multi, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto", vmin = -0.120, vmax = -0.110)
    plt.colorbar(img, ax = axes[1,1])
    axes[1,1].set_xlabel("X (mm)")
    axes[1,1].set_ylabel("Y (mm)")


    fig.tight_layout()
    plt.savefig(f"../plots/low_itr_bi214/mean_itr_rho2_z_MC_{ITR_LOW}_cut.pdf")
    plt.close()

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (12, 12))
    mean_itr_mc, x_edge, y_edge, _ = binned_statistic_2d(mc_z, mc_rho2, mc_itr, bins = [bins_z, bins_rho])
    mean_itr_data, x_edge, y_edge, _ = binned_statistic_2d(event_z, event_rho2, event_itr, bins = [bins_z, bins_rho])
    img = axes[0,0].imshow((mean_itr_data - mean_itr_mc)/mean_itr_data, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto")
    plt.colorbar(img, ax = axes[0,0])
    axes[0,0].set_xlabel("Rho (mm)")
    axes[0,0].set_ylabel("Z (mm)")
    axes[0,0].set_title("ITR: data - mc")

    mean_multi_mc, x_edge, y_edge, _ = binned_statistic_2d(mc_z, mc_rho2, mc_multi, bins = [bins_z, bins_rho])
    mean_multi_data, x_edge, y_edge, _ = binned_statistic_2d(event_z, event_rho2, event_multi, bins = [bins_z, bins_rho])
    img = axes[0,1].imshow((mean_multi_data - mean_multi_mc)/mean_multi_data, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto")
    plt.colorbar(img, ax = axes[0,1])
    axes[0,1].set_xlabel("Rho (mm)")
    axes[0,1].set_ylabel("Z (mm)")
    axes[0,1].set_title("Multisite: data - mc")

    mean_itr_mc, x_edge, y_edge, _ = binned_statistic_2d(mc_y, mc_x, mc_itr, bins = [bins_z, bins_z])
    mean_itr_data, x_edge, y_edge, _ = binned_statistic_2d(event_y, event_x, event_itr, bins = [bins_z, bins_z])
    img = axes[1,0].imshow((mean_itr_data - mean_itr_mc)/mean_itr_data, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto")
    plt.colorbar(img, ax = axes[1,0])
    axes[1,0].set_xlabel("X (mm)")
    axes[1,0].set_ylabel("Y (mm)")

    mean_multi_mc, x_edge, y_edge, _ = binned_statistic_2d(mc_y, mc_x, mc_multi, bins = [bins_z, bins_z])
    mean_multi_data, x_edge, y_edge, _ = binned_statistic_2d(event_y, event_x, event_multi, bins = [bins_z, bins_z])
    print(mean_itr.shape)
    img = axes[1,1].imshow((mean_multi_data - mean_multi_mc)/mean_multi_data, extent = [y_edge[0], y_edge[-1], x_edge[0], x_edge[-1]], aspect = "auto")
    plt.colorbar(img, ax = axes[1,1])
    axes[1,1].set_xlabel("X (mm)")
    axes[1,1].set_ylabel("Y (mm)")


    fig.tight_layout()
    plt.savefig(f"../plots/low_itr_bi214/mc_data_diffs_{ITR_LOW}_cut.pdf")
    plt.close()

    plt.figure()
    bins_z = np.arange(3000, 6000, 500)
    mc_average_binned, bins, _ = binned_statistic(mc_r, mc_itr, bins = bins_z)
    mc_std_binned, bins, _ = binned_statistic(mc_r, mc_itr, bins = bins_z, statistic="std")

    data_average_binned, bins, _ = binned_statistic(event_r, event_itr, bins = bins_z)
    data_std_binned, bins, _ = binned_statistic(event_r, event_itr, bins = bins_z, statistic="std")

    plt.errorbar(bins_z[:-1] + np.diff(bins_z)[0]/2, mc_average_binned, yerr = mc_std_binned, xerr = np.diff(bins_z) / 2, linestyle = "", color = "black")
    plt.errorbar(bins_z[:-1] + np.diff(bins_z)[0]/2, data_average_binned, yerr = data_std_binned, xerr = np.diff(bins_z) / 2, linestyle = "", color = "red")
    plt.xlabel("Reconstructed Radius (mm)")
    plt.ylabel("ITR")
    plt.savefig(f"../plots/low_itr_bi214/mc_data_r_vs_itr_scatter_{ITR_LOW}cut.pdf")
    plt.close()

    plt.figure()
    print(np.isnan(mc_r).any())
    print(np.isnan(mc_multi).any())
    mc_average_binned, bins, _ = binned_statistic(mc_r, mc_multi, bins = bins_z)
    print(mc_average_binned)
    mc_std_binned, bins, _ = binned_statistic(mc_r, mc_multi, bins = bins_z, statistic="std")

    data_average_binned, bins, _ = binned_statistic(event_r, event_multi, bins = bins_z)
    print(data_average_binned)
    data_std_binned, bins, _ = binned_statistic(event_r, event_multi, bins = bins_z, statistic="std")

    plt.errorbar(bins_z[:-1] + np.diff(bins_z)[0]/2, mc_average_binned, yerr = mc_std_binned, xerr = np.diff(bins_z) / 2, linestyle = "", color = "red", label = "MC")
    plt.errorbar(bins_z[:-1] + np.diff(bins_z)[0]/2, data_average_binned, yerr = data_std_binned, xerr = np.diff(bins_z) / 2, alpha = 0.7, linestyle = "", color = "black", label = "Data")
    plt.legend()
    plt.xlabel("Reconstructed Radius (mm)")
    plt.ylabel("Multisite")
    plt.ylim((-0.123, -0.11))
    plt.savefig(f"../plots/low_itr_bi214/mc_data_r_vs_multi_scatter_{ITR_LOW}cut.pdf")
    plt.close()
    
    plt.figure()
    perc_error_per_FV = (data_average_binned - mc_average_binned) / data_average_binned
    plt.errorbar(bins_z[:-1]+ np.diff(bins_z)[0]/2, perc_error_per_FV, xerr = np.diff(bins_z)/2, capsize = 2, marker = "o", linestyle = "", color = "black")
    plt.xlabel("Reconstructed Radius (mm)")
    plt.ylabel("Mean % Error in Multisite (Data - MC)/data")
    plt.ylim((0, 0.01))
    plt.tight_layout()
    plt.savefig("../plots/low_itr_bi214/presentation/perc_error_with_FV.pdf")
    plt.close()

    plt.figure()
    plt.hist(removed_itr_r, bins = np.linspace(0, 6000,40), density = False, histtype = "step")
    plt.axvline(np.mean(removed_itr_r), label = f"mean: {np.mean(removed_itr_r):.3f} +- {np.std(removed_itr_r):.3f} mm", color = "black")
    plt.axvline(np.mean(removed_itr_r) - np.std(removed_itr_r), color = "red")
    plt.axvline(np.mean(removed_itr_r) + np.std(removed_itr_r), color = "red")
    plt.legend()
    plt.xlabel("Radius (mm)")
    plt.ylabel("Counts")
    plt.title(f"Data event Radius Removed by ITR Cut: ITR < {ITR_LOW}")
    plt.savefig(f"../plots/low_itr_bi214/data_radius_removed_ITR_{ITR_LOW}cut.pdf")

    plt.figure()
    plt.hist(removed_itr_x, bins = np.linspace(-6000, 6000,40), density = False, histtype = "step")
    plt.axvline(np.mean(removed_itr_x), label = f"mean: {np.mean(removed_itr_x):.3f} +- {np.std(removed_itr_x):.3f} mm", color = "black")
    plt.axvline(np.mean(removed_itr_x) - np.std(removed_itr_x), color = "red")
    plt.axvline(np.mean(removed_itr_x) + np.std(removed_itr_x), color = "red")
    plt.legend()
    plt.xlabel("X (mm)")
    plt.ylabel("Counts")
    plt.title(f"Data event Radius Removed by ITR Cut: ITR < {ITR_LOW}")
    plt.savefig(f"../plots/low_itr_bi214/data_X_removed_ITR_{ITR_LOW}cut.pdf")

    plt.figure()
    plt.hist(removed_itr_y, bins = np.linspace(-6000, 6000,40), density = False, histtype = "step")
    plt.axvline(np.mean(removed_itr_y), label = f"mean: {np.mean(removed_itr_y):.3f} +- {np.std(removed_itr_y):.3f} mm", color = "black")
    plt.axvline(np.mean(removed_itr_y) - np.std(removed_itr_y), color = "red")
    plt.axvline(np.mean(removed_itr_y) + np.std(removed_itr_y), color = "red")
    plt.legend()
    plt.xlabel("Y (mm)")
    plt.ylabel("Counts")
    plt.title(f"Data event Radius Removed by ITR Cut: ITR < {ITR_LOW}")
    plt.savefig(f"../plots/low_itr_bi214/data_Y_removed_ITR_{ITR_LOW}cut.pdf")

    plt.figure()
    plt.hist(removed_itr_z, bins = np.linspace(-6000, 6000,40), density = False, histtype = "step")
    plt.axvline(np.mean(removed_itr_z), label = f"mean: {np.mean(removed_itr_z):.3f} +- {np.std(removed_itr_z):.3f} mm", color = "black")
    plt.axvline(np.mean(removed_itr_z) - np.std(removed_itr_r), color = "red")
    plt.axvline(np.mean(removed_itr_z) + np.std(removed_itr_r), color = "red")
    plt.legend()
    plt.xlabel("Z (mm)")
    plt.ylabel("Counts")
    plt.title(f"Data event Radius Removed by ITR Cut: ITR < {ITR_LOW}")
    plt.savefig(f"../plots/low_itr_bi214/data_Z_removed_ITR_{ITR_LOW}cut.pdf")

    # fig = px.density_heatmap(mean_itr, text_auto = True)
    # fig.show()
    

def check_pdfs():

    file = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_pdf/full_analysis2_B8_solar_nue/300005.root")
    branch = file.Get("PDF_event_by_event_Information")

    x = []
    y = []
    z = []
    r = []
    for ientry in branch:
        x = ientry.x
        y = ientry.y
        z = ientry.z

        radius = np.sqrt(x**2 + y**2 + z**2)
        r.append(radius)

    plt.hist(r, bins = np.arange(0, 6000, 10), histtype = "step")
    plt.savefig("../plots/low_itr_bi214/pdf_radius_check_fullanalysis2.png")

def simple_multisite_comparison():
    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    # file = ROOT.TFile.Open(working_dir + "/extracted_data/bi214/bismuth214_data_discriminants/total.root")
    # branch = file.Get("1p25_3p0")
    file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/total.root")
    branch = file.Get("bi214")
    
    FV_CUT   = 3000
    ITR_LOW  = 0.18
    ITR_HIGH = 0.3
    data_multisite = []
    for ientry in branch:

        if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
            continue
        if ientry.itr < 0.18:
            print(ientry.itr)
        x = ientry.x
        y = ientry.y
        z = ientry.z
        r = np.sqrt(x**2 + y**2 + z**2)
        if r > FV_CUT:
            continue

        data_multisite.append(ientry.dlogL)
    
    # plot the MC versions
    mc       = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/Bi214/total.root")
    mc_tree  = mc.Get("1p25_3p0")
    
    mc_multi = []    
    for ientry in mc_tree:

        if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
            continue
        itr       = ientry.itr
        multi     = ientry.dlogL
        x         = ientry.x
        y         = ientry.y
        z         = ientry.z
        r         = np.sqrt(x**2 + y**2 + z**2)
        if r > 5250:
            continue
        
        mc_multi.append(multi)
    
    
    binning = np.arange(-0.130, -0.105, 0.0005)
    binning = np.linspace(-0.125, -0.11, 20)
    mids = binning[:-1] + np.diff(binning)[0] / 2 
    counts_data, bins_data = np.histogram(data_multisite, bins = binning)
    counts_mc, bins_mc     = np.histogram(mc_multi,  bins = binning)
    
    # scale the counts in the MC to match the counts in the data
    scale = np.sum(counts_data) / np.sum(counts_mc)
    error_mc = np.sqrt(counts_mc) * scale
    error_data = np.sqrt(counts_data)

    plt.figure()
    plt.errorbar(mids, counts_data, yerr = error_data, linestyle = "", color = "black", capsize = 2, label = f"data | % Difference: {(np.mean(data_multisite)-np.mean(mc_multi))/np.mean(data_multisite):.3f}")
    plt.hist(mc_multi, color = "red", bins = binning, weights = np.full_like(mc_multi, scale), histtype="step", label = r"$^{214}Bi$ MC")
    plt.xlabel("Multisite")
    plt.title(f"FV: {FV_CUT} mm | {ITR_LOW} < ITR < {ITR_HIGH}")
    plt.legend()
    plt.savefig(f"../plots/low_itr_bi214/presentation/with_FV/multisite_agreement_7.0.15_{ITR_LOW}_{ITR_HIGH}ITR_{FV_CUT}FV.pdf")
    plt.close()

def excluding_old_data():
    
    run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    
    FV_CUT   = 4000
    ITR_LOW  = 0.18
    ITR_HIGH = 0.3
    data_multisite = []
    for irun in run_list:
        if irun < 301122:
            continue
        file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_PDF/{irun}.root")
        branch = file.Get("bi214")
        for ientry in branch:

            if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
                continue
            if ientry.itr < 0.18:
                print(ientry.itr)
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue

            data_multisite.append(ientry.dlogL)
    
    # plot the MC versions
    mc       = ROOT.TFile.Open("/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test/Bi214/total.root")
    mc_tree  = mc.Get("1p25_3p0")
    
    mc_multi = []    
    for ientry in mc_tree:

        if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
            continue
        itr       = ientry.itr
        multi     = ientry.dlogL
        x         = ientry.x
        y         = ientry.y
        z         = ientry.z
        r         = np.sqrt(x**2 + y**2 + z**2)
        if r > 5250:
            continue
        
        mc_multi.append(multi)
    
    
    binning = np.arange(-0.130, -0.105, 0.0005)
    binning = np.linspace(-0.125, -0.11, 20)
    mids = binning[:-1] + np.diff(binning)[0] / 2 
    counts_data, bins_data = np.histogram(data_multisite, bins = binning)
    counts_mc, bins_mc     = np.histogram(mc_multi,  bins = binning)
    
    # scale the counts in the MC to match the counts in the data
    scale = np.sum(counts_data) / np.sum(counts_mc)
    error_mc = np.sqrt(counts_mc) * scale
    error_data = np.sqrt(counts_data)

    plt.figure()
    plt.errorbar(mids, counts_data, yerr = error_data, marker = "o", markersize = 2, linestyle = "", color = "black", capsize = 2, label = f"data | % Difference: {(np.mean(data_multisite)-np.mean(mc_multi))/np.mean(data_multisite):.3f}")
    plt.hist(mc_multi, color = "red", bins = binning, weights = np.full_like(mc_multi, scale), histtype="step", label = r"$^{214}Bi$ MC")
    plt.xlabel("Multisite")
    plt.title(f"FV: {FV_CUT} mm | {ITR_LOW} < ITR < {ITR_HIGH}")
    plt.legend()
    plt.savefig(f"../plots/low_itr_bi214/presentation/excluding_early_data/multisite_agreement_7.0.15_{ITR_LOW}_{ITR_HIGH}ITR_{FV_CUT}FV.pdf")
    plt.close()

def plot_time_residual_agreement(mc_residuals):

    run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    working_dir     = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies"
    
    
    FV_CUT   = 5250
    ITR_LOW  = 0.18
    ITR_HIGH = 0.3
    log_flag = False
    data_tres       = []
    quiet_data_tres = []
    dirty_data_tres = []
    for irun in run_list:
        print(irun)
        file = ROOT.TFile.Open(working_dir + f"/extracted_data/bi214/reprocessed_7.0.15_discriminants_2.5_3.125_tres/{irun}.root")
        branch = file.Get("hit_level")
        for ientry in branch:

            if ientry.itr < ITR_LOW or ientry.itr > ITR_HIGH:
                continue
            if ientry.itr < 0.18:
                print(ientry.itr)
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            tres = list(ientry.tRes)
            if irun < 301122:
                dirty_data_tres += tres
            else:
                quiet_data_tres += tres
            data_tres += tres

    
    binning = np.arange(-10, 81, 1)
    plt.figure()
    plt.hist(dirty_data_tres, bins = binning, density = True, histtype = "step", label = "Dirty", color = "black")
    plt.hist(quiet_data_tres, bins = binning, density = True, histtype = "step", label = "Clean", color = "blue")
    plt.hist(mc_residuals, bins = binning, density = True, histtype = "step", label = "MC", color = "red")

    plt.xlabel("Time Residual (ns)")
    plt.ylabel("Counts")
    plt.legend()
    if log_flag == True:
        plt.yscale("log")
    plt.title(rf"Time Residuals | FV {FV_CUT} m | ${ITR_LOW} \leq ITR \leq {ITR_HIGH}$ | ${1.25} \leq E \leq {3.00}$ MeV")
    plt.savefig(f"../plots/low_itr_bi214/presentation/time_residuals/dirty_vs_clean_{FV_CUT}_FV_log_{log_flag}_peakBins.pdf")
    plt.close()

def extract_mc_residuals():
    
    run_list = np.loadtxt("../runlists/contains_solar_candidates.txt", dtype = int)
    input_fpath  = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ratds/simulationBiPo214"
    E_LOW = 1.25
    E_HIGH = 3.00
    FV_CUT = 5250
    ITR_LOW = 0.18
    ITR_HIGH = 0.3
    residualsRECON= []
    for irun in [1]:
        for ientry, _ in rat.dsreader(f"{input_fpath}_307919.root"):
            # light path calculator and point3D stuff loaded after ratds constructor
            # timeResCalc = rat.utility().GetTimeResidualCalculator()
            PMTCalStatus = RAT.DU.Utility.Get().GetPMTCalStatus()
            light_path = rat.utility().GetLightPathCalculator()
            group_velocity = rat.utility().GetGroupVelocity()
            pmt_info = rat.utility().GetPMTInfo()
            psup_system_id = RAT.DU.Point3D.GetSystemId("innerPMT")
            av_system_id = RAT.DU.Point3D.GetSystemId("av")
            
            # entry = ds.GetEntry(i)
            if ientry.GetEVCount() == 0:
                continue

            #### RECONSTRUCTION INFORMATION EXTRACTED ####
            reconEvent = ientry.GetEV(0)
            
            # did event get reconstructed correctly?
            fit_name = reconEvent.GetDefaultFitName()
            if not reconEvent.FitResultExists(fit_name):
                continue

            vertex = reconEvent.GetFitResult(fit_name).GetVertex(0)
            if (not vertex.ContainsPosition() or
                not vertex.ContainsTime() or
                not vertex.ValidPosition() or
                not vertex.ValidTime() or
                not vertex.ContainsEnergy() or
                not vertex.ValidEnergy()):
                continue
            # print("Reconstruction checks PASSED!")
            # reconstruction valid so get reconstructed position and energy
            reconPosition  = vertex.GetPosition() # returns in PSUP coordinates
            reconEnergy    = vertex.GetEnergy()        
            reconEventTime = vertex.GetTime()
            
            # apply AV offset to position
            event_point = RAT.DU.Point3D(psup_system_id, reconPosition)
            event_point.SetCoordinateSystem(av_system_id)
            if event_point.Mag() > FV_CUT:
                continue
            # convert back to PSUP coordinates
            event_point.SetCoordinateSystem(psup_system_id)
            # apply energy tagging cuts the same as that in data
            if reconEnergy < E_LOW or reconEnergy > E_HIGH:
                continue
            
            itr = reconEvent.GetClassifierResult("ITR:scintFitter").GetClassification("ITR");
            if itr < ITR_LOW or itr > ITR_HIGH:
                continue
            # event has passed all the cuts so we can extract the time residuals
            calibratedPMTs = reconEvent.GetCalPMTs()
            pmtCalStatus = rat.utility().GetPMTCalStatus()
            for j in range(calibratedPMTs.GetCount()):
                pmt = calibratedPMTs.GetPMT(j)
                if pmtCalStatus.GetHitStatus(pmt) != 0:
                    continue
                
                # residual_recon = timeResCalc.CalcTimeResidual(pmt, reconPosition, reconEventTime, True)
                pmt_point = RAT.DU.Point3D(psup_system_id, pmt_info.GetPosition(pmt.GetID()))
                light_path.CalcByPosition(event_point, pmt_point)
                inner_av_distance = light_path.GetDistInInnerAV()
                av_distance = light_path.GetDistInAV()
                water_distance = light_path.GetDistInWater()
                transit_time = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)
                residual_recon = pmt.GetTime() - transit_time - reconEventTime
                
                residualsRECON.append(residual_recon)
            
    
    return residualsRECON

# simple_multisite_comparison()
# plot_itr_vs_multi()
# plot_time_residuals()
# write_gtids_of_weird_events()
# check_pdfs()
# excluding_old_data()
mc_residuals = extract_mc_residuals()
plot_time_residual_agreement(mc_residuals)

    
    

