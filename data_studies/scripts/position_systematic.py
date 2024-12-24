import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
Script is to evaluate the position systematic on tagged BiPo214. Will use my
16000 odd sample of tagged BiPo214 used for the MC vs Data multisite agreement 
test.

Comparison of the dR distributions between MC and data. Applying a gaussian 
smearing factor on the data of X mm, and finding which smear gives the best
data / MC agreement.

re-run analysis (I guess BiPo tagging? That will take ages ...) with +- that FV
and propagate the expectation value...
"""

def get_positions(tree, bi_x, bi_y, bi_z, po_x, po_y, po_z, bi_energy, po_energy):
    print("Total entries: ", tree.GetEntries())
    for ientry in tree:
        x_bi = ientry.x_bi
        y_bi = ientry.y_bi
        z_bi = ientry.z_bi

        x_po = ientry.x_po
        y_po = ientry.y_po
        z_po = ientry.z_po

        # check if the Bi214 event reconstructs inside the FV
        r_bi = np.sqrt(x_bi**2 + y_bi**2 + z_bi**2)
        r_po = np.sqrt(x_po**2 + y_po**2 + z_po**2)
        if r_bi > FV_CUT and r_po > FV_CUT:
            continue

        bi_x.append(x_bi)
        bi_y.append(y_bi)
        bi_z.append(z_bi)
        po_x.append(x_po)
        po_y.append(y_po)
        po_z.append(z_po)
        bi_energy.append(ientry.energy_bi)
        po_energy.append(ientry.energy_po)

def systematic_scan():
    # load in the MC dR
    directory = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/better_tagging"
    file_data = ROOT.TFile.Open(f"{directory}/raw_tagging_output/total.root")
    file_mc   = ROOT.TFile.Open(f"/data/snoplus3/hunt-stokes/nuPhysPoster/mc_tagged_bipo/214/total.root")

    tree_data = file_data.Get("tag_info")
    tree_mc   = file_mc.Get("tag_info")

    FV_CUT = 6000

    data_pos_bi_x = []
    data_pos_bi_y = []
    data_pos_bi_z = []
    data_pos_po_x = []
    data_pos_po_y = []
    data_pos_po_z = []
    data_bi_e     = []
    data_po_e     = []

    mc_pos_bi_x   = []
    mc_pos_bi_y   = []
    mc_pos_bi_z   = []
    mc_pos_po_x   = []
    mc_pos_po_y   = []
    mc_pos_po_z   = []
    mc_bi_e       = []
    mc_po_e       = []

    get_positions(tree_data, data_pos_bi_x, data_pos_bi_y, data_pos_bi_z, data_pos_po_x, data_pos_po_y, data_pos_po_z, data_bi_e, data_po_e)
    get_positions(tree_mc, mc_pos_bi_x, mc_pos_bi_y, mc_pos_bi_z, mc_pos_po_x, mc_pos_po_y, mc_pos_po_z, mc_bi_e, mc_po_e)

    # np.save("./data_positions_bi.npy", data_positions_bi)
    # np.save("./data_positions_po.npy", data_positions_po)
    # np.save("./mc_positions_bi.npy", mc_positions_bi)
    # np.save("./mc_positions_po.npy", mc_positions_po)

    # data_positions_bi = np.load("./data_positions_bi.npy")
    # data_positions_po = np.load("./data_positions_po.npy")
    # mc_positions_bi   = np.load("./mc_positions_bi.npy")
    # mc_positions_po   = np.load("./mc_positions_po.npy")

    num_data = len(data_pos_bi_x)
    num_mc   = len(mc_pos_bi_x)

    mc_positions   = np.zeros((num_mc, 6))
    data_positions = np.zeros((num_data, 6)) 

    mc_positions[:,0] = mc_pos_bi_x
    mc_positions[:,1] = mc_pos_bi_y
    mc_positions[:,2] = mc_pos_bi_z
    mc_positions[:,3] = mc_pos_po_x
    mc_positions[:,4] = mc_pos_po_y
    mc_positions[:,5] = mc_pos_po_z

    data_positions[:,0] = data_pos_bi_x
    data_positions[:,1] = data_pos_bi_y
    data_positions[:,2] = data_pos_bi_z
    data_positions[:,3] = data_pos_po_x
    data_positions[:,4] = data_pos_po_y
    data_positions[:,5] = data_pos_po_z

    print("MC Events: ", mc_positions.shape)
    print("Data Events: ", data_positions.shape)

    """
    Now we apply a gaussian smearing to every bi and every po MC position following gaussian of given width. We assume the same
    smear in each direction (although this could be improved upon).

    Following application of each smear, we recalculate the dR distribution in MC and compare to data.

    The chi2 between the dR MC and data as calculated. We return the smear that gives the closest match to data dR distribution.
    """

    standard_deviations = np.arange(1, 101, 1) # in mm
    bin_width           = 30
    dr_bins             = np.arange(0, 2000 + bin_width, bin_width)
    chi2                = []
    data_dr             = np.sqrt((data_positions[:, 0] - data_positions[:, 3])**2 + (data_positions[:, 1] - data_positions[:, 4])**2 + (data_positions[:, 2] - data_positions[:, 5])**2)

    for istd in standard_deviations:

        # do not edit / smear the original positions
        copy_mc_positions = np.copy(mc_positions)

        # generate the smear factors for each Bi-Po pair in the MC
        gaussian_smears   = np.random.normal(loc = 0, scale = istd, size = mc_positions.shape)

        # add the smears to the mc positions x y z for Bi and Po
        copy_mc_positions = copy_mc_positions + gaussian_smears

        # calculate the dR distribution from this smeared distribution
        smeared_dr = np.sqrt((copy_mc_positions[:, 0] - copy_mc_positions[:, 3])**2 + (copy_mc_positions[:, 1] - copy_mc_positions[:, 4])**2 + (copy_mc_positions[:, 2] - copy_mc_positions[:, 5])**2)

        # work out the uncertainty on the model and data for the chi2
        counts_mc, _   = np.histogram(smeared_dr, bins = dr_bins)
        counts_data, _ = np.histogram(data_dr, bins = dr_bins)

        # work out the normalisation factor (integral / area of the histogram)
        norm_mc   = np.sum(bin_width * counts_mc)
        norm_data = np.sum(bin_width * counts_data)
        err_mc    = np.sqrt(counts_mc)   / norm_mc
        err_data  = np.sqrt(counts_data) / norm_data
        print(err_data)
        # tot_err   = np.sqrt(err_mc**2 + err_data**2)
        tot_err   = err_data

        # normalise the counts
        counts_data = counts_data / norm_data
        counts_mc   = counts_mc   / norm_mc

        idx_nonzero = np.nonzero(tot_err)
        # find the sum of the chi2 for this smearing
        chi2_val = np.sum((counts_data[idx_nonzero] - counts_mc[idx_nonzero])**2 / tot_err[idx_nonzero]**2)
        # normalise by NDF = N.bins - N_params - 1
        # print("num bins: ", len(idx_nonzero[0]))
        # chi2_val /= (len(idx_nonzero[0]) - 1 - 1)
        
        chi2.append(chi2_val)

        # create a plot of the data vs MC dR distribution
        plt.hist(smeared_dr, density = True, bins = dr_bins, histtype = "step", color = "red", linewidth = 2, label = f"MC | Sys. {istd} mm")
        plt.hist(data_dr, density = True, bins = dr_bins, histtype = "step", color = "black", linewidth = 2, label = f"Data")
        plt.xlim((0, 1000))
        plt.legend()
        plt.xlabel(r"$\Delta R$ [mm]")
        plt.xlabel("Normalised Counts")
        plt.savefig(f"../plots/position_systematic_plots/{istd}.pdf")
        plt.close()

    # fit the systematics to a quadratic
    def quadradtic(x, A, B, C, D):
        return A*x**3 + B*x**2 + C*x + D

    def expo(x, A, B, C):
        return A*np.exp(x*B) + C
    plt.figure()

    popt, pcov = curve_fit(quadradtic, standard_deviations, chi2)
    fitted_x = np.linspace(0, 100, 1000)
    fitted_y = quadradtic(fitted_x, *popt)
    plt.plot(fitted_x, fitted_y, linewidth = 1, color = "red", alpha = 0.5, label = f"Cubices Fit")

    popt, pcov = curve_fit(expo, standard_deviations, chi2, p0 = [1, 0.05, 50])
    fitted_x = np.linspace(0, 100, 1000)
    fitted_y = expo(fitted_x, *popt)
    plt.plot(fitted_x, fitted_y, linewidth = 1, color = "blue", label = f"Exponential Fit", alpha = 0.5)

    # plot the chi2 per smear factor
    plt.scatter(standard_deviations, chi2, s = 1, color = "black")
    plt.xlabel(r"$\sigma _{sys}$ [mm]")
    plt.ylabel(r"$\sum_i \frac{(D_i - MC_i)^2}{\sigma_i^2}$")
    plt.legend(frameon = False)
    # plt.yscale("log")
    plt.savefig("../plots/position_systematic_plots/chi2_2s.pdf")
    plt.close()

    ### energy systematic ###
    """
    This section of the code evaluates the energy systematics. The shift in mean and median of
    Po214 and Bi214 energy distributions is calculated, respectively.

    Then, accounting for this shift, the distributions are shifted so their centred on 0 and 
    we decided if we need an energy spearing factor as well.
    """

    data_mean_po   = np.mean(data_po_e)
    mc_mean_po     = np.mean(mc_po_e)
    data_median_bi = np.median(data_bi_e)
    mc_median_bi   = np.median(mc_bi_e)

    ratio_po = data_mean_po/mc_mean_po
    ratio_bi = data_median_bi/mc_median_bi
    print(f"Po214 data / MC mean: {data_mean_po/mc_mean_po}")
    print(f"Bi214 data / MC median: {data_median_bi/mc_median_bi}")
    po_bins = np.arange(-5, 5.1, 0.025)
    bi_bins = np.arange(-5, 5.00, 0.1)
    plt.hist(data_po_e, bins = po_bins, histtype = "step", density = True, label = "Po214 Data", color = "black")
    plt.hist(mc_po_e, bins = po_bins, histtype = "step", density = "True", label = "Po214 MC", color = "red")
    plt.xlabel("Reconstructed Energy [MeV]")
    plt.legend()
    plt.xlim((0, 2.0))
    plt.savefig("../plots/position_systematic_plots/po214_energy.pdf")
    plt.close()

    plt.hist(data_bi_e, bins = bi_bins,  histtype = "step", density = True, label = "Bi214 Data", color = "black")
    plt.hist(mc_bi_e, bins = bi_bins, histtype = "step", density = "True", label = "Bi214 MC", color = "red")
    plt.xlabel("Reconstructed Energy [MeV]")
    plt.xlim((0, 5.0))
    plt.legend()
    plt.savefig("../plots/position_systematic_plots/bi214_energy.pdf")
    plt.close()

    # now shift all the MC by the ratios calculated and check the widths
    mc_po_e = np.array(mc_po_e) * ratio_po
    plt.hist(data_po_e, bins = po_bins, histtype = "step", density = True, label = "Po214 Data", color = "black")
    plt.hist(mc_po_e, bins = po_bins, histtype = "step", density = "True", label = "Po214 MC", color = "red")
    plt.xlabel("Reconstructed Energy [MeV]")
    plt.savefig("../plots/position_systematic_plots/po214_energy_shifted.pdf")
    plt.close()

    mc_bi_e = np.array(mc_bi_e) * ratio_bi

    plt.hist(data_bi_e, bins = bi_bins,  histtype = "step", density = True, label = "Bi214 Data", color = "black")
    plt.hist(mc_bi_e, bins = bi_bins, histtype = "step", density = "True", label = "Bi214 MC", color = "red")
    plt.xlabel("Reconstructed Energy [MeV]")
    plt.savefig("../plots/position_systematic_plots/bi214_energy_shifted.pdf")
    plt.close()

    data_mean_po   = np.mean(data_po_e)
    mc_mean_po     = np.mean(mc_po_e)
    data_median_bi = np.median(data_bi_e)
    mc_median_bi   = np.median(mc_bi_e)
    ratio_po = data_mean_po/mc_mean_po
    ratio_bi = data_median_bi/mc_median_bi

    print("Shifted Ratios: ", ratio_po, ratio_bi)

    # now we apply gaussian smearing to the reconstructed energy of each one, and find the chi2 of each smearing
    # to find a width systematic on the energy spectra
    standard_deviations = np.arange(0.01, 0.25, 0.01) # in MeV
    chi2_po             = []
    chi2_bi             = []
    for istd in standard_deviations:

        # do not edit / smear the original energies
        copy_mc_po_e = np.copy(mc_po_e)
        copy_mc_bi_e = np.copy(mc_bi_e)

        # generate the smear factors for each Bi-Po pair in the MC
        gaussian_smears   = np.random.normal(loc = 0, scale = istd, size = copy_mc_bi_e.shape)

        # add the smears to the mc energies for both Bi and Po
        copy_mc_po_e = copy_mc_po_e + gaussian_smears
        copy_mc_bi_e = copy_mc_bi_e + gaussian_smears

        # work out the uncertainty on the model and data for the chi2
        counts_mc_po, _   = np.histogram(copy_mc_po_e, bins = po_bins)
        counts_mc_bi, _   = np.histogram(copy_mc_bi_e, bins = bi_bins)

        counts_data_po, _ = np.histogram(data_po_e, bins = po_bins)
        counts_data_bi, _ = np.histogram(data_bi_e, bins = bi_bins)

        # work out the normalisation factor (integral / area of the histogram)
        norm_mc_po   = np.sum(bin_width * counts_mc_po)
        norm_mc_bi   = np.sum(bin_width * counts_mc_bi)
        norm_data_po = np.sum(bin_width * counts_data_po)
        norm_data_bi = np.sum(bin_width * counts_data_bi)
        err_mc_po    = np.sqrt(counts_mc_po)   / norm_mc_po
        err_mc_bi    = np.sqrt(counts_mc_bi)   / norm_mc_bi
        err_data_po  = np.sqrt(counts_data_po) / norm_data_po
        err_data_bi  = np.sqrt(counts_data_bi) / norm_data_bi
        tot_err_po   = err_data_po#np.sqrt(err_mc_po**2 + err_data_po**2)
        tot_err_bi   = err_data_bi#np.sqrt(err_mc_bi**2 + err_data_bi**2)
        # print(tot_err_po, tot_err_bi, istd)
        # normalise the counts
        counts_data_po = counts_data_po / norm_data_po
        counts_data_bi = counts_data_bi / norm_data_bi
        counts_mc_po   = counts_mc_po   / norm_mc_po
        counts_mc_bi   = counts_mc_bi   / norm_mc_bi

        nonzero_po = np.nonzero(tot_err_po)
        nonzero_bi = np.nonzero(tot_err_bi)
        # find the sum of the chi2 for this smearing
        chi2_val_po = np.sum((counts_data_po[nonzero_po] - counts_mc_po[nonzero_po])**2 / tot_err_po[nonzero_po]**2)
        chi2_val_bi = np.sum((counts_data_bi[nonzero_bi] - counts_mc_bi[nonzero_bi])**2 / tot_err_bi[nonzero_bi]**2)
        chi2_po.append(chi2_val_po)
        chi2_bi.append(chi2_val_bi)

        # create a plot of the data vs MC dR distribution
        plt.hist(copy_mc_po_e, density = True, bins = po_bins, histtype = "step", color = "red", linewidth = 2, label = f"MC | Sys. {istd} mm")
        plt.hist(data_po_e, density = True, bins = po_bins, histtype = "step", color = "black", linewidth = 2, label = f"Data")
        
        plt.legend()
        plt.xlabel(r"E [MeV]")
        plt.ylabel("Normalised Counts")
        plt.savefig(f"../plots/position_systematic_plots/{istd}_energy.pdf")
        plt.close()

    # fit the systematics to a quadratic
    def quadradtic(x, A, B, C, D):
        return A*x**3 + B*x**2 + C*x + D

    def expo(x, A, B, C):
        return A*np.exp(x*B) + C
    plt.figure()

    # popt, pcov = curve_fit(quadradtic, standard_deviations, chi2)
    # fitted_x = np.linspace(0, 100, 1000)
    # fitted_y = quadradtic(fitted_x, *popt)
    # plt.plot(fitted_x, fitted_y, linewidth = 1, color = "red", alpha = 0.5, label = f"Cubices Fit")

    # popt, pcov = curve_fit(expo, standard_deviations, chi2, p0 = [1, 0.05, 50])
    # fitted_x = np.linspace(0, 100, 1000)
    # fitted_y = expo(fitted_x, *popt)
    # plt.plot(fitted_x, fitted_y, linewidth = 1, color = "blue", label = f"Exponential Fit", alpha = 0.5)

    # plot the chi2 per smear factor
    # plt.scatter(standard_deviations, chi2_po, s = 2, color = "red", label = "Po214")
    plt.scatter(standard_deviations, chi2_bi, s = 2, color = "blue", label = "Bi214")
    plt.xlabel(r"$\sigma _{sys}$ [MeV]")
    plt.ylabel(r"$\sum_i \frac{(D_i - MC_i)^2}{\sigma_i^2}$")
    plt.legend(frameon = False)
    # plt.yscale("log")
    plt.tight_layout()
    plt.savefig("../plots/position_systematic_plots/chi2_energy_bi.pdf")
    plt.close()

def build_pdfs_position_shifted(FV_CUT, systematic):
    """
    Function is used to create energy and multisite PDFs that are built from events
    within +- position systematic FV.

    Ignoring the position dependence of time residuals for now.
    """

    # load up the PDF run list
    