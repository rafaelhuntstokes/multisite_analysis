import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT

"""
This script is the *correct* way to find the improved sensitivity to B8 flux with an 
Asimov dataset, created from *real* MC.

1. B8 and Tl208 MC is loaded and their multisite discriminants and energy saved to arrays.
2. 2D histograms of M vs E for signal and background are created and normalised such that 
   the sum of counts = 1.
3. An Asimov dataset is created by scaling the normalised PDFs to a desired event number (
   and ratio). The scaled PDFs are added together to create the Asimov dataset.
4. The log-likelihood function is calculated for a set of different signal hypotheses.
    - chi2
    - Poisson only
    - Multisite only
    - Multisite + Poisson
5. Frequentist confidence intervals are plotted to show the minimum & uncertainty in each 
   energy bin.
"""
def loop_over_entries(ntuple, fv_cut):
    """
    Simple loop function over entries in TTree to avoid duplicating code. Extracts the 
    energy and multisite discriminant value for the given ntuple and FV cut.
    """

    energy = []
    dlogL  = []
    for ientry in ntuple:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            e = ientry.energy
            L = ientry.dlogL
            
            # apply FV cut
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > fv_cut:
                continue

            # energy cut already applied on the branch 2.5 --> 5 MeV
            energy.append(e)
            dlogL.append(L)
    
    return energy, dlogL

def create_pdfs_and_asimov(fv_cut, energy_pdf_bins, multi_pdf_bins, num_sig, num_back):
    """
    Function loads the Tl208 and B8 MC - this MC is different to that which produced the 
    time residual PDFs (that in turn are used to produce the multisite discriminant).
    
    Normalised PDFs in energy and M are created for signal (B8) and background (Tl208).

    An Asimov dataset is created by scaling the normalised PDFs for a given signal and 
    background normalisation and adding the two distrubtions together.

    The PDFs and the Asimov dataset are plotted as 2D histograms, and the arrays returned. 
    """

    # load up the mc run list
    runlist = np.loadtxt("../runlists/full_test.txt", dtype = int)

    # arrays to store the PDF energy and multisite discriminant values
    energy_b8    = []
    energy_tl208 = []
    multi_b8     = []
    multi_tl208  = []
    failedB8     = 0   # count the number of ntuples that do not exist
    failedTl208  = 0

    # loop over every run in the run list and extract the multsite discriminant and energy
    # within the FV
    for irun in runlist:

        # some runs exist for B8 or Tl208 only, so have flags to check if the files exist 
        # for each isotope of the given run
        exists_tl = True
        exists_b  = True
        try:
            file_b8    = ROOT.TFile.Open(f"../run_by_run_test/B8_solar_nue/{irun}.root")
        except:
            exists_b = False
        try:
            file_tl208 = ROOT.TFile.Open(f"../run_by_run_test/Tl208/{irun}.root")
        except:
            exists_tl = False
        
        # these arrays are initialised empty. If the file exists and has entries inside, 
        # they are replaced by the filled arrays and concatenated with the distributions
        # for the entire run list. However, if the file does not exist, add these empty
        # lists to the full distributions, which is the same as not adding at all.
        b8_e = []
        b8_l = [] 
        tl208_e = []
        tl208_l = []
        # if the file exists, get the ntuple and try and extract the entries' information
        # file can exist but be corrupted / empty, which leads to failures here
        if exists_b == True:
            # get all events in ROI 2.5 --> 5.0 MeV
            ntuple_b8    = file_b8.Get("2p5_5p0")
            try:
                b8_e, b8_l    = loop_over_entries(ntuple_b8, fv_cut)
            except:
                failedB8+=1
                
        if exists_tl ==  True:
            ntuple_tl208 = file_tl208.Get("2p5_5p0")
            try:
                tl208_e, tl208_l = loop_over_entries(ntuple_tl208, fv_cut)
            except:
                
                failedTl208 +=1

        # add the extracted information to the full distributions
        energy_b8     += b8_e
        multi_b8      += b8_l
        energy_tl208  += tl208_e
        multi_tl208   += tl208_l
    print(f"Failed runs B8: {failedB8}\nFailed runs Tl208: {failedTl208}.")
    energy_b8 = np.array(energy_b8)
    energy_tl208 = np.array(energy_tl208)
    multi_b8 = np.array(multi_b8)
    multi_tl208 = np.array(multi_tl208)
    
    # create the 4 sub arrays in each energy range
    pdf_b8_1 = multi_b8[(energy_b8 >= energy_pdf_bins[0]) & (energy_b8 < energy_pdf_bins[1])]
    pdf_b8_2 = multi_b8[(energy_b8 >= energy_pdf_bins[1]) & (energy_b8 < energy_pdf_bins[2])]
    pdf_b8_3 = multi_b8[(energy_b8 >= energy_pdf_bins[2]) & (energy_b8 < energy_pdf_bins[3])]
    pdf_b8_4 = multi_b8[(energy_b8 >= energy_pdf_bins[3]) & (energy_b8 < energy_pdf_bins[4])]

    pdf_tl208_1 = multi_tl208[(energy_tl208 >= energy_pdf_bins[0]) & (energy_tl208 < energy_pdf_bins[1])]
    pdf_tl208_2 = multi_tl208[(energy_tl208 >= energy_pdf_bins[1]) & (energy_tl208 < energy_pdf_bins[2])]
    pdf_tl208_3 = multi_tl208[(energy_tl208 >= energy_pdf_bins[2]) & (energy_tl208 < energy_pdf_bins[3])]
    pdf_tl208_4 = multi_tl208[(energy_tl208 >= energy_pdf_bins[3]) & (energy_tl208 < energy_pdf_bins[4])]

    # bin them and return the total counts
    # bin this information in E and M and create the normalised PDFs
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (17, 12))

    counts_b8_1, _, _ = axes[0,0].hist(pdf_b8_1, linewidth = 2, bins = multi_pdf_bins, color = "black", histtype = "step", label = r"$^8$B | $N_{events}$: " + f"{len(pdf_b8_1)}", density = False, weights = [1/len(pdf_b8_1)] * len(pdf_b8_1))
    counts_tl_1, _, _ = axes[0,0].hist(pdf_tl208_1, linewidth = 2, bins = multi_pdf_bins, color = "orange", histtype = "step", label = r"$^{208}$Tl | $N_{events}$: " + f"{len(pdf_tl208_1)}", weights = [1/len(pdf_tl208_1)] * len(pdf_tl208_1), density = False)
    axes[0,0].legend(frameon = False, fontsize = 15)
    axes[0,0].set_title(f"{energy_pdf_bins[0]}" + r" $\rightarrow$ " + f"{energy_pdf_bins[1]} MeV", fontsize = 20)
    axes[0,0].set_xlabel("Multisite Discriminant", fontsize = 20)
    axes[0,0].set_ylabel("Normalised Counts", fontsize = 20)

    counts_b8_2, _, _ = axes[0,1].hist(pdf_b8_2, linewidth = 2, bins = multi_pdf_bins, color = "black", histtype = "step", label = r"$^8$B | $N_{events}$: " + f"{len(pdf_b8_2)}", weights = [1/len(pdf_b8_2)] * len(pdf_b8_2), density = False)
    counts_tl_2, _, _ = axes[0,1].hist(pdf_tl208_2, linewidth = 2, bins = multi_pdf_bins, color = "orange", histtype = "step", label = r"$^{208}$Tl | $N_{events}$: " + f"{len(pdf_tl208_2)}", weights = [1/len(pdf_tl208_2)] * len(pdf_tl208_2), density = False)
    axes[0,1].legend(frameon = False, fontsize = 15)
    axes[0,1].set_title(f"{energy_pdf_bins[1]}" + r" $\rightarrow$ " + f"{energy_pdf_bins[2]} MeV", fontsize = 20)
    axes[0,1].set_xlabel("Multisite Discriminant", fontsize = 20)
    axes[0,1].set_ylabel("Normalised Counts", fontsize = 20)

    counts_b8_3, _, _ = axes[1,0].hist(pdf_b8_3, linewidth = 2, bins = multi_pdf_bins, color = "black", histtype = "step", label = r"$^8$B | $N_{events}$: " + f"{len(pdf_b8_3)}", weights = [1/len(pdf_b8_3)] * len(pdf_b8_3), density = False)
    counts_tl_3, _, _ = axes[1,0].hist(pdf_tl208_3, linewidth = 2, bins = multi_pdf_bins, color = "orange", histtype = "step", label = r"$^{208}$Tl | $N_{events}$: " + f"{len(pdf_tl208_3)}", weights = [1/len(pdf_tl208_3)] * len(pdf_tl208_3), density = False)
    axes[1,0].legend(frameon = False, fontsize = 15)
    axes[1,0].set_title(f"{energy_pdf_bins[2]}" + r" $\rightarrow$ " + f"{energy_pdf_bins[3]} MeV", fontsize = 20)
    axes[1,0].set_xlabel("Multisite Discriminant", fontsize = 20)
    axes[1,0].set_ylabel("Normalised Counts", fontsize = 20)

    counts_b8_4, _, _ = axes[1,1].hist(pdf_b8_4, linewidth = 2, bins = multi_pdf_bins, color = "black", histtype = "step", label = r"$^8$B | $N_{events}$: " + f"{len(pdf_b8_4)}", weights = [1/len(pdf_b8_4)] * len(pdf_b8_4), density = False)
    counts_tl_4, _, _ = axes[1,1].hist(pdf_tl208_4, linewidth = 2, bins = multi_pdf_bins, color = "orange", histtype = "step", label = r"$^{208}$Tl | $N_{events}$: " + f"{len(pdf_tl208_4)}", weights = [1/len(pdf_tl208_4)] * len(pdf_tl208_4), density = False)
    axes[1,1].legend(frameon = False, fontsize = 15)
    axes[1,1].set_title(f"{energy_pdf_bins[3]}" + r" $\rightarrow$ " + f"{energy_pdf_bins[4]} MeV", fontsize = 20)
    axes[1,1].set_xlabel("Multisite Discriminant", fontsize = 20)
    axes[1,1].set_ylabel("Normalised Counts", fontsize = 20)
    fig.tight_layout()
    plt.savefig("../plots/asimov_study/real_mc/multisite_pdfs.png")
    plt.close()

    print(f"Sum of bin counts: {np.round(np.sum(counts_b8_1), 3)}, {np.round(np.sum(counts_b8_2), 3)}, {np.round(np.sum(counts_b8_3), 3)}, {np.round(np.sum(counts_b8_4), 3)}, {np.round(np.sum(counts_tl_1), 3)}, {np.round(np.sum(counts_tl_2), 3)}, {np.round(np.sum(counts_tl_3), 3)}, {np.round(np.sum(counts_tl_4), 3)}")
    
    # take these normalised PDFs and add them to arrays for easy movement
    normalised_pdfs_b8 = np.zeros((4, len(multi_pdf_bins)-1)) # dims (N_energy_bins, N_multisite_bins)
    normalised_pdfs_tl = np.zeros((4, len(multi_pdf_bins)-1))
    normalised_pdfs_b8[0, :] = counts_b8_1
    normalised_pdfs_b8[1, :] = counts_b8_2
    normalised_pdfs_b8[2, :] = counts_b8_3
    normalised_pdfs_b8[3, :] = counts_b8_4
    normalised_pdfs_tl[0, :] = counts_tl_1
    normalised_pdfs_tl[1, :] = counts_tl_2
    normalised_pdfs_tl[2, :] = counts_tl_3
    normalised_pdfs_tl[3, :] = counts_tl_4

    # pad the zero values
    normalised_pdfs_b8 = normalised_pdfs_b8.astype(np.float32)
    normalised_pdfs_tl = normalised_pdfs_tl.astype(np.float32)
    normalised_pdfs_b8[normalised_pdfs_b8 == 0] = 1e-6
    normalised_pdfs_tl[normalised_pdfs_tl == 0] = 1e-6
    
    # create an asimov dataset in each energy bin
    asimov_datasets = np.zeros((4, len(multi_pdf_bins) - 1))

    # asimov dataset is the signal and background pdfs in each energy bin scaled by the desired number of counts
    asimov_datasets[0, :] = num_sig[0] * counts_b8_1 + num_back[0] * counts_tl_1
    asimov_datasets[1, :] = num_sig[1] * counts_b8_2 + num_back[1] * counts_tl_2
    asimov_datasets[2, :] = num_sig[2] * counts_b8_3 + num_back[2] * counts_tl_3
    asimov_datasets[3, :] = num_sig[3] * counts_b8_4 + num_back[3] * counts_tl_4

    return normalised_pdfs_b8, normalised_pdfs_tl, asimov_datasets

def evaluate_likelihood(signal_pdfs, backg_pdfs, asimov_datasets, signal_hypothesis, true_back_num):
    """
    Evaluate the Poisson only, multisite only, poisson + multisite and chi2 for the asimov
    dataset and a range of signal hypothesis.

    Produce curves for each energy bin
    """

    # create arrays to save the poisson, chi2, multisite -2log(L) for each energy bin
    
    # find the longest signal hypothesis list (ie the one with the most events in it)
    idx_longest = np.argmax(true_back_num) # i.e. the one with the most background events in it
    chi2_data    = np.zeros((4, len(signal_hypothesis[idx_longest]))) # each row is the curve for an energy bin
    poisson_data = np.zeros((4, len(signal_hypothesis[idx_longest])))
    multi_data   = np.zeros((4, len(signal_hypothesis[idx_longest])))

    for ienergy_bin in range(4):

        # load up the asimov dataset and the corresponding signal and background pdfs
        signal_pdf = signal_pdfs[ienergy_bin, :]
        backg_pdf  = backg_pdfs[ienergy_bin, :]
        asimov     = asimov_datasets[ienergy_bin, :]
        tot_counts_asimov = np.sum(asimov)
        
        # loop over the signal hypothesis
        hypothesis = signal_hypothesis[ienergy_bin]
        for isig in range(len(hypothesis)):
            if tot_counts_asimov - hypothesis[isig] < 0:
                print("OOps! ", hypothesis[isig], f" for energy bin {ienergy_bin} with total counts {tot_counts_asimov}")
            # calculate the chi2
            chi2    = (tot_counts_asimov - (hypothesis[isig] + true_back_num[ienergy_bin]))**2 / (hypothesis[isig] + true_back_num[ienergy_bin])
            
            # calculate poisson likelihood
            poisson = 2 * ((hypothesis[isig] + true_back_num[ienergy_bin]) - tot_counts_asimov * np.log(hypothesis[isig]+true_back_num[ienergy_bin]))

            # calculate the multisite likelihood
            multi   = 2 * ((tot_counts_asimov) - np.sum(asimov * np.log(hypothesis[isig] * signal_pdf + (tot_counts_asimov - hypothesis[isig]) * backg_pdf)))

            chi2_data[ienergy_bin, isig]    = chi2
            poisson_data[ienergy_bin, isig] = poisson
            multi_data[ienergy_bin, isig]   = multi
    
    # calculate the full likelihood function of poisson + multisite parts
    full_data = poisson_data + multi_data

    # mask out all of the zeros that aren't updated in each curve
    chi2_mask = (chi2_data == 0)
    chi2_masked = np.ma.masked_array(chi2_data, chi2_mask)
    poisson_mask = (poisson_data == 0)
    poisson_masked = np.ma.masked_array(poisson_data, poisson_mask)
    multi_mask = (multi_data == 0)
    multi_masked = np.ma.masked_array(multi_data, multi_mask)
    full_mask = (full_data == 0)
    full_masked = np.ma.masked_array(full_data, full_mask)
    # shift all of the curves so the minima are at zero
    min_idx_chi2     = np.argmin(chi2_masked, axis = 1)
    min_idx_poisson  = np.argmin(poisson_masked, axis = 1)
    min_idx_multi    = np.argmin(multi_masked, axis = 1)
    min_idx_full     = np.argmin(full_masked, axis = 1)
    
    for i in range(4):
        zero_diff_chi2    = 0 - chi2_masked[i, min_idx_chi2[i]]
        zero_diff_poisson = 0 - poisson_masked[i, min_idx_poisson[i]]
        zero_diff_multi   = 0 - multi_masked[i, min_idx_multi[i]]
        zero_diff_full    = 0 - full_masked[i, min_idx_full[i]]

        # shift the curves by this amount
        chi2_masked[i, :]    = chi2_masked[i, :] + zero_diff_chi2
        poisson_masked[i, :] = poisson_masked[i, :] + zero_diff_poisson
        multi_masked[i, :]   = multi_masked[i, :] + zero_diff_multi
        full_masked[i, :]    = full_masked[i, :] + zero_diff_full
    return full_masked, multi_masked, poisson_masked, chi2_masked

fv_cut            = 4500                                        # FV cut applied to MC in mm
pdf_bins_energy   = np.arange(2.5, 5.0, 0.5)                    # energy bin edges in MeV
pdf_bins_multi    = np.arange(-1.13, -1.08, 0.0005)             # multisite discriminant bin edges
num_sig_events    = 87.8 * np.array([0.28, 0.26, 0.24, 0.22])   # number of B8 events in Asimov dataset
num_back_events   = 830.5 * np.array([0.048, 0.35, 0.44, 0.16]) # number of Tl208 events in Asimov dataset
step              = 0.1                                         # step size in signal hypothesis grid search

# create a signal hypothesis grid search array in each dimension, ensuring the scanned number of signal events does not exceed the
# total number of counts observed in that bin (otherwise get -ve number in multisite likelihood log!)
signal_hypothesis = []
for i in range(4):
    hypothesis_scan = np.arange(0, num_sig_events[i]+num_back_events[i], step)
    signal_hypothesis.append(hypothesis_scan)

# 1. create the PDFs, and use them to create an Asimov dataset
sig_pdfs, back_pdfs, asimov_datasets = create_pdfs_and_asimov(fv_cut, pdf_bins_energy, pdf_bins_multi, num_sig_events, num_back_events)

# 2. evaluate the likelihood in each bin
full, multi, poisson, chi2 = evaluate_likelihood(sig_pdfs, back_pdfs, asimov_datasets, signal_hypothesis, num_back_events)
print(full.shape)
# 3. create the output plots per energy bin

"""
Want to find the minimum of each curve and the 1sigma confidence interval. No way to do this
analytically (it's a nightmare!), so I will do this 'numerically'.
"""

# create an array (N_energy_bins, 2) where each column is 1sigma interval start, 1sigma interval end
confidence_intervals = np.zeros((4, 6))

# loop over energy bins and find the start and end of each confidence interval
for i in range(4):
    # find an array of absolute differences to the sigma level and then return idx of array element corresponding to it
    idx_1sig_full    = np.argmin(np.abs(full[i, :][full[i, :] != 0] - 1.0))  # 1.0 is the 1 sigma level of the -2log(L) plots
    idx_1sig_poisson = np.argmin(np.abs(poisson[i, :][poisson[i, :] != 0] - 1.0))
    idx_1sig_multi   = np.argmin(np.abs(multi[i, :][multi[i, :] != 0] - 1.0))

    # find the second smallest difference (i.e. the other end of the confidence interval)
    masked_1sig_full = np.ma.masked_array(full[i, :][full[i, :] != 0], mask = False).copy()

    # now remove the smallest value
    masked_1sig_full.mask[idx_1sig_full] = True # True means value is masked out and not considered in subsequent operations
    
    # repeat finding the smallest difference
    idx2_1sig_full = np.argmin(np.abs(masked_1sig_full-1.0))

    # repeat for the poisson likelihood
    masked_1sig_poisson = np.ma.masked_array(poisson[i, :][poisson[i, :] != 0], mask = False).copy()
    masked_1sig_poisson.mask[idx_1sig_poisson] = True
    idx2_1sig_poisson = np.argmin(np.abs(masked_1sig_poisson-1.0))

    # repeat for multisite likelihood
    masked_1sig_multi = np.ma.masked_array(multi[i, :][multi[i, :] != 0], mask = False).copy()
    masked_1sig_multi.mask[idx_1sig_multi] = True
    idx2_1sig_multi = np.argmin(np.abs(masked_1sig_multi-1.0))
    print(idx_1sig_full, idx2_1sig_full)
    confidence_intervals[i, 0] = signal_hypothesis[i][idx_1sig_full]
    confidence_intervals[i, 1] = signal_hypothesis[i][idx2_1sig_full]
    confidence_intervals[i, 2] = signal_hypothesis[i][idx_1sig_poisson]
    confidence_intervals[i, 3] = signal_hypothesis[i][idx2_1sig_poisson]
    confidence_intervals[i, 4] = signal_hypothesis[i][idx_1sig_multi]
    confidence_intervals[i, 5] = signal_hypothesis[i][idx2_1sig_multi]

print(confidence_intervals)

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (17, 12))
f_size_labels      = 20
f_size_tick_labels = 15
f_size_legend      = 15
x_axis_range_perc  = 0.95
axes[0,0].plot(signal_hypothesis[0], chi2[0, :].compressed(), color = "blue", label = r"$\chi ^2$")
axes[0,0].plot(signal_hypothesis[0], poisson[0, :].compressed(), color = "orange", label = f"Poisson likelihood | {np.round(float(signal_hypothesis[0][np.argmin(poisson[0, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[0, 2]-confidence_intervals[0, 3])/2, 2)}")
axes[0,0].vlines(x = confidence_intervals[0, 2], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[0,0].vlines(x = confidence_intervals[0, 3], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[0,0].plot(signal_hypothesis[0], multi[0, :].compressed(), color = "green", label = f"Multisite likelihood | {np.round(float(signal_hypothesis[0][np.argmin(multi[0, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[0, 4]-confidence_intervals[0, 5])/2, 2)}")
axes[0,0].vlines(x = confidence_intervals[0, 4], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[0,0].vlines(x = confidence_intervals[0, 5], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[0,0].plot(signal_hypothesis[0], full[0, :].compressed(), color = "black", label = f"Full likelihood | {np.round(float(signal_hypothesis[0][np.argmin(full[0, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[0, 0]-confidence_intervals[0, 1])/2, 2)}")
axes[0,0].vlines(x = confidence_intervals[0, 0], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[0,0].vlines(x = confidence_intervals[0, 1], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[0,0].axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$", linewidth = 2)
axes[0,0].set_xlabel(r"$^8$B Signal Events", fontsize = f_size_labels)
axes[0,0].set_ylabel(r"$-2log(\mathcal{L})$", fontsize = f_size_labels)
axes[0,0].legend(frameon = False, fontsize = f_size_legend)
axes[0,0].set_title("2.5 " + r"$\rightarrow$" + " 3.0 MeV", fontsize = f_size_labels)
axes[0,0].set_ylim((0, 3))
axes[0,0].set_xlim((num_sig_events[0]-num_sig_events[0]*x_axis_range_perc, num_sig_events[0]+ num_sig_events[0]*x_axis_range_perc))
axes[0,0].tick_params(axis = "both", labelsize = f_size_tick_labels)

axes[0,1].plot(signal_hypothesis[1], chi2[1, :].compressed(), color = "blue", label = r"$\chi^2$")
axes[0,1].plot(signal_hypothesis[1], poisson[1, :].compressed(), color = "orange", label = f"Poisson likelihood | {np.round(float(signal_hypothesis[1][np.argmin(poisson[1, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[1, 2]-confidence_intervals[1, 3])/2, 2)}")
axes[0,1].plot(signal_hypothesis[1], multi[1, :].compressed(), color = "green", label = f"Multisite likelihood | {np.round(float(signal_hypothesis[1][np.argmin(multi[1, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[1, 4]-confidence_intervals[1, 5])/2, 2)}")
axes[0,1].plot(signal_hypothesis[1], full[1, :].compressed(), color = "black", label = f"Full likelihood | {np.round(float(signal_hypothesis[1][np.argmin(full[1, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[1, 0]-confidence_intervals[1, 1])/2, 2)}")
axes[0,1].vlines(x = confidence_intervals[1, 2], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[0,1].vlines(x = confidence_intervals[1, 3], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[0,1].vlines(x = confidence_intervals[1, 4], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[0,1].vlines(x = confidence_intervals[1, 5], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[0,1].vlines(x = confidence_intervals[1, 0], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[0,1].vlines(x = confidence_intervals[1, 1], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[0,1].axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$", linewidth = 2)
axes[0,1].legend(frameon = False, fontsize = f_size_legend)
axes[0,1].set_title("3.0 " + r"$\rightarrow$" + " 3.5 MeV", fontsize = f_size_labels)
axes[0,1].set_xlabel(r"$^8$B Signal Events", fontsize = f_size_labels)
axes[0,1].set_ylabel(r"$-2log(\mathcal{L})$", fontsize = f_size_labels)
axes[0,1].set_ylim((0, 3))
axes[0,1].set_xlim((num_sig_events[1]-num_sig_events[1]*x_axis_range_perc, num_sig_events[1]+ num_sig_events[1]*x_axis_range_perc))
axes[0,1].tick_params(axis = "both", labelsize = f_size_tick_labels)

axes[1,0].plot(signal_hypothesis[2], chi2[2, :].compressed(), color = "blue", label = r"$\chi^2$")
axes[1,0].plot(signal_hypothesis[2], poisson[2, :].compressed(), color = "orange", label = f"Poisson likelihood | {np.round(float(signal_hypothesis[2][np.argmin(poisson[2, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[2, 2]-confidence_intervals[2, 3])/2, 2)}")
axes[1,0].plot(signal_hypothesis[2], multi[2, :].compressed(), color = "green", label = f"Multisite likelihood | {np.round(float(signal_hypothesis[2][np.argmin(multi[2, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[2, 4]-confidence_intervals[2, 5])/2, 2)}")
axes[1,0].plot(signal_hypothesis[2], full[2, :].compressed(), color = "black", label = f"Full likelihood | {np.round(float(signal_hypothesis[2][np.argmin(full[2, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[2, 0]-confidence_intervals[2, 1])/2, 2)}")
axes[1,0].vlines(x = confidence_intervals[2, 2], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[1,0].vlines(x = confidence_intervals[2, 3], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[1,0].vlines(x = confidence_intervals[2, 4], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[1,0].vlines(x = confidence_intervals[2, 5], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[1,0].vlines(x = confidence_intervals[2, 0], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[1,0].vlines(x = confidence_intervals[2, 1], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[1,0].axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$", linewidth = 2)
axes[1,0].legend(frameon = False, fontsize = f_size_legend)
axes[1,0].set_title("3.5 " + r"$\rightarrow$" + " 4.0 MeV", fontsize = f_size_labels)
axes[1,0].set_xlabel(r"$^8$B Signal Events", fontsize = f_size_labels)
axes[1,0].set_ylabel(r"$-2log(\mathcal{L})$", fontsize = f_size_labels)
axes[1,0].set_ylim((0, 3))
axes[1,0].set_xlim((num_sig_events[2]-num_sig_events[2]*x_axis_range_perc, num_sig_events[2]+ num_sig_events[2]*x_axis_range_perc))
axes[1,0].tick_params(axis = "both", labelsize = f_size_tick_labels)

axes[1,1].plot(signal_hypothesis[3], chi2[3, :].compressed(), color = "blue", label = r"$\chi^2$")
axes[1,1].plot(signal_hypothesis[3], poisson[3, :].compressed(), color = "orange", label = f"Poisson likelihood | {np.round(float(signal_hypothesis[3][np.argmin(poisson[3, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[3, 2]-confidence_intervals[3, 3])/2, 2)}")
axes[1,1].plot(signal_hypothesis[3], multi[3, :].compressed(), color = "green", label = f"Multisite likelihood | {np.round(float(signal_hypothesis[3][np.argmin(multi[3, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[3, 4]-confidence_intervals[3, 5])/2, 2)}")
axes[1,1].plot(signal_hypothesis[3], full[3, :].compressed(), color = "black", label = f"Full likelihood | {np.round(float(signal_hypothesis[3][np.argmin(full[3, :].compressed())]), 2)}" + r" $\pm$ " + f"{np.round(np.abs(confidence_intervals[3, 0]-confidence_intervals[3, 1])/2, 2)}")
axes[1,1].vlines(x = confidence_intervals[3, 2], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[1,1].vlines(x = confidence_intervals[3, 3], ymin = 0, ymax = 1, color = "orange", linestyle = "dashed")
axes[1,1].vlines(x = confidence_intervals[3, 4], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[1,1].vlines(x = confidence_intervals[3, 5], ymin = 0, ymax = 1, color = "green", linestyle = "dashed")
axes[1,1].vlines(x = confidence_intervals[3, 0], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[1,1].vlines(x = confidence_intervals[3, 1], ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axes[1,1].axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$", linewidth = 2)
axes[1,1].legend(frameon = False, fontsize = f_size_legend)
axes[1,1].set_title("4.0 " + r"$\rightarrow$" + " 4.5 MeV", fontsize = f_size_labels)
axes[1,1].set_xlabel(r"$^8$B Signal Events", fontsize = f_size_labels)
axes[1,1].set_ylabel(r"$-2log(\mathcal{L})$", fontsize = f_size_labels)
axes[1,1].set_ylim((0, 3))
axes[1,1].set_xlim((num_sig_events[3]-num_sig_events[3]*x_axis_range_perc, num_sig_events[3]+ num_sig_events[3]*x_axis_range_perc))
axes[1,1].tick_params(axis = "both", labelsize = f_size_tick_labels)
fig.tight_layout()
plt.savefig("../plots/asimov_study/real_mc/result_test.png")
plt.close()