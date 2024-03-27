import numpy as np
import matplotlib
matplotlib.use("Agg")
import ROOT
import matplotlib.pyplot as plt

class AsimovLikelihoodExtraction():
    """
    1. Create Asimov dataset from Tl208 and B8 MC, with a given fraction of B8 in the 2.5 --> 5 MeV region of interest.
    2. Chi2 minimisation of B8 number in each energy bin
    3. Multisite Discrimination to return optimal number B8 in each energy bin
    4. Comparison plots of each method discriminant vs B8 number vs true number of B8 in each bin.
    """
    
    def apply_scaling(self, energy_spectrum, expected_number, binning):
        """
        Function takes in the unbinned energy spectrum from MC, bins and scales it such that the sum of
        the counts is equal to the expected number of events given by the background model calculations.

        The poisson error on these scaled counts is calculated.

        INPUTS : raw energy spectrum of MC background
                expected number of events given by background model
                energy spectrum binning to be applied

        RETURNS: scaled counts, poisson error on counts, midpoint of bins.
        """

        # bin the raw energy spectrum
        unscaled_counts, bins = np.histogram(energy_spectrum, bins = binning)

        # scale the counts such that the sum of counts equals expected number of events
        scaled_counts = ( unscaled_counts / sum(unscaled_counts) ) * expected_number

        # find the midpoint of the bins for plotting
        mid_points = binning + np.diff(binning)[0] / 2

        # calculate the poisson error on the counts of each bin # 
        ################## TO BE IMPLEMENTED ####################

        return scaled_counts, mid_points[:-1] 

    def loop_over_entries(self, ntuple, fv_cut):
        
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
    
    def create_asimov_and_pdfs(self, fraction_b8, fv_cut, binning, pdf_binning):
        """
        Function loads up Tl208 and B8 MC and creates a combined data spectrum. This spectrum is then binned in
        the ROI (2.5 --> 5 MeV) range, and scaled so the B8 is a given fraction of the spectrum.

        Also builds the dlogL PDFs for use in the multisite likelihood.

        fraction_b8, float, between 0 --> 1, % of asimov spectrum B8 events
        fv_Cut, float, FV cut to apply in mm
        binning, np.array, defined bin edges to bin data spectrum in
        """

        # load up the mc run list
        runlist = np.loadtxt("../runlists/additional_test_runlist.txt", dtype = int)

        # use half the MC for the asimov dataset - the other half will be used to create the dlog(L) PDFs
        split_model = int(len(runlist) *0.25)
        print(split_model)
        asimov_runs = runlist[:split_model]  # runs to build dataset
        pdf_runs   = runlist[split_model:]  # different runs to build dlog(L) PDFs

        energy_b8        = [] # for asimov data spectrum
        energy_tl208     = [] # for asimov data spectrum
        energy_dataset   = [] # combined energy spectrum of dataset
        dlogL_dataset    = [] # combined dlogL spectrum of dataset
        dlogL_b8_pdf     = [] # for dlogL PDF creation
        dlogL_tl208_pdf  = [] # for dlogL PDF creation
        energy_b8_pdf    = [] # for dlogL PDF creation
        energy_tl208_pdf = [] # for dlogL PDF creation
        failedB8 = 0
        failedTl208 = 0
        # loop over the asimov runs
        for irun in asimov_runs:
            file_b8    = ROOT.TFile.Open(f"../run_by_run_test/B8_solar_nue/{irun}.root")
            file_tl208 = ROOT.TFile.Open(f"../run_by_run_test/Tl208/{irun}.root")

            # extract information from the 2.0 --> 5.0 MeV branch
            ntuple_b8    = file_b8.Get("2p5_5p0")
            ntuple_tl208 = file_tl208.Get("2p5_5p0")

            try:
                b8_e, b8_l    = self.loop_over_entries(ntuple_b8, fv_cut)
            except:
                failedB8+=1 
                pass
            try:
                tl208_e, tl208_l = self.loop_over_entries(ntuple_tl208, fv_cut)
            except:
                failedTl208 +=1 
                pass
            energy_b8     += b8_e
            energy_tl208  += tl208_e
            dlogL_dataset += b8_l
            dlogL_dataset += tl208_l
            energy_dataset += b8_e
            energy_dataset += tl208_e
        print(f"Failed B8: {failedB8}\nFailed Tl208: {failedTl208}")
        energy_dataset = np.array(energy_dataset)
        dlogL_dataset  = np.array(dlogL_dataset)

        # create dlogL distributions for each energy bin
        DLOGL_DATA_PER_BIN = [] # list of distributions
        for i in range(len(binning) -1):
            min_energy = binning[i]
            max_energy = binning[i+1]

            dlog = dlogL_dataset[(energy_dataset >= min_energy) & (energy_dataset < max_energy)]
            DLOGL_DATA_PER_BIN.append(dlog)

        # now we have the energy spectrum between 2.5 --> 5 MeV we can bin it
        asimov_b8, _    = np.histogram(energy_b8, bins = binning)
        asimov_tl208, _ = np.histogram(energy_tl208, bins = binning)

        error_b8    = np.sqrt(asimov_b8) 
        error_tl208 = np.sqrt(asimov_tl208)
        error_total = np.sqrt(asimov_b8 + asimov_tl208)
        total_spectrum = np.sum(asimov_b8 + asimov_tl208)
        total_b8       = np.sum(asimov_b8)
        total_tl208    = np.sum(asimov_tl208)
        print(f"Before scaling Numbers:\nTl208: {total_tl208} | {round(total_tl208/total_spectrum, 2)}\nB8: {total_b8} | {round(total_b8/total_spectrum, 2)}\nTotal: {total_spectrum}")

        # find the desired number of B8 events in the spectrum, given the input % B8
        desired_b8_number    = total_spectrum * fraction_b8
        desired_tl208_number = total_spectrum * (1-fraction_b8)

        # apply this scaling to the asimov energy spectrum of each component
        scaled_b8, _    = self.apply_scaling(energy_b8, desired_b8_number, binning)
        scaled_tl208, _ = self.apply_scaling(energy_tl208, desired_tl208_number, binning) 
        
        total_asimov_counts  = np.sum(scaled_b8 + scaled_tl208)
        print(f"Fraction of events B8 after scaling: {np.sum(scaled_b8)/total_asimov_counts}")

        # add the counts together to get the total asimov dataset
        asimov_data  = scaled_b8 + scaled_tl208
        asimov_data  = asimov_b8 + asimov_tl208
        scaled_b8 = asimov_b8
        scaled_tl208 = asimov_tl208
        # plot the asimov dataset
        plt.figure(figsize = (12, 6))
        plt.errorbar(binning[0:-1] + np.diff(binning)[0]/2, asimov_data, yerr = np.sqrt(asimov_data), xerr = np.diff(binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "black", label = "Total Spectrum")
        plt.errorbar(binning[0:-1] + np.diff(binning)[0]/2, scaled_b8, yerr = np.sqrt(scaled_b8), xerr = np.diff(binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "red", label = "B8")
        plt.errorbar(binning[0:-1] + np.diff(binning)[0]/2, scaled_tl208, yerr = np.sqrt(scaled_tl208), xerr = np.diff(binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "green", label = "Tl208")
        plt.xlim((2.5, 4.5))
        plt.ylim((0, np.max(asimov_data)* 1.1))
        plt.legend()
        plt.xlabel("Reconstructed Energy (MeV)", fontsize = 20)
        plt.ylabel(f"Counts per {np.diff(binning)[0]} MeV Bin", fontsize = 20)
        plt.title("Asimov Dataset", fontsize = 20)
        plt.savefig("../plots/asimov_study/asimov_dataset.png")
        plt.close()

        #############################
        # create the multisite PDFs #
        #############################

        # create some output plots of the PDFS for each energy bin
        fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (17, 12))

        # loop over all the runs in the PDF runlist
        # fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 6))
        # counter = 0
        for irun in pdf_runs:
            # if counter > 9: 
            #     break
            file_b8    = ROOT.TFile.Open(f"../run_by_run_test/B8_solar_nue/{irun}.root")
            file_tl208 = ROOT.TFile.Open(f"../run_by_run_test/Tl208/{irun}.root")

            # extract information from the 2.0 --> 5.0 MeV branch
            ntuple_b8    = file_b8.Get("2p5_5p0")
            ntuple_tl208 = file_tl208.Get("2p5_5p0")

            try:
                b8_e, b8_l       = self.loop_over_entries(ntuple_b8, fv_cut)
            except:
                pass
            try:
                tl208_e, tl208_l = self.loop_over_entries(ntuple_tl208, fv_cut)
            except:
                pass
            
            # save the energy and dlogL of each event
            dlogL_b8_pdf     += b8_l
            dlogL_tl208_pdf  += tl208_l

            
            # axes[0].hist(b8_l, bins = 10, density = True, histtype = "step", label = f"Run {irun}")
            # axes[1].hist(tl208_l, bins = 10, density=True, histtype = "step", label = f"Run {irun}")
            # axes[0].legend()
            # axes[1].legend()
            # counter += 1
            energy_b8_pdf    += b8_e
            energy_tl208_pdf += tl208_e
        # plt.tight_layout()
        # plt.savefig("../plots/asimov_study/understanding/run_id_dlogL.pdf")
        # plt.close()
        # cast to a np array for numpy conditional indexing
        dlogL_b8_pdf     = np.array(dlogL_b8_pdf)
        dlogL_tl208_pdf  = np.array(dlogL_tl208_pdf)
        energy_tl208_pdf = np.array(energy_tl208_pdf)
        energy_b8_pdf    = np.array(energy_b8_pdf)
        
        # for each bin in the energy spectrum, create a dlog(L) PDF
        LIST_OF_DLOGL_B8_PDFS    = []
        LIST_OF_DLOGL_TL208_PDFS = []
        row = 0
        col = 0
        count = 0
        for i in range(len(binning) - 1):
            if count == 2:
                row = 1
                col = 0
            
            # define the energy bins for PDF
            min_bin_energy = binning[i]
            max_bin_energy = binning[i+1]

            # apply energy cuts to extract dlogL values corresponding to this energy bin
            pdf_vals_b8    = dlogL_b8_pdf[(energy_b8_pdf >= min_bin_energy) & (energy_b8_pdf < max_bin_energy)]
            pdf_vals_tl208 = dlogL_tl208_pdf[(energy_tl208_pdf >= min_bin_energy) & (energy_tl208_pdf < max_bin_energy)] 

            # update the PDF plot
            axes[row, col].hist(pdf_vals_b8, bins = pdf_binning, density = True,  color = "black", linewidth = 2, histtype = "step", label = f"B8 | Events: {len(pdf_vals_b8)}")
            axes[row, col].hist(pdf_vals_tl208, bins = pdf_binning, density = True, color = "orange", linewidth = 2, histtype = "step", label = f"Tl208 | Events: {len(pdf_vals_tl208)}")
            axes[row, col].set_title(f"{min_bin_energy}" + r"$\rightarrow$" +  f"{max_bin_energy} MeV", fontsize = 20)
            axes[row, col].set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = 20)
            axes[row, col].set_ylabel("Counts", fontsize = 20)
            axes[row, col].legend(fontsize = 15)
            axes[row, col].tick_params(axis = "x", rotation = 45, labelsize = 15)
            axes[row, col].tick_params(axis = "y", labelsize = 15)
            # axes.hist(pdf_vals_b8, bins = pdf_binning, density = True,  histtype = "step", label = f"B8 | Events: {len(pdf_vals_b8)}")
            # axes.hist(pdf_vals_tl208, bins = pdf_binning, density = True, histtype = "step", label = f"Tl208 | Events: {len(pdf_vals_tl208)}")
            # axes.set_title(f"{min_bin_energy}" + r"$\rightarrow$" +  f"{max_bin_energy} MeV", fontsize = 20)
            # axes.set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = 20)
            # axes.set_ylabel("Counts", fontsize = 20)
            # axes.legend()
            # axes.tick_params(axis = "x", rotation = 45)

            # save the PDF array into a list. Each entry in the list is a PDF
            LIST_OF_DLOGL_B8_PDFS.append(pdf_vals_b8)
            LIST_OF_DLOGL_TL208_PDFS.append(pdf_vals_tl208)

            count +=1
            col += 1
        
        fig.tight_layout()
        plt.savefig("../plots/asimov_study/dlogL_pdfs.png")
        plt.close()
        
        # return the total data spectrum and the number of Tl208 and B8 counts in each bin
        # along with the arrays of corresponding dlogL PDFs
        return asimov_data, scaled_b8, scaled_tl208, DLOGL_DATA_PER_BIN, LIST_OF_DLOGL_B8_PDFS, LIST_OF_DLOGL_TL208_PDFS

    def perform_loglikelihood_minimisation(self, number_data, number_model, values_b8):
        """
        For a given energy bin, a simple grid search over number of B8 events
        is performed. For each number of B8 tested, the chi2 (-2log(L)) is 
        evaluated and saved. 
        
        The chi2 vs number of B8 events in the model is returned.
        """

        chi2 = []
        for i in range(len(values_b8)):

            # increase the number in the model by the number of B8 events assumed
            updated_model = number_model + values_b8[i]

            # calculate the chi2 between model and data for this number of B8
            fit_goodness = (number_data - updated_model)**2 / updated_model

            chi2.append(fit_goodness)

        return chi2
    
    def perform_multisite_minimisation(self, values_b8, value_tl208, dlogL_data, pdf_b8, pdf_tl208, pdf_binning, scaled_counts):
        """
        Function performs the B8 fit taking into account the multisite PDFs.
        
        INPUTS: 
                values_b8, list of assumed numbers of B8 events in given energy bin
                value_tl208, float, fixed number of Tl208 events assumed by model in energy bin
                dlogL_data, array, value of multisite discriminant for every event in this energy bin
                pdf_b8, array, normalised PDF for given energy bin for single-site events (counts)
                pdf_tl208, array, normalised PDF for given energy bin for multi-site events (counts)
                pdf_binning, the bin edges used to bin the dlogL pdfs
                chi2, array, the chi2 prefactor found previously for each assumed number of B8 events

        RETURNS: 
                loglikelihood, array, array of loglikelihood vs assumed B8 numbers
        """

        # calculate the probability each event is signal given dlogL
        multisite_bin_idx  = np.digitize(dlogL_data, bins = pdf_binning) - 1 # the bin IDX each event's dlogL value falls into
        
        plt.figure()
        plt.plot(pdf_binning[:-1] + np.diff(pdf_binning)[0]/2, pdf_tl208, label = "Tl208 PDF")
        plt.plot(pdf_binning[:-1] + np.diff(pdf_binning)[0]/2, pdf_b8, label = "B8 PDF")
        plt.hist(dlogL_data, bins = 50, histtype = "step", label = "dlogL Data", density = True, linewidth = 2)
        plt.xlabel("Multisite Discriminant Value")
        plt.legend()
        plt.savefig(f"../plots/asimov_study/understanding/dlogL_data_{i_energy_bin}.png")
        plt.close()
        # find the signal and background probability arrays from binned data dlogL and pdfs
        # these arrays do not change so no need to recompute the loop over events
        probability_signal     = pdf_b8[multisite_bin_idx]
        probability_background = pdf_tl208[multisite_bin_idx]
        print(probability_signal)
        print(value_tl208)
        print(f"Num data events in bin: {len(dlogL_data)}")
        # loop over each assumed number of B8 events
        loglikelihood  = []
        multisite_vals = []
        poisson_vals   = []
        output = []
        output_sig = []
        output_back = []
        for i in range(len(values_b8)):
            # if i >= 1000:
            #     break
            model_expectation = values_b8[i] + value_tl208
            
            # multisite part of likelihood
            fraction_b8 = values_b8[i] / model_expectation
            # print(model_expectation, fraction_b8)
            multisite   = fraction_b8 * probability_signal + (1-fraction_b8) * probability_background
            # multisite = fraction_b8 * probability_background + (1-fraction_b8) * probability_signal
            # multisite[multisite == 0] = 1e9
            # print(f"Prob Signal: {probability_signal}\nProb Backg: {probability_background}\nFrac B8: {fraction_b8}\n")
            # output.append(np.sum(np.log(multisite)))
            # output_sig.append(fraction_b8* np.sum(probability_signal))
            # output_back.append((1-fraction_b8)* np.sum(probability_background))
            multisite   = np.sum(-np.log(multisite))
            # print(multisite)
            # multisite = np.prod(multisite)
            #multisite   = np.sum(multisite)

            # poisson statistics part of likelihood
            poisson = model_expectation - scaled_counts * np.log(model_expectation)
            # print(f"Poisson Part: {poisson}\nMultisite Part: {multisite}")
            # loglikelihood.append(poisson-multisite)
            multisite_vals.append(2*multisite)
            poisson_vals.append(2*poisson)
            # loglikelihood.append(multisite)
            
            loglikelihood.append(2*poisson + 2*multisite)
            # print(multisite)

            # multiply by chi2 prefactor
            # val = chi2[i] * val

            # take log
            # loglikelihood.append(-2*np.log(val))

        return loglikelihood, poisson_vals

# define inputs to the analysis function
energy_binning = np.arange(2.5, 5.0, 0.5)
pdf_binning    = np.arange(-1.13, -1.08, 0.0005)
fraction_b8    = 0.1
fv_cut         = 4500
values_b8      = np.arange(0, 7000, 1)
# values_b8 = values_b8[:1000]
print(values_b8[-1])
print("Number B8 values: ", len(values_b8))
# create THE ANALYSIS OBJECT
X = AsimovLikelihoodExtraction()

# create asimov dataset, the true number of events in each energy bin, and the bin-by-bin dlogL PDFs for multisite
asimov_dataset, num_b8_per_bin, num_tl208_per_bin, dlogL_data, pdfs_b8, pdfs_tl208 = X.create_asimov_and_pdfs(fraction_b8, fv_cut, energy_binning, pdf_binning)
print(num_b8_per_bin, num_tl208_per_bin)
# print(asimov_dataset[1], len(dlogL_data[1]), num_tl208_per_bin[1])
# create an output plot showing the chi2 minimisation
fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (17, 12))
# fig, axes = plt.subplots(nrows = 1, ncols = len(asimov_dataset), figsize = (12, 6))
# run bin-by-bin chi2 fitting for number of B8 events in each bin
chi2_minimisation = [] # list of lists --> each contains the chi2 for a given number of B8 events
row = 0
col = 0
count = 0
fitted_b8    = []
error_fit_b8 = []
for i_energy_bin in range(len(asimov_dataset)):
    if count == 2:
        row = 1
        col = 0
    

    # do the chi2 minimisation
    chi2 = X.perform_loglikelihood_minimisation(asimov_dataset[i_energy_bin], num_tl208_per_bin[i_energy_bin], values_b8)
    

    # do the multisite minimisation

    # bin and normalise the PDFs
    normed_b8_pdf, _    = np.histogram(pdfs_b8[i_energy_bin], bins = pdf_binning, density = True)
    normed_tl208_pdf, _ = np.histogram(pdfs_tl208[i_energy_bin], bins = pdf_binning, density = True)
    print(len(normed_b8_pdf), len(normed_tl208_pdf))

    # replace any zero count bins in the PDFs with a tiny number
    normed_b8_pdf[normed_b8_pdf == 0] = 1e-6
    normed_tl208_pdf[normed_tl208_pdf ==0] = 1e-6

    logL, poisson = X.perform_multisite_minimisation(values_b8, num_tl208_per_bin[i_energy_bin], dlogL_data[i_energy_bin], normed_b8_pdf, normed_tl208_pdf, pdf_binning, asimov_dataset[i_energy_bin])

    min_chi2_idx    = np.argmin(chi2)
    min_logl_idx    = np.argmin(logL)
    min_poisson_idx = np.argmin(poisson)

    diff_from_zero_chi2    = 0 - chi2[min_chi2_idx]
    diff_from_zero_logl    = 0 - logL[min_logl_idx]
    diff_from_zero_poisson = 0 - poisson[min_poisson_idx]
    chi2    = chi2 + diff_from_zero_chi2
    logL    = logL + diff_from_zero_logl
    poisson = poisson + diff_from_zero_poisson  

    print(min(logL), max(logL))
    print(min(chi2), max(chi2))
    print(min(poisson), max(poisson))

    # plt.figure()
    # plt.plot(values_b8 / (values_b8 + num_tl208_per_bin[i_energy_bin]), full, color = "black", label = "Total")
    # # plt.plot(values_b8 / (values_b8 + num_tl208_per_bin[i_energy_bin]), sig, label = "Signal Term")
    # # plt.plot(values_b8 / (values_b8 + num_tl208_per_bin[i_energy_bin]), back, label = "Background Term")
    # plt.xlabel(r"$F_s$")
    # plt.legend()
    # plt.savefig(f"../plots/asimov_study/understanding/{i_energy_bin}.png")
    # plt.close()

    # calculate the 1 sigma confidence intervals
    sig1 = 0.5
    idx_1sig_full    = np.argmin(np.abs(logL - sig1))
    idx_1sig_poisson = np.argmin(np.abs(poisson - sig1))

    # find the second smallest difference (i.e. the other end of the confidence interval)
    masked_1sig_full    = np.ma.masked_array(logL, mask = False).copy()

    # now remove the smallest value
    masked_1sig_full[idx_1sig_full] = True # True means value is masked out and not considered in subsequent operations

    # repeat finding the smallest difference
    idx2_1sig_full = np.argmin(np.abs(masked_1sig_full-sig1))

    # repeat for the poisson likelihood
    masked_1sig_poisson = np.ma.masked_array(poisson, mask = False).copy()
    masked_1sig_poisson[idx_1sig_poisson] = True
    idx2_1sig_poisson = np.argmin(np.abs(masked_1sig_poisson-sig1))

    signal_1sig_full     = values_b8[idx_1sig_full]
    signal_1sig_full2    = values_b8[idx2_1sig_full]
    signal_1sig_poisson  = values_b8[idx_1sig_poisson]
    signal_1sig_poisson2 = values_b8[idx2_1sig_poisson]

    # create the output plot
    print(row, col, count)
    axes[row, col].vlines(x = signal_1sig_full, ymin = 0, ymax = 0.5, color = "black", linestyle = "dashed")
    axes[row, col].vlines(x = signal_1sig_full2, ymin = 0, ymax = sig1, color = "black", linestyle = "dashed")
    axes[row, col].vlines(x = signal_1sig_poisson, ymin = 0, ymax = sig1, color = "orange", linestyle = "dashed")
    axes[row, col].vlines(x = signal_1sig_poisson2, ymin = 0, ymax = sig1, color = "orange", linestyle = "dashed")
    axes[row, col].plot(values_b8, chi2, label = r"$\chi^2$")
    axes[row, col].plot(values_b8, logL, color = "black")
    axes[row, col].plot(values_b8, poisson, color = "orange")
    axes[row, col].plot([], [], color = "black", label = r"Full $-log(\mathcal{L})$")
    axes[row, col].plot([], [], color = "orange", label = r"Poisson $-log(\mathcal{L})$")
    axes[row, col].set_xlabel("Number of B8", fontsize = 20)
    axes[row, col].axhline(0.5, color = "red", linestyle = "dotted", label = r"$1 \sigma$")
    # axes[row, col].set_ylabel(r"$\chi ^2$", fontsize = 20)
    axes[row, col].set_title(f"{energy_binning[i_energy_bin]}" + r"$\rightarrow$" + f"{energy_binning[i_energy_bin + 1]} MeV", fontsize = 20)
    axes[row, col].axvline(num_b8_per_bin[i_energy_bin], color = "red", label = f"True Num. B8: {round(num_b8_per_bin[i_energy_bin], 2)}")
    axes[row, col].plot([], [], linestyle = "", label = r"$\chi ^2$" + f" Fitted Num. B8: {round(values_b8[np.argmin(chi2)], 2)}\n" + r"Full $-log(\mathcal{L})$" + f" Fitted Num. B8: {round(values_b8[np.argmin(logL)], 2)}" + r" $\pm$ " + f"{np.abs(signal_1sig_full-signal_1sig_full2)/2}\nPoisson " + r"$-log(\mathcal{L})$ Fitted Num. B8: " + f"{round(values_b8[np.argmin(poisson)], 2)}" + r" $\pm$ " + f"{np.abs(signal_1sig_poisson-signal_1sig_poisson2)/2}")
    # axes[i_energy_bin].plot(values_b8, logL, color = "orange")

    axes[row, col].legend(frameon = True, facecolor = "lightgrey", loc = "upper left", fontsize = 10)
    # set some axis boundaries to visualise likelihood around the true minima
    # find idx of likelihood functions corresponding to the true value of B8 in that bin
    print(num_b8_per_bin[i_energy_bin])
    idx_x  = np.where(values_b8 == num_b8_per_bin[i_energy_bin])[0]
    length = 300
    print(idx_x)
    print(values_b8[idx_x-length], values_b8[idx_x + length])
    axes[row,col].set_xlim((values_b8[idx_x-length], values_b8[idx_x+length]))

    # cut the y-axis off at 3 sigma
    axes[row, col].set_ylim(0, 3**2/2)
    axes[row, col].tick_params(axis = "x", labelsize = 15)
    axes[row, col].tick_params(axis = "y", labelsize = 15)
    
    # find the minima of the chi2 and the logL and scale them both so the minima is at 0
    
    # ax2 = axes[row, col].twinx()
    # axes.plot(values_b8, chi2, label = r"$\chi^2$")
    
    # axes.plot([], [], color = "orange", label = r"$-log(\mathcal{L})$")
    # axes.set_xlabel("Number of B8", fontsize = 20)
    # axes.set_ylabel(r"$\chi ^2$", fontsize = 20)
    # axes.set_title(f"{energy_binning[i_energy_bin]}" + r"$\rightarrow$" + f"{energy_binning[i_energy_bin + 1]} MeV", fontsize = 20)
    # axes.axvline(num_b8_per_bin[i_energy_bin], color = "red", label = f"True Num. B8: {round(num_b8_per_bin[i_energy_bin], 2)}")
    # axes.plot([], [], linestyle = "", label = r"$\chi ^2$" + f" Fitted Num. B8: {round(values_b8[np.argmin(chi2)], 2)}\n" + r"$-log(\mathcal{L})$" + f" Fitted Num. B8: {round(values_b8[np.argmin(logL)], 2)}")
    
    
    
    # axes.legend(frameon = False, loc = "upper left")

    # ax2 = axes.twinx()
    
    axes[row, col].set_ylabel(r"$-log(\mathcal{L})$", fontsize = 20)
    count +=1
    col += 1
    
    fitted_b8.append(values_b8[np.argmin(logL)])
    error_fit_b8.append(np.abs(signal_1sig_full-signal_1sig_full2)/2)
    # axes.set_xlim((0,1000))
    # axes.set_ylim((-100, 200))
    # ax2.set_ylim((-7400, -7300))

fig.tight_layout()
plt.savefig("../plots/asimov_study/minimisation_both.pdf")
plt.savefig("../plots/asimov_study/minimisation_both.png")
plt.close()

# create a picture of the asimov dataset with the fitted model and the observed data
plt.figure(figsize = (12, 6))
print(len(energy_binning[0:-1]), len(asimov_dataset[:-1]))
plt.errorbar(energy_binning[0:-1] + np.diff(energy_binning)[0]/2, asimov_dataset, yerr = np.sqrt(asimov_dataset), xerr = np.diff(energy_binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "black", label = "Observed Counts")
# plt.errorbar(energy_binning[0:-1] + np.diff(energy_binning)[0]/2, num_tl208_per_bin, yerr = np.sqrt(num_tl208_per_bin), xerr = np.diff(energy_binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "red", label = "Tl208")
# plt.errorbar(energy_binning[0:-1] + np.diff(energy_binning)[0]/2, fitted_b8, yerr = error_fit_b8, xerr = np.diff(energy_binning)[0]/2, capsize = 2, linestyle = "", marker = "o", color = "green", label = "Fitted B8")
plt.errorbar(energy_binning[0:-1] + np.diff(energy_binning)[0]/2, fitted_b8+num_tl208_per_bin, xerr = 0.25, color = "orange", label = "Total Model")
plt.xlim((2.5, 4.5))
plt.ylim((0, np.max(asimov_dataset)* 1.1))
plt.legend()
plt.xlabel("Reconstructed Energy (MeV)", fontsize = 20)
plt.ylabel(f"Counts per {np.diff(energy_binning)[0]} MeV Bin", fontsize = 20)
plt.title("Asimov Dataset", fontsize = 20)
plt.yscale("log")
plt.savefig("../plots/asimov_study/asimov_dataset_fitted_onlyModelData_log.png")
plt.close()
