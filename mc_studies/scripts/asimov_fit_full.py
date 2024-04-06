import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.optimize
import ROOT

"""
This code is written to perform a likelihood fit from an Asimov dataset, using real MC.

The full background model in 4.5 m in ROI 2.5 --> 5.0 MeV is used as the Asimov 
normalisations for each of the following backgrounds:

    1. Bi214
    2. Bi212
    3. Tl210
    4. Tl208
    5. B8 (nue + nu mu)

    (Alpha, N) fixed at zero counts.

The background normalisations are constrained to vary around the background model
expectation, +- the uncertainty propagated via Poisson statistics from number of
tagged BiPo214/212.

I am writing the optimisation code in python, using scipy.minimize.

The output of the script should be:
    1. 2D -2log(L) heatmap of Nsig vs Ni for each normalisation
    2. Full 5D minimization result (just a number and uncertainty)

For: multisite ONLY, energy PDFs ONLY, multisite + energy PDFs.
"""

class Ahab():
    """
    “...to the last I grapple with thee; from hell's heart I stab at thee; for hate's sake I spit my last breath at thee.” 
    """

    def evaluate_loglikelihood(self, pdfs, normalisations, model_expectations, model_uncertainty, events_per_bin):
        """
        Evaluate the binned extended log-likelihood (BELL) for a given set of PDFs, normalisations, background model
        predictions and model uncertainties (gaussian constraints).

        We assume the data has been binned according to PDF schema previously, and we evaluate the BELL for a given
        set of normalisations. It's all 'static' in this function --> a single constant comes out.

        pdfs              : float array, (N_PDFS, N_bins), each row is the normalised PDF counts for each input PDF
        normalisations    : float array, (N_PDFS), each float is normalisation applied to ith PDF 
        model_expectations: float array, (N_PDFS), float expected number of events for each type from background model [penalty term]
        model_uncertainty : float array, (N_PDFS), poisson statistics (from bipo tagging) constraint on each normalisation [penalty term]
        events_per_bin    : int array, [N_bins], the data counts binned according to the same schema as the PDFs used  
        """

        # norm sum is a constant but for some reason it's included ...
        # should always equal the total number of expected events ...
        norm_sum = np.sum(normalisations)

        """
        For the full derivation of this FULLY OPTIMIZED method, see fancy journal dated 3rd April 2024.
        """

        # multiply the counts of each PDF by the appropriate normalisation factor
        normalisations = normalisations[:, None] # stretch / copy the normalisation factors so it matches number of cols in PDFs
        pdf_norm       = normalisations * pdfs   # elementwise multiplication results in (N_pdfs, N_bins)

        # sum each column to get the inside of each log term contributing to outer sum over bins
        inner_log      = np.sum(pdf_norm, axis = 0) # remove the rows so axis = 0
        
        log            = np.log(inner_log)          # 1D array of where each element is log(Nj Pj) for every bin in PDF

        # multiply each element by the number of data events in each PDF bin and sum all the terms
        scaled_log     = np.sum(events_per_bin * log)

        return norm_sum - scaled_log
        # return -scaled_log
    
    def obtain_pdf(self, location, fv_cut, multisite_bins, energy_bins, run_list, energy_range):
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

        # bin extracted data to create unnormalised PDFs
        energy_counts,    _ = np.histogram(energies, bins = energy_bins)
        multisite_counts, _ = np.histogram(multisites, bins = multisite_bins)

        # pad the zero bins in the PDFs
        energy_counts    = energy_counts.astype(np.float32)    # np arrays contain 1 type of data and counts is int
        multisite_counts = multisite_counts.astype(np.float32) # need to set to float if want to pad with tiny values
        energy_counts[energy_counts == 0]       = 1e-6
        multisite_counts[multisite_counts == 0] = 1e-6

        print(f"Number of counts in energy and multisite PDFs, respectively: {np.sum(energy_counts)}, {np.sum(multisite_counts)}")
        # normalise histogram counts so sum of counts = 1
        energy_counts    = energy_counts    / np.sum(energy_counts)
        multisite_counts = multisite_counts / np.sum(multisite_counts)

        return energy_counts, multisite_counts 

    def rescale_ll(self, ll):
        
        min_idx   = np.argmin(ll)
        diff_zero = 0 - np.ravel(ll)[min_idx]
        ll        = ll + diff_zero

        return ll 

    def grid_search(self, hypothesis_sig, hypothesis_backg, multisite_pdf_array, energy_pdf_array, asimov_multisite, asimov_energy):
        """
        Function performs a 2D grid search over the -2log(L) space for a given
        set of 2 event processes. All other event classes are kept fixed at 
        their background model expectated number.
        """

        # 2D arrays to keep hold of LL calculated at each point
        ll_multisite = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))
        ll_energy    = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))
        ll_full      = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))

        # 2D grid search over LL space
        for isig in range(len(hypothesis_sig)):
            for ibackg in range(len(hypothesis_backg)):
                
                # evaluate 2 x log-likelihood from energy and multisite PDF shapes
                normalisations = np.array([hypothesis_sig[isig], hypothesis_backg[ibackg]])
                multi  = 2 * self.evaluate_loglikelihood(multisite_pdf_array, normalisations, 0, 0, asimov_multisite)
                energy = 2 * self.evaluate_loglikelihood(energy_pdf_array, normalisations, 0, 0, asimov_energy) 
                full   = multi + energy

                # save result for this signal hypothesis
                ll_multisite[isig, ibackg] = multi
                ll_energy[isig, ibackg]    = energy
                ll_full[isig, ibackg]      = full

        return ll_multisite, ll_energy, ll_full

    def minimisation(self):
        pass

    def create_plots(self, sig_hypothesis, backg_hypothesis, ll_full, ll_multi, ll_energy):
        """
        Function creates all the output plots. 
        ---> to be completed once minimisation code is working.
        """
        
        # for making the colourbar the right size
        im_ratio = ll_full.shape[1] / ll_full.shape[0]        
        
        # reverse order of columns so number of B8 increases from bottom left of array upwards
        sig_hypothesis = sig_hypothesis[::-1]
        ll_full        = np.flipud(ll_full)
        ll_multi       = np.flipud(ll_multi)
        ll_energy      = np.flipud(ll_energy)

        # find the minimum value of each likelihood array
        flat_min_idx_full       = np.argmin(ll_full)
        reshaped_min_idx_full   = np.unravel_index(flat_min_idx_full, ll_full.shape)
        flat_min_idx_multi      = np.argmin(ll_multi)
        reshaped_min_idx_multi  = np.unravel_index(flat_min_idx_multi, ll_multi.shape)
        flat_min_idx_energy     = np.argmin(ll_full)
        reshaped_min_idx_energy = np.unravel_index(flat_min_idx_energy, ll_energy.shape)

        # work out the 1 sigma frequentist contour region
        delta_logl = 1 # the variance from minimum for 2log(l) is 1
        
        # find all the points which differ from the minimum by less than or equal to 1
        contour_points = np.where(ll_full <= delta_logl)
        full_vals      = ll_full <= delta_logl
        multi_vals     = ll_multi <= delta_logl 
        energy_vals    = ll_energy <= delta_logl
        # create the plot
        fig, axes = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 10))
        
        img = axes[0].imshow(ll_full, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[0], fraction=0.072*im_ratio)
        axes[0].scatter(backg_hypothesis[reshaped_min_idx_full[1]], sig_hypothesis[reshaped_min_idx_full[0]],  color = "blue", label = f"Min: {sig_hypothesis[reshaped_min_idx_full[0]]} ")
        axes[0].contour(backg_hypothesis, sig_hypothesis, ll_full, levels = [1], colors = "red")
        axes[0].set_title("Full Log-Likelihood")
        axes[0].set_ylabel(r"$N_{^8B}$")
        axes[0].set_xlabel(r"$N_{^{208}Tl}$")
        axes[0].legend()
        
        img = axes[1].imshow(ll_multi, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[1], fraction=0.072*im_ratio)
        axes[1].scatter(backg_hypothesis[reshaped_min_idx_multi[1]], sig_hypothesis[reshaped_min_idx_multi[0]], color = "blue", label = f"Min: {sig_hypothesis[reshaped_min_idx_multi[0]]} ")
        axes[1].contour(backg_hypothesis, sig_hypothesis, ll_multi, levels = [1], colors = "red")
        axes[1].set_title("Multisite Log-Likelihood")
        axes[1].set_ylabel(r"$N_{^8B}$")
        axes[1].set_xlabel(r"$N_{^{208}Tl}$")
        axes[1].legend()
        
        img = axes[2].imshow(ll_energy, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[2], fraction=0.072*im_ratio)
        axes[2].scatter(backg_hypothesis[reshaped_min_idx_energy[1]], sig_hypothesis[reshaped_min_idx_energy[0]], color = "blue", label = f"Min: {sig_hypothesis[reshaped_min_idx_energy[0]]} ")
        axes[2].contour(backg_hypothesis, sig_hypothesis, ll_energy, levels = [1], colors = "red")
        axes[2].set_title("Energy Log-Likelihood")
        axes[2].set_ylabel(r"$N_{^8B}$")
        axes[2].set_xlabel(r"$N_{^{208}Tl}$")
        axes[2].legend()

        fig.tight_layout()
        plt.savefig("../plots/asimov_study/real_mc/advanced/test_2d.png")
        plt.close()

    def run_analysis(self, expected_signal, expected_backg, analyse_real_data, data_path = "", fv_cut = 4500.0, energy_range = "2p5_5p0"):
        """
        Function that actually calls all the other functions and gets
        the analysis done.

        INPUTS: analyse_real_data, bool, if TRUE - do not create Asimov dataset
                
                data_path, str, path to dataset directory - only used if 
                analyse_real_data set to TRUE

                expected_signal, float, expected number of B8 solar events

                expected_backg, float list, expected number of each background 
                considered. list MUST be 1D (4,):

                    [Tl208, Tl210, Bi212, Bi214]

                Pass None to not include normalisation in the fit.

                fv_cut, float, fiducial volume cut to apply when generating PDFs
                and running analyses in mm.

                energy_range: str, must match the name of TTree branch names in 
                dataset ntuples. Determines multisite PDFs used for analysis.
        """
        
        # run checks on the input arguments
        if type(expected_signal) != float or expected_signal < 0:
            print("Expected signal must be a positive float.")
            return 0
        elif len(expected_backg) != 4:
            print("Expected background must be a length 4 array (one entry per background):\n\
                  [Tl208, Tl210, Bi212, Bi214].\n Pass 'None' to exclude normalisation from fit.")
            return 0
        elif expected_backg[0] == None:
            print("Cannot not-include the Tl208 normalisation --> This is the dominant background.")
            return 0
        elif type(analyse_real_data) != bool:
            print("Must enter boolean flag (True/False) whether analysis is run on real data or Asimov MC.")
            return 0
        elif analyse_real_data == True and data_path == "":
            print("Specified to analyse real data! Must enter a path to the dataset-containing folder.")
            return 0
        elif energy_range != "2p5_3p0" or energy_range != "3p0_3p5" or energy_range != "3p5_4p0"\
        or energy_range != "4p0_4p5" or energy_range != "4p5_5p0" or energy_range != "2p5_5p0":
            print("Energy range not recognised. Has previous steps been run for this range?\
                  If So, need to update Ahab to allow this value.")
            return 0
        
        # count the number of backgrounds included in this fit
        num_norms         = sum(1 for item in expected_backg if item != None)
        print("Found ", num_norms, " backgrounds to include in the fit.")

        # binning for the energy and multisite discriminant PDFs
        energy_bins    = np.arange(2.5, 5.1, 0.1)
        multisite_bins = np.arange(-1.4, -1.2, 0.0005)

        """
        Extract information and create the binned PDFS for energy shape and multisite,
        for each isotope of interest.
        """
        # we will create an Asimov dataset using the multisite and energy PDFs of every included normalisation
        pdf_runlist    = np.loadtxt("../runlists/test_runlist.txt", dtype = int)
        mc_path        = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
        sig_mc_path    = f"{mc_path}/full_analysis_B8_solar_nue" # signal path is always the same
        
        # loop over the included normalisations
        print("Creating PDFs from signal and backgrounds in energy domain: ", energy_range)
        backg_names       = ["Tl208", "Tl208", "Bi212", "Bi214"]
        
        # create array containing the signal and background pdfs for multisite and energy
        multisite_pdf_array = np.zeros((num_norms + 1, multisite_bins.size - 1))  # (number backgrounds + 1, number of bins)
        energy_pdf_array    = np.zeros((num_norms + 1 , energy_bins.size - 1))    # (number backgrounds + 1, number of bins)

        # obtain the signal PDF
        sig_energy_pdf, sig_multi_pdf = self.obtain_pdf(sig_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist, energy_range)

        # add the signal pdfs as the first row in the pdfs arrays
        multisite_pdf_array[0, :] = sig_multi_pdf
        energy_pdf_array[0, :]    = sig_energy_pdf

        # obtain background pdfs and add to subsequent rows
        idx_counter = 1
        for iname in range(len(backg_names)):
            if expected_backg[iname] != None:
                
                # we have an expected rate and want to include this normalisation in the Asimov dataset
                backg_mc_path  = f"{mc_path}/full_analysis_{backg_names[iname]}"

                # obtain the normalised PDFs for this background and add to the respective arrays
                backg_energy_pdf, backg_multi_pdf   = self.obtain_pdf(backg_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist, energy_range)
                multisite_pdf_array[idx_counter, :] = backg_multi_pdf
                energy_pdf_array[idx_counter, :]    = backg_energy_pdf

                # incremement idx counter ... not the same idx as the number of backg names!
                idx_counter += 1 
         
        # decide whether to create an Asimov dataset out of the signal + background PDFs, or load in real data
        if analyse_real_data == False:
            print("Creating Asimov dataset from normalisation-scaled PDFs.")

            # add the signal normalisation and to the front of the backg norm list
            expected_backg = np.array(expected_backg)
            expected_backg = expected_backg[expected_backg != None].tolist()
            normalisations = [expected_signal] + expected_backg

            # create the asimov dataset - scale each PDF by the expected rate and sum
            normalisations = normalisations[:, None]
            
            dataset_energy    = normalisations * energy_pdf_array    # multiply every bin by corresponding normalisation
            dataset_multisite = normalisations * multisite_pdf_array
            dataset_energy    = np.sum(dataset_energy, axis = 0)     # remove the rows so axis = 0
            dataset_multisite = np.sum(dataset_multisite, axis = 0)

        else:

            # open the dataset --> saved as a single hadded TTree containing ITR,
            # multisite, position and energy of each extracted event
            datafile = ROOT.TFile.Open(data_path)

            # extract the correct energy range dataset --> load the correct tree
            datatree = datafile.Get(energy_range)

            # get number of entries in tree to setup numpy dataset arrays
            num_entries       = datatree.GetEntries()
            dataset_energy    = np.zeros(num_entries)
            dataset_multisite = np.zeros(num_entries)
            
            # fill the dataset arrays
            entry_count = 0
            for ientry in datatree:
                
                # apply the fv cut
                x = ientry.x
                y = ientry.y
                z = ientry.z

                r = np.sqrt(x**2 + y**2 + z**2)
                if r > fv_cut:
                    continue

                dataset_energy[entry_count]    = ientry.energy
                dataset_multisite[entry_count] = ientry.dlogL

            # now we need to create the binned dataset
            dataset_energy, _    = np.histogram(dataset_energy, bins = energy_bins)
            dataset_multisite, _ = np.histogram(dataset_multisite, bins = multisite_bins)

        """
        Binned enegy / multisite discriminant datasets obtained. Now we run the fits.
        """

        # perform 2D grid searches over each background and the signal
        # grid search around the expected values given for each isotope and signal
        for iname in range(len(backg_names)):
            if expected_backg[iname] != None:

                # include this normalisation and perform 2 D scan
                # signal and backg hypothesis are expected val +- 50 %
                signal_range = expected_signal       * 0.5
                backg_range  = expected_backg[iname] * 0.5 

                sig_hypothesis   = np.arange(expected_signal - signal_range, expected_signal + signal_range + 1, 1)
                backg_hypothesis = np.arange(expected_backg[iname] - backg_range, expected_backg[iname] _ backg_range + 1, 1)

                # perform 2D scan over given signal and background normalisation
                ll_multisite, ll_energy, ll_full = self.grid_search(sig_hypothesis, backg_hypothesis, multisite_pdf_array, energy_pdf_array, dataset_multisite, dataset_energy)
        
        # rescale all the log-likelihood functions so minimum value is at zero --> for plotting & comparisons
        ll_full = self.rescale_ll(ll_full)
        ll_multisite = self.rescale_ll(ll_multisite)
        ll_energy = self.rescale_ll(ll_energy)
        
        # create plots
        self.create_plots(sig_hypothesis, backg_hypothesis, ll_full, ll_multisite, ll_energy)

Ahab().run_analysis()