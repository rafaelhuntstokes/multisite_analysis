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

    def evaluate_loglikelihood(self, normalisations, pdfs, model_expectations, model_uncertainty, events_per_bin):
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
        
        # we allow the total normalisation to float so see what it is
        norm_sum = np.sum(normalisations)

        # have to do this inside the function for scipy.optimize.minimize doesn't
        # respect the reshaping being done outside its own calls ... 
        normalisations = normalisations[:, None]

        """
        For the full derivation of this FULLY OPTIMIZED method, see fancy journal dated 3rd April 2024.
        """

        # multiply the counts of each PDF by the appropriate normalisation factor
        pdf_norm       = normalisations * pdfs   # elementwise multiplication results in (N_pdfs, N_bins)

        # sum each column to get the inside of each log term contributing to outer sum over bins
        inner_log      = np.sum(pdf_norm, axis = 0) # remove the rows so axis = 0
        
        log            = np.log(inner_log)          # 1D array of where each element is log(Nj Pj) for every bin in PDF

        # multiply each element by the number of data events in each PDF bin and sum all the terms
        scaled_log     = np.sum(events_per_bin * log)

        return 2 * (norm_sum - scaled_log) # convert to -2 x log(L)
        # return -scaled_log
    
    def combined_loglikelihood(self, normalisations, energy_pdfs, multisite_pdfs, model_expectations, model_uncertainty, counts_energy, counts_multisite):
        """
        Combined log-likelihood evaluation of energy and multisite discriminant.
        """

        # we allow the total normalisation to float so see what it is
        norm_sum = np.sum(normalisations)

        """
        For the full derivation of this FULLY OPTIMIZED method, see fancy journal dated 3rd April 2024.
        """

        # have to do this inside the function for scipy.optimize.minimize doesn't
        # respect the reshaping being done outside its own calls ... 
        normalisations = normalisations[:, None]

        # multiply the counts of each PDF by the appropriate normalisation factor
        energy_pdf_norm       = normalisations * energy_pdfs         # elementwise multiplication results in (N_pdfs, N_bins)
        multisite_pdf_norm    = normalisations * multisite_pdfs
        
        # sum each column to get the inside of each log term contributing to outer sum over bins
        energy_inner_log      = np.sum(energy_pdf_norm, axis = 0)    # remove the rows so axis = 0
        multisite_inner_log   = np.sum(multisite_pdf_norm, axis = 0)
        
        energy_log            = np.log(energy_inner_log)             # 1D array of where each element is log(Nj Pj) for every bin in PDF
        multisite_log         = np.log(multisite_inner_log)
        
        # multiply each element by the number of data events in each PDF bin and sum all the terms
        energy_scaled_log     = np.sum(counts_energy * energy_log)
        multisite_scaled_log  = np.sum(counts_multisite * multisite_log)

        energy_component      = 2 * (norm_sum - energy_scaled_log)
        multisite_component   = 2 * (norm_sum - multisite_scaled_log)

        # combine the information
        full_loglikelihood    = energy_component + multisite_component
        return full_loglikelihood
    
    def obtain_pdf(self, location, fv_cut, multisite_bins, energy_bins, run_list, energy_range, plot_name):
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
        energy_counts    = energy_counts.astype(np.float64)    # np arrays contain 1 type of data and counts is int
        multisite_counts = multisite_counts.astype(np.float64) # need to set to float if want to pad with tiny values
        energy_counts[energy_counts == 0]       = 1e-6
        multisite_counts[multisite_counts == 0] = 1e-6

        print(f"Number of counts in energy and multisite PDFs, respectively: {np.sum(energy_counts)}, {np.sum(multisite_counts)}")
        # normalise histogram counts so sum of counts = 1
        energy_counts    = energy_counts    / np.sum(energy_counts)
        multisite_counts = multisite_counts / np.sum(multisite_counts)

        ## plot the pdfs ##
        fig, axes = plt.subplots(nrows = 1, ncols = 2)
        axes[0].hist(energies, bins = energy_bins, histtype= "step", color = "black", linewidth = 2)
        axes[0].set_xlabel("Reconstructed Energy (MeV)")
        axes[0].set_ylabel(f"Counts per {round(np.diff(energy_bins)[0], 2)} MeV Bin")
        axes[0].set_title(f"Energy PDF for {plot_name}")

        axes[1].hist(multisites, bins = multisite_bins, histtype= "step", color = "black", linewidth = 2)
        axes[1].set_xlabel("Multisite Discriminant")
        axes[1].set_ylabel(f"Counts per {round(np.diff(multisite_bins)[0], 2)} MeV Bin")
        axes[1].set_title(f"Multisite PDF for {plot_name}")
        fig.tight_layout()
        plt.savefig(f"../plots/asimov_study/real_mc/advanced/pdf_{plot_name}.png")
        plt.close()

        return energy_counts, multisite_counts 

    def rescale_ll(self, ll):
        
        min_idx   = np.argmin(ll)
        diff_zero = 0 - np.ravel(ll)[min_idx]
        ll        = ll + diff_zero

        return ll 

    def grid_search(self, hypothesis_sig, hypothesis_backg, normalisations, backg_idx, multisite_pdf_array, energy_pdf_array, dataset_multisite, dataset_energy):
        """
        Function performs a 2D grid search over the -2log(L) space for a given
        set of 2 event processes. All other event classes are kept fixed at 
        their background model expectated number.
        """

        # copy the normalisation array to avoid updating the main array
        normalisations = normalisations.copy()

        # 2D arrays to keep hold of LL calculated at each point
        ll_multisite = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))
        ll_energy    = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))
        ll_full      = np.zeros((len(hypothesis_sig), len(hypothesis_backg)))

        # 2D grid search over LL space
        for isig in range(len(hypothesis_sig)):
            for ibackg in range(len(hypothesis_backg)):
                
                # update the normalisations at the backg_idx and signal_idx (rest kept constant)
                normalisations[0]         = hypothesis_sig[isig]
                normalisations[backg_idx] = hypothesis_backg[ibackg]

                multi  = self.evaluate_loglikelihood(normalisations, multisite_pdf_array, 0, 0, dataset_multisite)
                energy = self.evaluate_loglikelihood(normalisations, energy_pdf_array, 0, 0, dataset_energy) 
                full   = multi + energy

                # save result for this signal hypothesis
                ll_multisite[isig, ibackg] = multi
                ll_energy[isig, ibackg]    = energy
                ll_full[isig, ibackg]      = full

        return ll_multisite, ll_energy, ll_full

    def minimisation(self, normalisations, energy_pdfs, multisite_pdfs, dataset_multisite, dataset_energy):
        """
        Function performs a minimisation optimisation by varying the normalisation
        of the signal and a set of backgrounds, constrained to lie within
        some confidence bound.

        INPUTS: 
            normalisations, float array, array of normalisations to optimise set at their
                            expected values

        RETURNS:
            energy_result, multisite_result, comined_result : objects containing best fit norms (.x) and termination flag .success (bool)
        """
        normalisations = normalisations.copy()

        # first run optimisation using energy PDFs - returns OptimizeResult object containing best fit norms
        print("Performing minimisation using energy PDF information ...")
        energy_result    = scipy.optimize.minimize(self.evaluate_loglikelihood, x0 = normalisations, method = "BFGS", tol = 1e-4, args = (energy_pdfs, 0, 0, dataset_energy,))
        
        print("Performing minimisation using multisite PDF information ...")
        multisite_result = scipy.optimize.minimize(self.evaluate_loglikelihood, x0 = normalisations, method = "BFGS", tol = 1e-4, args = (multisite_pdfs, 0, 0, dataset_multisite,))
        
        # run minimisation using combination of energy and multisite PDFs
        print("Performing minimisation using multisite + energy PDF information ...")
        combined_result  = scipy.optimize.minimize(self.combined_loglikelihood, x0 = normalisations, method = "BFGS", tol = 1e-4, args = (energy_pdfs, multisite_pdfs, 0, 0, dataset_energy, dataset_multisite,))

        return energy_result, multisite_result, combined_result
    
    def profile_likelihood_scan(self, sig_hypothesis, normalisations, energy_pdfs, multisite_pdfs, dataset_multisite, dataset_energy):
        """
        Function performs a profile likelihood scan: grid search over B8 signal
        and fit of the other normalisations performed at every point.

        Overall minimum log-likelihood returned at each fixed B8 signal norm.
        """

        def wrapped_combined(other_norms, sig_norm, energy_pdfs, multisite_pdfs, dataset_energy, dataset_multisite):
            """
            Function handles the fact now we fix norms[0] and vary the others.
            """

            # insert the fixed signal norm into the norms array
            norms = [sig_norm] + other_norms.tolist()
            norms = np.array(norms)

            return self.combined_loglikelihood(norms, energy_pdfs, multisite_pdfs, 0, 0, dataset_energy, dataset_multisite)
        
        def wrapped_multisite(other_norms, sig_norm, multisite_pdfs, dataset_multisite):
            """
            Function handles the fact now we fix norms[0] and vary the others.
            """

            # insert the fixed signal norm into the norms array
            norms = [sig_norm] + other_norms.tolist()
            norms = np.array(norms)

            return self.evaluate_loglikelihood(norms, multisite_pdfs, 0, 0, dataset_multisite)

        def wrapped_energy(other_norms, sig_norm, energy_pdfs, dataset_energy):
            """
            Function handles the fact now we fix norms[0] and vary the others.
            """

            # insert the fixed signal norm into the norms array
            norms = [sig_norm] + other_norms.tolist()
            norms = np.array(norms)

            return self.evaluate_loglikelihood(norms, energy_pdfs, 0, 0, dataset_energy)
        
        norms = normalisations.copy() # don't update / mess with the input normalisations
        
        # loop over each signal hypothesis
        profile_ll      = np.zeros((3, len(sig_hypothesis)))             # (num. optimisations, num. signal vals)
        optimised_norms = np.zeros((3, len(sig_hypothesis), len(norms))) # (num. optimisations, num. sig vals, num. norms)

        for isig in range(len(sig_hypothesis)):

            # set the norm of the signal
            signal_norm = sig_hypothesis[isig]

            # define the optimisation to be run
            combined_result  = scipy.optimize.minimize(wrapped_combined, x0 = norms[1:], bounds = [(0, np.inf)] * len(norms[1:]), method = "L-BFGS-B", tol = 1e-4, args = (signal_norm, energy_pdfs, multisite_pdfs, dataset_energy, dataset_multisite,))
            multisite_result = scipy.optimize.minimize(wrapped_multisite, x0 = norms[1:], bounds = [(0, np.inf)] * len(norms[1:]), method = "L-BFGS-B", tol = 1e-4, args = (signal_norm, multisite_pdfs, dataset_multisite,))
            energy_result    = scipy.optimize.minimize(wrapped_energy, x0 = norms[1:], bounds = [(0, np.inf)] * len(norms[1:]), method = "L-BFGS-B", tol = 1e-4, args = (signal_norm, energy_pdfs, dataset_energy,))

            # save the minimum ll for this minimisation
            profile_ll[0, isig]   = combined_result.fun
            profile_ll[1, isig]   = multisite_result.fun
            profile_ll[2, isig]   = energy_result.fun

            # save the set of normalsations for this minimisation
            optimised_norms[0, isig, :] = [signal_norm] + combined_result.x.tolist()
            optimised_norms[1, isig, :] = [signal_norm] + multisite_result.x.tolist()
            optimised_norms[2, isig, :] = [signal_norm] + energy_result.x.tolist()
        
        return profile_ll, optimised_norms

    def create_gridsearch_plots(self, sig_hypothesis, backg_hypothesis, name_idx, names, ll_full, ll_multi, ll_energy):
        """
        Create plot showing result of pair-wise grid search of signal and background (whilst keeping other
        normalisations fixed to their expected values... suppose this is similar to a contour likelihood...).
        """
        
        # for making the colourbar the right size
        im_ratio = ll_full.shape[1] / ll_full.shape[0]        
        
        # work out what background isotope is being used
        name = names[name_idx]
        if name == "Tl208":
            label_name = r"$N_{^{208}Tl}$"
        if name == "Tl210":
            label_name = r"$N_{^{210}Tl}$"
        if name == "BiPo212":
            label_name = r"$N_{^{212}Bi}$"
        if name == "BiPo214":
            label_name = r"$N_{^{214}Bi}$"

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
        
        # create the plot
        fig, axes = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 10))
        
        img = axes[0].imshow(ll_full, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[0], fraction=0.072*im_ratio)
        axes[0].scatter(backg_hypothesis[reshaped_min_idx_full[1]], sig_hypothesis[reshaped_min_idx_full[0]],  color = "blue", label = f"Min: {round(sig_hypothesis[reshaped_min_idx_full[0]], 2)} ")
        axes[0].contour(backg_hypothesis, sig_hypothesis, ll_full, levels = [1.15], colors = "red")
        axes[0].set_title("Full Log-Likelihood")
        axes[0].set_ylabel(r"$N_{^8B}$")
        axes[0].set_xlabel(label_name)
        axes[0].legend()
        
        img = axes[1].imshow(ll_multi, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[1], fraction=0.072*im_ratio)
        axes[1].scatter(backg_hypothesis[reshaped_min_idx_multi[1]], sig_hypothesis[reshaped_min_idx_multi[0]], color = "blue", label = f"Min: {round(sig_hypothesis[reshaped_min_idx_multi[0]], 2)} ")
        axes[1].contour(backg_hypothesis, sig_hypothesis, ll_multi, levels = [1.15], colors = "red")
        axes[1].set_title("Multisite Log-Likelihood")
        axes[1].set_ylabel(r"$N_{^8B}$")
        axes[1].set_xlabel(label_name)
        axes[1].legend()
        
        img = axes[2].imshow(ll_energy, origin = "lower", aspect = "auto", cmap = "magma", extent = [backg_hypothesis[0], backg_hypothesis[-1], sig_hypothesis[0], sig_hypothesis[-1]])
        # plt.colorbar(img, ax = axes[2], fraction=0.072*im_ratio)
        axes[2].scatter(backg_hypothesis[reshaped_min_idx_energy[1]], sig_hypothesis[reshaped_min_idx_energy[0]], color = "blue", label = f"Min: {round(sig_hypothesis[reshaped_min_idx_energy[0]], 2)} ")
        axes[2].contour(backg_hypothesis, sig_hypothesis, ll_energy, levels = [1.15], colors = "red")
        axes[2].set_title("Energy Log-Likelihood")
        axes[2].set_ylabel(r"$N_{^8B}$")
        axes[2].set_xlabel(label_name)
        axes[2].legend()

        fig.tight_layout()
        plt.savefig(f"../plots/asimov_study/real_mc/advanced/test_2d_{name}.png")
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
        elif (energy_range != "2p5_3p0" and energy_range != "3p0_3p5" and energy_range != "3p5_4p0"\
        and energy_range != "4p0_4p5" and energy_range != "4p5_5p0" and energy_range != "2p5_5p0"):
            print("Energy range not recognised. Has previous steps been run for this range?\nIf So, need to update Ahab to allow this value.")
            return 0
        
        # count the number of backgrounds included in this fit
        num_norms         = sum(1 for item in expected_backg if item != None)
        print("Found ", num_norms, " backgrounds to include in the fit.")

        # binning for the energy and multisite discriminant PDFs
        energy_bins    = np.arange(2.5, 5.1, 0.1)
        multisite_bins = np.arange(-1.375, -1.325, 0.0005)

        """
        Extract information and create the binned PDFS for energy shape and multisite,
        for each isotope of interest.
        """
        # we will create an Asimov dataset using the multisite and energy PDFs of every included normalisation
        pdf_runlist    = np.loadtxt("../runlists/full_test.txt", dtype = int)
        mc_path        = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
        sig_mc_path    = f"{mc_path}/full_analysis_B8_solar_nue" # signal path is always the same
        
        # loop over the included normalisations
        print("Creating PDFs from signal and backgrounds in energy domain: ", energy_range)
        backg_names       = ["Tl208", "Tl210", "BiPo212", "BiPo214"]
        
        # create array containing the signal and background pdfs for multisite and energy
        multisite_pdf_array = np.zeros((num_norms + 1, multisite_bins.size - 1))  # (number backgrounds + 1, number of bins)
        energy_pdf_array    = np.zeros((num_norms + 1 , energy_bins.size - 1))    # (number backgrounds + 1, number of bins)

        # obtain the signal PDF
        print("Obtaining signal PDF.")
        sig_energy_pdf, sig_multi_pdf = self.obtain_pdf(sig_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist, energy_range, "B8_nue")

        # add the signal pdfs as the first row in the pdfs arrays
        multisite_pdf_array[0, :] = sig_multi_pdf
        energy_pdf_array[0, :]    = sig_energy_pdf

        # obtain background pdfs and add to subsequent rows
        idx_counter = 1
        for iname in range(len(backg_names)):
            if expected_backg[iname] != None:
                print("Obtaining background pdf: ", backg_names[iname])
                # we have an expected rate and want to include this normalisation in the Asimov dataset
                backg_mc_path  = f"{mc_path}/full_analysis_{backg_names[iname]}"

                # obtain the normalised PDFs for this background and add to the respective arrays
                backg_energy_pdf, backg_multi_pdf   = self.obtain_pdf(backg_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist, energy_range, backg_names[iname])
                multisite_pdf_array[idx_counter, :] = backg_multi_pdf
                energy_pdf_array[idx_counter, :]    = backg_energy_pdf

                # incremement idx counter ... not the same idx as the number of backg names!
                idx_counter += 1 
        
        # add the signal normalisation and to the front of the backg norm list
        backg_norms    = [element for element in expected_backg if element is not None]
        normalisations = [expected_signal] + backg_norms

        # create the asimov dataset - scale each PDF by the expected rate and sum
        normalisations = np.array(normalisations, dtype = np.float64)
        print("Normalisations Applied: ", normalisations, normalisations.shape)
        
        # stretch / copy the normalisation factors so it matches number of cols in PDFs
        normalisations_stretched = normalisations[:, None]

        # decide whether to create an Asimov dataset out of the signal + background PDFs, or load in real data
        if analyse_real_data == False:
            
            print("Creating Asimov dataset from normalisation-scaled PDFs.")
            dataset_energy    = normalisations_stretched * energy_pdf_array    # multiply every bin by corresponding normalisation
            
            # plot Asimov dataset as a stacked histogram
            labels = ["B8"] + backg_names
            mids_energy = energy_bins[:-1] + np.diff(energy_bins)[0] / 2
            mids_multi = multisite_bins[:-1] + np.diff(multisite_bins)[0] / 2
            fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 8))
            color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
            for i in range(dataset_energy.shape[0]):
                col = next(color_cycle)
                axes[0].step(energy_bins, dataset_energy[i,:].tolist() + [0], where = 'post', color = col, label = f"{labels[i]}| Num: {round(np.sum(dataset_energy[i,:]), 2)}")
            
            # sum 
            dataset_energy    = np.sum(dataset_energy, axis = 0)     # remove the rows so axis = 0
            axes[0].step(energy_bins, dataset_energy.tolist() + [0], where = 'post', color="black", label = f"Sum: {round(np.sum(dataset_energy), 2)}")
            axes[0].legend(frameon=False)
            axes[0].set_xlabel("Reconstructed Energy (MeV)")
            axes[0].set_ylabel(f"Counts per {round(np.diff(energy_bins)[0],2)} MeV")
            axes[0].set_title("Asimov Dataset: Energy")
            axes[0].set_ylim((0, 65))
            axes[0].set_xlim((2.5, 5.0))
            
            # repeat for the multisite discriminant
            color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
            dataset_multisite = normalisations_stretched * multisite_pdf_array
            for i in range(dataset_multisite.shape[0]):
                col = next(color_cycle)
                axes[1].step(multisite_bins, dataset_multisite[i,:].tolist() + [0], where = 'post', color=col, label = f"{labels[i]}| Num: {round(np.sum(dataset_multisite[i,:]), 2)}")
            dataset_multisite = np.sum(dataset_multisite, axis = 0)
            axes[1].step(multisite_bins, dataset_multisite.tolist() + [0], where = 'post', color="black", label = f"Sum: {round(np.sum(dataset_multisite), 2)}")
            axes[1].legend(frameon=False)
            axes[1].set_xlabel("Multisite Discriminant")
            axes[1].set_ylabel(f"Counts per {round(np.diff(multisite_bins)[0],2)}")
            axes[1].set_title("Asimov Dataset: Multisite")
            axes[1].set_xlim((-1.355, -1.335))
            axes[1].set_ylim((0, 65))
            fig.tight_layout()
            plt.savefig("../plots/asimov_study/real_mc/advanced/asimov_dataset_energy.png")
            plt.close()
            
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

            # now we need to create the binned dataset --> just counts in each bin
            dataset_energy, _    = np.histogram(dataset_energy, bins = energy_bins)
            dataset_multisite, _ = np.histogram(dataset_multisite, bins = multisite_bins)

        """
        Binned energy / multisite discriminant datasets obtained. Now we run the fits.
        """

        # perform 2D grid searches over each background and the signal
        # grid search around the expected values given for each isotope and signal
        print("Performing grid search over signal and backg pairs.")

        # not necessarily the same location in normalisation array if masked out some backgrounds
        # starts at 1 because idx 0 of normalisations is the signal
        norm_idx = 1 

        for iname in range(len(backg_names)):
            if expected_backg[iname] != None:
                print("Grid search for backg ", norm_idx, expected_backg[iname])
                # include this normalisation and perform 2 D scan
                # signal and backg hypothesis are expected val +- 50 %
                signal_range = expected_signal       * 0.5
                backg_range  = expected_backg[iname] * 0.5 

                # sig_hypothesis   = np.arange(expected_signal - signal_range, expected_signal + signal_range + 1, 1)
                # backg_hypothesis = np.arange(expected_backg[iname] - backg_range, expected_backg[iname] + backg_range + 1, 1)
                sig_hypothesis   = np.linspace(expected_signal - signal_range, expected_signal + signal_range + 1, 100)
                backg_hypothesis = np.linspace(expected_backg[iname] - backg_range, expected_backg[iname] + backg_range + 1, 100)
                
                print(backg_hypothesis)
                # perform 2D scan over given signal and background normalisation
                ll_multisite, ll_energy, ll_full = self.grid_search(sig_hypothesis, backg_hypothesis, normalisations, norm_idx, multisite_pdf_array, energy_pdf_array, dataset_multisite, dataset_energy)

                # rescale all the log-likelihood functions so minimum value is at zero --> for plotting & comparisons
                ll_full      = self.rescale_ll(ll_full)
                ll_multisite = self.rescale_ll(ll_multisite)
                ll_energy    = self.rescale_ll(ll_energy)
        
                # create plots
                self.create_gridsearch_plots(sig_hypothesis, backg_hypothesis, iname, backg_names, ll_full, ll_multisite, ll_energy)

                # next pair so apply hypothesis to next pair in normalisations array
                norm_idx += 1 

        # now we run the full minimisation scripts for an arbitrary number of backgrounds
        res_energy, res_multisite, res_combined = self.minimisation(normalisations, energy_pdf_array, multisite_pdf_array, dataset_multisite, dataset_energy)

        # calculate the uncertainty by taking inverse of hessian --> covariance matrix and finding sqrt
        # of the diagonal elements
        err_energy    = np.sqrt(np.diag((res_energy.hess_inv)))
        err_multisite = np.sqrt(np.diag((res_multisite.hess_inv)))
        err_combined  = np.sqrt(np.diag((res_combined.hess_inv)))

        print("Best Fit Results:")
        print(f"Energy: {res_energy.x} | Error: {err_energy}\nTerminate status {res_energy.success} | Terminate reason: {res_energy.message}\n")
        print(f"Multisite: {res_multisite.x} | Error: {err_multisite}\nTerminate status {res_multisite.success} | Terminate reason: {res_multisite.message}\n")
        print(f"Combined: {res_combined.x} | Error: {err_combined}\nTerminate status {res_combined.success} | Terminate reason: {res_combined.message}\n")

        # profile log-likelihood scans
        print("Performing profile likelihood scan ... ")

        # row 0 - combined; row 1 multisite; row 2 energy
        profile_ll, optimised_norms = self.profile_likelihood_scan(sig_hypothesis, normalisations, energy_pdf_array, multisite_pdf_array, dataset_multisite, dataset_energy)
        
        # rescale to stick minimum at zero
        profile_ll[0, :] = self.rescale_ll(profile_ll[0, :])
        profile_ll[1, :] = self.rescale_ll(profile_ll[1, :])
        profile_ll[2, :] = self.rescale_ll(profile_ll[2, :])
        
        # find 1 sigma frequentist confidence interval #
        
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
        combined_min   = sig_hypothesis[np.argmin(profile_ll[0,:])]
        combined_upper = abs(combined_min - sig_hypothesis[second_closest[0]]) 
        combined_lower = abs(combined_min - sig_hypothesis[closest_to_interval[0]])
        multi_min      = sig_hypothesis[np.argmin(profile_ll[1,:])]
        multi_upper    = abs(multi_min - sig_hypothesis[second_closest[1]]) 
        multi_lower    = abs(multi_min - sig_hypothesis[closest_to_interval[1]])
        energy_min     = sig_hypothesis[np.argmin(profile_ll[2,:])]
        energy_upper   = abs(energy_min - sig_hypothesis[second_closest[2]]) 
        energy_lower   = abs(energy_min - sig_hypothesis[closest_to_interval[2]])
        
        # create a plot of profile likelihood scan
        plt.figure()
        plt.plot(sig_hypothesis, profile_ll[0,:], color = "black", label = rf"${round(combined_min, 1)}^{{+{round(combined_upper, 1)}}}_{{-{round(combined_lower, 1)}}}$")
        plt.plot(sig_hypothesis, profile_ll[1,:], color = "green", label = rf"${round(multi_min, 1)}^{{+{round(multi_upper, 1)}}}_{{-{round(multi_lower, 1)}}}$")
        plt.plot(sig_hypothesis, profile_ll[2,:], color = "orange", label = rf"${round(energy_min, 1)}^{{+{round(energy_upper, 1)}}}_{{-{round(energy_lower, 1)}}}$")
        
        # confidence intervals
        plt.axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$ frequentist")
        # plt.vlines(x = sig_hypothesis[closest_to_interval[0]], ymin= 0, ymax = 1)
        # plt.vlines(x = sig_hypothesis[second_closest[0]], ymin = 0, ymax = 1)
        # plt.vlines(x = sig_hypothesis[closest_to_interval[1]], ymin = 0, ymax = 1)
        # plt.vlines(x = sig_hypothesis[second_closest[1]], ymin = 0, ymax = 1)
        # plt.vlines(x = sig_hypothesis[closest_to_interval[2]], ymin = 0, ymax = 1)
        # plt.vlines(x = sig_hypothesis[second_closest[2]], ymin = 0, ymax = 1)
        plt.xlabel(r"$N_{^8B}$ Events")
        plt.ylabel(r"$-2log(\mathcal{L})$")
        plt.title("Profile Log-Likelihood Scan")
        plt.ylim((0, 3))
        plt.axvline(x = expected_signal, color = "red", label = f"True Signal Number: {expected_signal}")
        plt.legend(loc = "upper right", fontsize = 11)
        plt.savefig("../plots/asimov_study/real_mc/advanced/profile_ll.png")
        print("All complete!")

        # create output 'fitted model' plot in terms of energy and multisite with profile likelihood result
        # find idx of log-likelihood with minimum value
        minimum_idx = np.argmin(profile_ll, axis = 1)

        # find the normalisations fom energy, multisite and combined LL giving that minimum
        norms_combined = optimised_norms[0, minimum_idx[0], :]
        norms_multi    = optimised_norms[1, minimum_idx[1], :]
        norms_energy   = optimised_norms[2, minimum_idx[2], :]

        print(norms_combined)
        print(norms_multi)
        print(norms_energy)

        # create 3 subplots showing the fitted spectrum using minimum of each method
        fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 15))
        font_s    = 14
        # need to scale the PDFs of each by the normalisations
        combined_energy_result    = norms_combined[:, None] * energy_pdf_array
        combined_multisite_result = norms_combined[:, None] * multisite_pdf_array
        energy_result             = norms_energy[:, None]   * energy_pdf_array
        multisite_result          = norms_multi[:, None]    * multisite_pdf_array

        # create the 4 subplots showing relative contributions of each normalisation, the sum and the data
        color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        for i in range(len(normalisations)):

            # plot each normalisation individually as a bar graph
            col = next(color_cycle)

            if col == '#d62728':
                col = next(color_cycle)
            axes[0,0].step(energy_bins, combined_energy_result[i, :].tolist() + [0], where = 'post', linewidth = 2, color=col, label = f"{labels[i]}| Num: {round(np.sum(combined_energy_result[i,:]), 2)}")
        
        # plot the sum of the model
        sum_combined_energy_model = np.sum(combined_energy_result, axis = 0)
        axes[0,0].step(energy_bins, sum_combined_energy_model.tolist() + [0], where = 'post', color = "red", linewidth = 2, label = f"Total Model: {round(np.sum(sum_combined_energy_model), 2)}")
        axes[0,0].set_title(r"Energy Fit: Combined $\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[0,0].set_xlabel("Reconstructed Energy (MeV)", fontsize = font_s)
        axes[0,0].set_ylabel("Counts", fontsize = font_s)
        axes[0,0].set_ylim((0, 65))
        axes[0,0].set_xlim((2.5, 5.0))
        
        # plot the dataset and error bars
        axes[0,0].errorbar(mids_energy, dataset_energy, yerr = np.sqrt(dataset_energy), color = "black", marker = "^", capsize = 2, linestyle = "", label = "Data")
        axes[0,0].legend(frameon = False, fontsize = 11)

        color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        for i in range(len(normalisations)):

            # plot each normalisation individually as a bar graph
            col = next(color_cycle)
            if col == '#d62728':
                col = next(color_cycle)
            axes[0,1].step(multisite_bins, combined_multisite_result[i, :].tolist() + [0], linewidth = 2, where = 'post', color=col, label = f"{labels[i]}| Num: {round(np.sum(combined_multisite_result[i,:]), 2)}")
        
        # plot the sum of the model
        sum_combined_multisite_model = np.sum(combined_multisite_result, axis = 0)
        axes[0,1].step(multisite_bins, sum_combined_multisite_model.tolist() + [0], where = 'post', linewidth = 2, color = "red", label = f"Total Model: {round(np.sum(sum_combined_multisite_model), 2)}")
        axes[0,1].set_title(r"Multisite Fit: Combined $\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[0,1].set_xlabel("Multisite Discriminant", fontsize = font_s)
        axes[0,1].set_ylabel("Counts", fontsize = font_s)
        axes[0,1].set_xlim((-1.355, -1.335))
        axes[0,1].set_ylim((0, 65))

        # plot the dataset and error bars
        axes[0,1].errorbar(mids_multi, dataset_multisite, yerr = np.sqrt(dataset_multisite), color = "black", marker = "^", capsize = 2, linestyle = "", label = "Data")
        axes[0,1].legend(frameon = False, fontsize = 11)

        color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        for i in range(len(normalisations)):

            # plot each normalisation individually as a bar graph
            col = next(color_cycle)
            if col == '#d62728':
                col = next(color_cycle)
            axes[1,0].step(energy_bins, energy_result[i, :].tolist() + [0], linewidth = 2, where = 'post', color=col, label = f"{labels[i]}| Num: {round(np.sum(energy_result[i,:]), 2)}")
        
        # plot the sum of the model
        sum_energy_model = np.sum(energy_result, axis = 0)
        axes[1,0].step(energy_bins, sum_energy_model.tolist() + [0], where = 'post',  linewidth = 2, color = "red", label = f"Total Model: {round(np.sum(sum_energy_model), 2)}")
        axes[1,0].set_title(r"Energy Fit: $\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[1,0].set_xlabel("Reconstructed Energy (MeV)", fontsize = font_s)
        axes[1,0].set_ylabel("Counts", fontsize = font_s)
        axes[1,0].set_ylim((0, 65))
        axes[1,0].set_xlim((2.5, 5.0))
        
        # plot the dataset and error bars
        axes[1,0].errorbar(mids_energy, dataset_energy, yerr = np.sqrt(dataset_energy), color = "black", marker = "^", capsize = 2, linestyle = "", label = "Data")
        axes[1,0].legend(frameon = False, fontsize = 11)
        
        color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        for i in range(len(normalisations)):

            # plot each normalisation individually as a bar graph
            col = next(color_cycle)
            if col == '#d62728':
                col = next(color_cycle)
            axes[1,1].step(multisite_bins, multisite_result[i, :].tolist() + [0], linewidth = 2, where = 'post', color=col, label = f"{labels[i]}| Num: {round(np.sum(multisite_result[i,:]), 2)}")
        
        # plot the sum of the model
        sum_multisite_model = np.sum(multisite_result, axis = 0)
        axes[1,1].step(multisite_bins, sum_multisite_model.tolist() + [0], where = 'post', linewidth = 2, color = "red",  label = f"Total Model: {round(np.sum(sum_multisite_model), 2)}")
        axes[1,1].set_title(r"Multisite Fit: $\Delta log(\mathcal{L})$", fontsize = font_s)
        axes[1,1].set_xlabel("Multisite Discriminant", fontsize = font_s)
        axes[1,1].set_ylabel("Counts", fontsize = font_s)
        axes[1,1].set_xlim((-1.355, -1.335))
        axes[1,1].set_ylim((0, 65))
        
        # plot the dataset and error bars
        axes[1,1].errorbar(mids_multi, dataset_multisite, yerr = np.sqrt(dataset_multisite), color = "black", marker = "^", capsize = 2, linestyle = "", label = "Data")
        axes[1,1].legend(frameon = False, fontsize = 11)

        fig.tight_layout()
        plt.savefig("../plots/asimov_study/real_mc/advanced/fitted_background_model.png")
        plt.close()




expected_signal = 101.2
expected_backg  = [495.3, 0.93, 65.3, 38.5]
# expected_backg  = [495.3, None, None, None]
Ahab().run_analysis(expected_signal, expected_backg, False)