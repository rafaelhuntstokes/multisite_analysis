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
        scaled_log    = np.sum(events_per_bin * log)

        return norm_sum - scaled_log
        # return -scaled_log
    
    def obtain_pdf(self, location, fv_cut, multisite_bins, energy_bins, run_list):
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
                ntuple = file.Get("2p5_5p0")
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

    def run_analysis(self):
        """
        Function that actually calls all the other functions and gets
        the analysis done.
        """
        
        # analysis settings
        expected_B8    = 101.2 # both these figures come from 4.5 m FV background model calculations in 2.5 --> 5.0 MeV ROI
        expected_Tl208 = 495.3 
        fv_cut         = 4500
        energy_bins    = np.arange(2.5, 5.1, 0.1)
        multisite_bins = np.arange(-1.4, -1.2, 0.0005)

        """
        Extract information and create the binned PDFS for energy shape and multisite,
        for each isotope of interest.
        """
        pdf_runlist    = np.loadtxt("../runlists/test_runlist.txt", dtype = int)
        mc_path        = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
        sig_mc_path    = f"{mc_path}/full_analysis_B8_solar_nue"
        backg_mc_path  = f"{mc_path}/full_analysis_Tl208"
        
        # signal and background pdfs extracted here
        sig_energy_pdf, sig_multi_pdf = self.obtain_pdf(sig_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist)
        backg_energy_pdf, backg_multi_pdf = self.obtain_pdf(backg_mc_path, fv_cut, multisite_bins, energy_bins, pdf_runlist)
        
        # put PDFs into 2D array (N_pdfs, N_bins) --> make this nicer and less hard-cody for arbitrary normalisations
        energy_pdf_array    = np.zeros((2, len(energy_bins) - 1))
        multisite_pdf_array = np.zeros((2, len(multisite_bins) - 1))
        energy_pdf_array[0, :]    = sig_energy_pdf
        energy_pdf_array[1, :]    = backg_energy_pdf
        multisite_pdf_array[0, :] = sig_multi_pdf
        multisite_pdf_array[1, :] = backg_multi_pdf

        # create the Asimov dataset from the PDFs themselves (add together the signal and background PDFs + apply normalisation)
        asimov_energy    = expected_B8 * sig_energy_pdf + expected_Tl208 * backg_energy_pdf
        asimov_multisite = expected_B8 * sig_multi_pdf  + expected_Tl208 * backg_multi_pdf

        # perform 2D grid searches over each background and the signal
        total_events     = expected_B8 + expected_Tl208  # we observe exactly the expected number in the Asimov dataset
        sig_hypothesis   = np.arange(80, 120, 1) # grid search in steps of 1
        backg_hypothesis = np.arange(400, 600, 2) 
        ll_multisite, ll_energy, ll_full = self.grid_search(sig_hypothesis, backg_hypothesis, multisite_pdf_array, energy_pdf_array, asimov_multisite, asimov_energy)
        
        # rescale all the log-likelihood functions so minimum value is at zero --> for plotting & comparisons
        ll_full = self.rescale_ll(ll_full)
        ll_multisite = self.rescale_ll(ll_multisite)
        ll_energy = self.rescale_ll(ll_energy)
        
        # create plots
        self.create_plots(sig_hypothesis, backg_hypothesis, ll_full, ll_multisite, ll_energy)

Ahab().run_analysis()