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
        
        ll        = np.array(ll)
        min_idx   = np.argmin(ll)
        diff_zero = 0 - ll[min_idx]
        ll        = ll + diff_zero

        return ll 

    def minimisation(self):
        pass

    def create_plots(self, signal_hypothesis, ll_full, ll_multi, ll_energy):
        """
        Function creates all the output plots. 
        ---> to be completed once minimisation code is working.
        """

        plt.figure()
        plt.plot(signal_hypothesis, ll_full, color = "black", label = f"Full | Min: {signal_hypothesis[np.argmin(ll_full)]}")
        plt.plot(signal_hypothesis, ll_multi, color = "green", label = f"Multisite | Min: {signal_hypothesis[np.argmin(ll_multi)]}")
        plt.plot(signal_hypothesis, ll_energy, color = "orange", label = f"Energy | Min: {signal_hypothesis[np.argmin(ll_energy)]}")
        plt.legend()
        plt.xlabel("Signal Hypothesis")
        plt.ylabel(r"$-2log(\mathcal{L})$")
        plt.xlim((80, 130))
        plt.ylim((0, 3))
        plt.savefig("../plots/asimov_study/real_mc/advanced/test.png")
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

        # perform gridsearch over Nsig, constrained such that Nsig + Nbackg = Ntotal
        total_events = expected_B8 + expected_Tl208      # we observe exactly the expected number in the Asimov dataset
        hypothesis   = np.arange(0, total_events, 1) # grid search in steps of 1
        
        ll_multisite = [] # track the ll for each BELL function for each norm_sig hypothesis
        ll_energy    = []
        ll_full      = []
        for norm_sig in hypothesis:
            
            # norm of backgrounds is constrained
            norm_backg     = total_events - norm_sig
            normalisations = np.array([norm_sig, norm_backg])

            # evaluate 2 x log-likelihood from energy and multisite PDF shapes
            multi  = 2 * self.evaluate_loglikelihood(multisite_pdf_array, normalisations, 0, 0, asimov_multisite)
            energy = 2 * self.evaluate_loglikelihood(energy_pdf_array, normalisations, 0, 0, asimov_energy) 
            full   = multi + energy

            # save result for this signal hypothesis
            ll_multisite.append(multi)
            ll_energy.append(energy)
            ll_full.append(full)

        # rescale all the log-likelihood functions so minimum value is at zero --> for plotting & comparisons
        ll_full = self.rescale_ll(ll_full)
        ll_multisite = self.rescale_ll(ll_multisite)
        ll_energy = self.rescale_ll(ll_energy)
        
        # create plots
        self.create_plots(hypothesis, ll_full, ll_multisite, ll_energy)

Ahab().run_analysis()