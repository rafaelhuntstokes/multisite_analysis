import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json

class LikelihoodExtraction():
    """
    Class contains a bucket of functions to run a simple log-likelihood fit of the 
    number of B8 solar neutrinos in my model.
    
    One function runs the fit without likelihood, using essentially a chi2 scan.
    
    The second runs the fit with multisite discrimination added. 
    
    The analysis is run in different energy bins between 2.5 --> 5.0 MeV.
    
    Plots of the likelihood minimisation in each energy bin, plus the final result 
    (model vs MC with best fit B8 number across all energy bins) are returned.
    """

    def obtain_data_counts_in_energy_bins(self, fv_cut, roi_start, roi_end, bin_width, data_path, run_list):
        """
        Function obtains the energy spectrum of the data for a given FV, and 
        bins it according to the the input bin_width, roi_start, roi_end 
        parameters.

        Also obtains the log-likelihood ratio discriminant of each event in the
        dataset.
        
        Returns the number of events in each energy bin in the data.
        """

        # check the FV desired exists
        if fv_cut != 3000 and fv_cut != 3500 and fv_cut != 4000 and fv_cut != 4500 and fv_cut != 5000 and fv_cut != 5500 and fv_cut != 6000:
            print("FV cut must match the previously applied analysis cuts: 3000, 3500, 4000, 4500, 5000, 5500, 6000")
            return 0

        if run_list != "directionality":
            print("Run list does not exist. Options are: directionality")
            return 0
        
        # obtain the dataset of interest based on run list
        data = np.load(f"{data_path}/energy_dist_{run_list}_{fv_cut}_counts.npy")
        print("Loaded the dataset.")

        # bin these data according to the energy binning schema
        binning        = np.arange(roi_start, roi_end + bin_width, bin_width)
        print(binning)
        counts_data, _ = np.histogram(data, bins = binning)
        print("Binned the dataset.")

        # obtain the dlogL for each event in the energy bins
        ##########################
        #### TO BE IMPLEMENTED ### 
        ##########################

        # return the number of data counts in each energy bin
        return counts_data
    
    def obtain_model_counts_in_energy_bins(self, fv_cut, roi_start, roi_end, bin_width, model_path):
        """
        Function creates a background model for a given FV from MC, binned 
        according to the input roi_start, roi_end, bin_width parameters.

        Background model expected counts are obtained from background model 
        calculations done on Google Drive (TM). The information is stored in a
        JSON dictionary file for each FV.

        Returns the number of background counts for each energy bin, excl. B8. 
        """

        def apply_scaling(energy_spectrum, expected_number, binning):
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

            return scaled_counts, mid_points[:-1] 

        # obtain all the energy distributions from MC for each component of the model
        energy_bi214     = np.load(model_path + "/BiPo214_mc_information/energy.npy")
        energy_po214     = np.load(model_path + "/Po214_mc_information/energy.npy")
        energy_bi212     = np.load(model_path + "/BiPo212_mc_information/energy.npy")
        energy_po212     = np.load(model_path + "/Po212_mc_information/energy.npy")
        energy_alphaN_p  = np.load(model_path + "/AlphaN_LAB_13C_mc_information/energy.npy") # prompt spectrum + pileup
        energy_alphaN_d  = np.load(model_path + "/gammas_2p2MeV_mc_information/energy.npy")  # delayed 2.2 MeV gamma spectrum
        energy_tl208     = np.load(model_path + "/Tl208_mc_information/energy.npy")
        energy_tl210     = np.load(model_path + "/Tl210_mc_information/energy.npy")
        energy_B8_nue    = np.load(model_path + "/B8_solar_nue_mc_information/energy.npy")
        energy_B8_numu   = np.load(model_path + "/B8_solar_numu_mc_information/energy.npy")
        energy_pa234m    = np.load(model_path + "/Pa234m_mc_information/energy.npy")

        radius_bi214     = np.load(model_path + "/BiPo214_mc_information/posR.npy")
        radius_po214     = np.load(model_path + "/Po214_mc_information/posR.npy")
        radius_bi212     = np.load(model_path + "/BiPo212_mc_information/posR.npy")
        radius_po212     = np.load(model_path + "/Po212_mc_information/posR.npy")
        radius_alphaN_p  = np.load(model_path + "/AlphaN_LAB_13C_mc_information/posR.npy") 
        radius_alphaN_d  = np.load(model_path + "/gammas_2p2MeV_mc_information/posR.npy")
        radius_tl208     = np.load(model_path + "/Tl208_mc_information/posR.npy")
        radius_tl210     = np.load(model_path + "/Tl210_mc_information/posR.npy")
        radius_B8_nue    = np.load(model_path + "/B8_solar_nue_mc_information/posR.npy")
        radius_B8_numu   = np.load(model_path + "/B8_solar_numu_mc_information/posR.npy")
        radius_pa234m    = np.load(model_path + "/Pa234m_mc_information/posR.npy")
        print("Loaded MC information.")

        # isolate MC events which fall within the FV cut
        idx_bi214     = np.where(radius_bi214 <fv_cut)
        idx_po214     = np.where(radius_po214 < fv_cut)
        idx_bi212     = np.where(radius_bi212 < fv_cut)
        idx_po212     = np.where(radius_po212 < fv_cut)
        idx_alphaN_p  = np.where(radius_alphaN_p < fv_cut)
        idx_alphaN_d  = np.where(radius_alphaN_d < fv_cut)
        idx_tl208     = np.where(radius_tl208 < fv_cut)
        idx_tl210     = np.where(radius_tl210 < fv_cut)
        idx_B8_nue    = np.where(radius_B8_nue < fv_cut)
        idx_B8_numu   = np.where(radius_B8_numu < fv_cut)
        idx_pa234m    = np.where(radius_pa234m < fv_cut)

        # apply idx
        energy_bi214     = energy_bi214[idx_bi214]
        energy_po214     = energy_po214[idx_po214]
        energy_bi212     = energy_bi212[idx_bi212]
        energy_po212     = energy_po212[idx_po212]
        energy_alphaN_p  = energy_alphaN_p[idx_alphaN_p]
        energy_alphaN_d  = energy_alphaN_d[idx_alphaN_d]
        energy_tl208     = energy_tl208[idx_tl208]
        energy_tl210     = energy_tl210[idx_tl210]
        energy_B8_nue    = energy_B8_nue[idx_B8_nue]
        energy_B8_numu   = energy_B8_numu[idx_B8_numu]
        energy_pa234m    = energy_pa234m[idx_pa234m]

        # obtain the number of expected events for each event type
        # load JSON file containing dictionary of isotopes and their calculated expected rates
        expected_backgrounds_rates = open("expected_rates_database.json")
        background_rates_db        = json.load(expected_backgrounds_rates)
        name_map  = {3000.0: "3m", 3500.0: "3p5m", 4000.0: "4m", 4500.0: "4p5m", 5000.0: "5m", 5500.0: "5p5m", 6000.0: "6m"}
        name_tree = name_map[float(fv_cut)]
        # find the expected number of events for each background from the database --> assuming will always want the coincidence
        # cut values and not the full distribution
        expected_bipo214   = background_rates_db["post_cuts"][name_tree]["BiPo214"] 
        expected_bipo212   = background_rates_db["post_cuts"][name_tree]["BiPo212"]
        expected_tl208     = background_rates_db["post_cuts"][name_tree]["Tl208"]
        expected_tl210     = background_rates_db["post_cuts"][name_tree]["Tl210"]
        expected_B8_nue    = background_rates_db["post_cuts"][name_tree]["B8_nue"]
        expected_B8_numu   = background_rates_db["post_cuts"][name_tree]["B8_numu"]
        expected_alphaN_p  = background_rates_db["post_cuts"][name_tree]["alphaN_prompt"]
        expected_alphaN_d  = background_rates_db["post_cuts"][name_tree]["alphaN_delayed"]
        expected_pa234m    = background_rates_db["post_cuts"][name_tree]["Pa234m"]

        # scale the MC to the expected number of events
        ####### careful! ######################################################
        ## might be broken if the bins here don't line up with the ROI bins ... 
        #######################################################################
        binning            = np.arange(2.0, 6 + bin_width, bin_width) # expected rates calculated for 2 --> 6 MeV range
        scaled_bi214, _    = apply_scaling(energy_bi214, expected_bipo214, binning)
        scaled_bi212, _    = apply_scaling(energy_bi212, expected_bipo212, binning)
        scaled_tl208, _    = apply_scaling(energy_tl208, expected_tl208, binning)
        scaled_tl210, _    = apply_scaling(energy_tl210, expected_tl210, binning)
        scaled_alphaN_p, _ = apply_scaling(energy_alphaN_p, expected_alphaN_p, binning)
        scaled_alphaN_d, _ = apply_scaling(energy_alphaN_d, expected_alphaN_d, binning)
        scaled_B8_nue, _   = apply_scaling(energy_B8_nue, expected_B8_nue, binning)
        scaled_B8_numu, _  = apply_scaling(energy_B8_numu, expected_B8_numu, binning)
        scaled_pa234m, _   = apply_scaling(energy_pa234m, expected_pa234m, binning)
        print("Applied scaling to MC.")

        # work out the number of counts in each bin within the ROI
        bin_idx_roi_start = np.digitize(roi_start, bins = binning) - 1 # -1 since digitize says bin 1 is idx 1 and we need to count from idx 0
        bin_idx_roi_end   = np.digitize(roi_end, bins = binning)   - 2 # since start counting from 0 (-1) and up to 5 (which is upper bin edge and excluded by this function) (-1)
        print(scaled_tl208, bin_idx_roi_start, bin_idx_roi_end, binning, binning[bin_idx_roi_start], binning[bin_idx_roi_end])

        # find the total counts in the model (excluding B8!) within each of the ROI bins
        binned_roi_model_contribution = scaled_bi214[bin_idx_roi_start:bin_idx_roi_end + 1]    + scaled_bi212[bin_idx_roi_start:bin_idx_roi_end + 1]    \
                                      + scaled_tl208[bin_idx_roi_start:bin_idx_roi_end + 1]    + scaled_tl210[bin_idx_roi_start:bin_idx_roi_end + 1]   \
                                      + scaled_alphaN_p[bin_idx_roi_start:bin_idx_roi_end + 1] + scaled_alphaN_d[bin_idx_roi_start:bin_idx_roi_end + 1] \
                                      + scaled_pa234m[bin_idx_roi_start:bin_idx_roi_end + 1]

        # create an ASIMOV dataset, with a known contribution of B8 in it to check my fitting framework works
        num_b8_asimov     = [10, 20, 25, 30, 5]#scaled_B8_nue[bin_idx_roi_start:bin_idx_roi_end + 1] + scaled_B8_numu[bin_idx_roi_start:bin_idx_roi_end + 1]
        asmimov_roi_model = binned_roi_model_contribution + num_b8_asimov 

        return binned_roi_model_contribution, asmimov_roi_model, num_b8_asimov

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

    def perform_loglikelihood_multisite_minimisation(self, number_data, number_background, loglikelihood_data):
        """
        For a given energy bin, the number of B8 events is grid searched over. In 
        combination with the dlog(L) values and PDFs of the dataset, the likelihood
        for each number of B8 events is returned.
        """

        pass
    
    def create_output_plots(self, chi2_dists, values_b8, roi_bins):
        """
        Function creates the output plots: -2log(L) vs B8 number in each energy bin
        and the full reconstructed model vs data with best fit B8 events.
        """

        for i in range(len(chi2_dists)):
            # create output plot for each energy bin
            plt.figure()
            plt.plot(values_b8, chi2_dists[i], color = "black")
            plt.title(f"{roi_bins[i]} --> {roi_bins[i+1]} MeV")
            plt.xlabel("Number B8 Events")
            plt.ylabel(r"$\chi ^2$")
            plt.savefig(f"../plots/likelihood_extraction/simple_chi2_{roi_bins[i]}_{roi_bins[i+1]}MeV.png")
            plt.close()

    def run_analysis(self):
        """
        Master function that specifies the energy binning to use and calls the other functions.
        """


        binned_data_counts = self.obtain_data_counts_in_energy_bins(4500, 2.5, 5.0, 0.5, "/home/hunt-stokes/multisite_analysis", "directionality")
        # if binned_data_counts == False:
        #     print("Failed to extract data counts. Check validity of FV and run list.")
        #     return 0

        binned_model_counts, asimov_model_counts, num_b8_asimov = self.obtain_model_counts_in_energy_bins(4500, 2.5, 5.0, 0.5, "/data/snoplus3/hunt-stokes/clean_multisite")

        print("Expected data counts per bin:  ", binned_data_counts)
        print("Expected model counts per bin: ", binned_model_counts)

        values_b8 = np.arange(0, 101, 1)
        roi_bins  = np.arange(2.5, 5.5, 0.5)
        # maximise the likelihood in each energy bin
        chi2_dists_data   = [] # array storing the chi2 vs B8 number for every energy bin [[bin1 vals], [bin2 vals], ... ]
        chi2_dists_asimov = [] # array storing the chi2 vs B8 number for each energy bin in the asmimov dataset
        for i_energy_bin in range(len(binned_model_counts)):
            chi2 = self.perform_loglikelihood_minimisation(binned_data_counts[i_energy_bin], binned_model_counts[i_energy_bin], values_b8)
            chi2_dists_data.append(chi2)

            chi2 = self.perform_loglikelihood_minimisation(asimov_model_counts[i_energy_bin], binned_model_counts[i_energy_bin], values_b8)
            chi2_dists_asimov.append(chi2)

        # create output plots
        fig, axes = plt.subplots(nrows = 1, ncols = 5, figsize = (20, 6))
        for i in range(len(chi2_dists_data)):
            # create output plot for each energy bin
            axes[i].plot(values_b8, chi2_dists_data[i], color = "black")
            axes[i].set_title(f"{roi_bins[i]} --> {roi_bins[i+1]} MeV")
            axes[i].set_xlabel("Number B8 Events")
            axes[i].set_ylabel(r"$\chi ^2$")
        plt.tight_layout()
        plt.savefig(f"../plots/likelihood_extraction/simple_chi2.png")
        plt.close()
        
        fig, axes = plt.subplots(nrows = 1, ncols = 5, figsize = (20, 6))
        for i in range(len(chi2_dists_asimov)):
             # create output plot for each energy bin
            axes[i].plot(values_b8, chi2_dists_asimov[i], color = "black")
            axes[i].set_title(f"{roi_bins[i]} --> {roi_bins[i+1]} MeV")
            axes[i].set_xlabel("Number B8 Events")
            axes[i].set_ylabel(r"$\chi ^2$")

            # find most likely value of B8 in this bin and the true asimov value of B8
            most_likely_b8 = values_b8[np.argmin(chi2_dists_asimov[i])]
            asimov_value   = num_b8_asimov[i]

            axes[i].axvline(most_likely_b8, color = "red", linestyle = "dashed", label = f"Fitted value: {round(most_likely_b8, 3)}")
            axes[i].axvline(asimov_value, color = "black", label = f"Asimov Value: {round(asimov_value,3)}")
            axes[i].legend()
        plt.tight_layout()
        plt.savefig(f"../plots/likelihood_extraction/asimov_chi2.png")
        plt.close()

X = LikelihoodExtraction()
X.run_analysis()

