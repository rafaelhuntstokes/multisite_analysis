import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.optimize
import ROOT

"""
Simple profile likelihood fitter using unconstrained Tl208 number, and fixed
normalisations for the backgrounds Bi214, Bi212 and Tl210.
"""

def evaluate_loglikelihood(variable_background, fixed_backgrounds, signal_number, dataset, pdfs):
    """
    Function evaluates the extended binned log-likelihood value for a given
    number of signal events (B8).
    """

    # create array of total normalisations
    normalisations   = np.array([signal_number] + [variable_background] + fixed_backgrounds)

    # we allow the total normalisation to float so see what it is
    norm_sum         = np.sum(normalisations) # common to both 

    # have to do this inside the function for scipy.optimize.minimize doesn't
    # respect the reshaping being done outside its own calls ... 
    normalisations   = normalisations[:, None]

    """
    For the full derivation of this FULLY OPTIMIZED method, see fancy journal dated 3rd April 2024.
    """

    # multiply the counts of each PDF by the appropriate normalisation factor
    pdf_norm         = normalisations * pdfs   # elementwise multiplication results in (N_pdfs, N_bins)

    # sum each column to get the inside of each log term contributing to outer sum over bins
    inner_log        = np.sum(pdf_norm, axis = 0) # remove the rows so axis = 0

    log              = np.log(inner_log)          # 1D array of where each element is log(Nj Pj) for every bin in PDF

    # multiply each element by the number of data events in each PDF bin and sum all the terms
    scaled_log       = np.sum(dataset * log)

    # convert to 2 * -log(L)
    loglikelihood    = 2 * (norm_sum - scaled_log)


    return loglikelihood

def evaluate_combined_loglikelihood(variable_background, fixed_backgrounds, signal_number, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs):
    """
    Same as above, but calculate the combined loglikelihood.
    """

    # create array of total normalisations
    normalisations   = np.array([signal_number] + [variable_background] + fixed_backgrounds)

    # we allow the total normalisation to float so see what it is
    norm_sum         = np.sum(normalisations) # common to both 

    # have to do this inside the function for scipy.optimize.minimize doesn't
    # respect the reshaping being done outside its own calls ... 
    normalisations   = normalisations[:, None]

    """
    For the full derivation of this FULLY OPTIMIZED method, see fancy journal dated 3rd April 2024.
    """

    # multiply the counts of each PDF by the appropriate normalisation factor
    energy_pdf_norm         = normalisations * energy_pdfs   # elementwise multiplication results in (N_pdfs, N_bins)
    multisite_pdf_norm      = normalisations * multisite_pdfs

    # sum each column to get the inside of each log term contributing to outer sum over bins
    energy_inner_log        = np.sum(energy_pdf_norm, axis = 0) # remove the rows so axis = 0
    multisite_inner_log     = np.sum(multisite_pdf_norm, axis = 0)

    energy_log              = np.log(energy_inner_log)          # 1D array of where each element is log(Nj Pj) for every bin in PDF
    multisite_log           = np.log(multisite_inner_log)

    # multiply each element by the number of data events in each PDF bin and sum all the terms
    energy_scaled_log       = np.sum(energy_dataset * energy_log)
    multisite_scaled_log    = np.sum(multisite_dataset * multisite_log)

    # convert to 2 * -log(L)
    energy_loglikelihood    = 2 * (norm_sum - energy_scaled_log)
    multisite_loglikelihood = 2 * (norm_sum - multisite_scaled_log)

    combined_loglikelihood  = energy_loglikelihood + multisite_loglikelihood

    return combined_loglikelihood

def profile_likelihood_scan(fixed_backgrounds, initial_guess, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs):
    """
    Function performs a profile likelihood scan, calling the loglikelihood functions
    defined above via the scipy minimizer.
    """

    # the fixed number of B8 solar neutrinos assumed for each minimisation
    

    # arrays store the loglikelihood value for each assumed B8 number
    loglikelihood_array = np.zeros((3, len(signal_hypothesis)))
    normalisations      = np.zeros((3, len(signal_hypothesis), len(fixed_backgrounds) + 2))
    
    # loop over each assumed signal number and perform minimisation
    counter = 0
    for isig in signal_hypothesis:

        energy_result    = scipy.optimize.minimize(evaluate_loglikelihood, x0 = initial_guess, method = "L-BFGS-B", tol = 1e-4, bounds = [(0, np.inf)], args = (fixed_backgrounds, isig, energy_dataset, energy_pdfs))
        multisite_result = scipy.optimize.minimize(evaluate_loglikelihood, x0 = initial_guess, method = "L-BFGS-B", tol = 1e-4, bounds = [(0, np.inf)], args = (fixed_backgrounds, isig, multisite_dataset, multisite_pdfs))

        # use minimizer to get combined LL result
        combined_result  = scipy.optimize.minimize(evaluate_combined_loglikelihood, x0 = initial_guess, method = "L-BFGS-B", tol = 1e-4, bounds = [(0, np.inf)], args = (fixed_backgrounds, isig, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs))

        loglikelihood_array[0, counter] = combined_result.fun
        loglikelihood_array[1, counter] = multisite_result.fun
        loglikelihood_array[2, counter] = energy_result.fun

        normalisations[0, counter, :] = [isig] + combined_result.x.tolist() + fixed_backgrounds
        normalisations[1, counter, :] = [isig] + multisite_result.x.tolist() + fixed_backgrounds
        normalisations[2, counter, :] = [isig] + energy_result.x.tolist() + fixed_backgrounds

        counter += 1

    # return the best normalisations
    return loglikelihood_array, normalisations

def obtain_pdf(location, fv_cut, multisite_bins, energy_bins, run_list, energy_range, plot_name):
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

        return energy_counts, multisite_counts

def rescale_ll(ll):
        
    min_idx   = np.argmin(ll)
    diff_zero = 0 - np.ravel(ll)[min_idx]
    ll        = ll + diff_zero

    return ll 

def obtain_dataset():

    run_list = np.loadtxt("../../data_studies/runlists/quiet_period.txt", dtype = int)

    energy_vals = []
    multi_vals  = []
    count = 0
    for irun in run_list:

        file = ROOT.TFile.Open(f"../../data_studies/extracted_data/full_analysis/processed_dataset/{irun}.root")
        ntuple = file.Get("2p5_5p0")

        for ientry in ntuple:
            energy_vals.append(ientry.energy)
            multi_vals.append(ientry.dlogL)

        count += 1
        print(count)

    return energy_vals, multi_vals

# first obtain the PDFs for each background
# binning for the energy and multisite discriminant PDFs
energy_bins       = np.arange(2.5, 5.1, 0.05)
multisite_bins    = np.arange(-1.375, -1.325, 0.0005)
backg_names       = ["Tl208", "Tl210", "Bi212", "Bi214"] 
labels            = ["B8"] + backg_names
mids_energy       = energy_bins[:-1] + np.diff(energy_bins)[0] / 2
mids_multi        = multisite_bins[:-1] + np.diff(multisite_bins)[0] / 2
signal_hypothesis = np.arange(0, 201, 1)
analyse_real_data = True
expected_signal   = 101.2
expected_backg    = [495.3, 0.93, 65.3, 38.5]
normalisations    = np.array([expected_signal] + expected_backg)

"""
Extract information and create the binned PDFS for energy shape and multisite,
for each isotope of interest.
"""
# we will create an Asimov dataset using the multisite and energy PDFs of every included normalisation
# pdf_runlist    = np.loadtxt("../runlists/full_test.txt", dtype = int)
# mc_path        = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
# sig_mc_path    = f"{mc_path}/full_analysis_B8_solar_nue" # signal path is always the same

# loop over the included normalisations

# create array containing the signal and background pdfs for multisite and energy
# multisite_pdf_array = np.zeros((5, multisite_bins.size - 1))  # (number backgrounds + 1, number of bins)
# energy_pdf_array    = np.zeros((5, energy_bins.size - 1))    # (number backgrounds + 1, number of bins)

# obtain the signal PDF
# print("Obtaining signal PDF.")
# sig_energy_pdf, sig_multi_pdf = obtain_pdf(sig_mc_path, 4500.0, multisite_bins, energy_bins, pdf_runlist, "2p5_5p0", "B8_nue")

# add the signal pdfs as the first row in the pdfs arrays
# multisite_pdf_array[0, :] = sig_multi_pdf
# energy_pdf_array[0, :]    = sig_energy_pdf

# obtain background pdfs and add to subsequent rows
# idx_counter = 1
# for iname in range(len(backg_names)):
    
    # print("Obtaining background pdf: ", backg_names[iname])
    # we have an expected rate and want to include this normalisation in the Asimov dataset
    # backg_mc_path  = f"{mc_path}/full_analysis_{backg_names[iname]}"

    # obtain the normalised PDFs for this background and add to the respective arrays
    # backg_energy_pdf, backg_multi_pdf   = obtain_pdf(backg_mc_path, 4500.0, multisite_bins, energy_bins, pdf_runlist, "2p5_5p0", backg_names[iname])
    # multisite_pdf_array[idx_counter, :] = backg_multi_pdf
    # energy_pdf_array[idx_counter, :]    = backg_energy_pdf

    # incremement idx counter ... not the same idx as the number of backg names!
    # idx_counter += 1

# np.save("./energy_pdf_array.npy", energy_pdf_array)
# np.save("./multisite_pdf_array.npy", multisite_pdf_array)
energy_pdf_array    = np.load("./energy_pdf_array.npy")
multisite_pdf_array = np.load("./multisite_pdf_array.npy")

# create the dataset --> either by building an Asimov dataset, or using real data
# if analyse_real_data == False:
#     normalisations_stretched = normalisations.copy() 
#     normalisations_stretched = normalisations_stretched[:, None]

#     print("Creating Asimov dataset from normalisation-scaled PDFs.")
#     dataset_energy    = normalisations_stretched * energy_pdf_array       # multiply every bin by corresponding normalisation
#     dataset_energy    = np.sum(dataset_energy, axis = 0)                  # remove the rows so axis = 0

#     dataset_multisite = normalisations_stretched * multisite_pdf_array    # multiply every bin by corresponding normalisation
#     dataset_multisite = np.sum(dataset_multisite, axis = 0)               # remove the rows so axis = 0
# else:
#     # use real data
#     energy_vals, multisite_vals = obtain_dataset()
#     dataset_energy, _    = np.histogram(energy_vals, bins = energy_bins)
#     dataset_multisite, _ = np.histogram(multisite_vals, bins = multisite_bins)

# np.save("./energy_dataset_asimov.npy", dataset_energy)
# np.save("./multisite_dataset_asimov.npy", dataset_multisite)
# np.save("./energy_dataset_real.npy", dataset_energy)
# np.save("./multisite_dataset_real.npy", dataset_multisite)

if analyse_real_data == False:
    dataset_energy    = np.load("./energy_dataset_asimov.npy")
    dataset_multisite = np.load("./multisite_dataset_asimov.npy")
    plot_name         = "asimov"
else:
    dataset_energy    = np.load("./energy_dataset_real.npy")
    dataset_multisite = np.load("./multisite_dataset_real.npy")
    plot_name         = "real"

# now, whatever data we have, perform the profile likelihood scan
profile_ll, norms = profile_likelihood_scan(expected_backg[1:], expected_backg[0], dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array)

profile_ll[0, :] = rescale_ll(profile_ll[0, :])
profile_ll[1, :] = rescale_ll(profile_ll[1, :])
profile_ll[2, :] = rescale_ll(profile_ll[2, :])

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
combined_min   = signal_hypothesis[np.argmin(profile_ll[0,:])]
combined_upper = abs(combined_min - signal_hypothesis[second_closest[0]]) 
combined_lower = abs(combined_min - signal_hypothesis[closest_to_interval[0]])
multi_min      = signal_hypothesis[np.argmin(profile_ll[1,:])]
multi_upper    = abs(multi_min - signal_hypothesis[second_closest[1]]) 
multi_lower    = abs(multi_min - signal_hypothesis[closest_to_interval[1]])
energy_min     = signal_hypothesis[np.argmin(profile_ll[2,:])]
energy_upper   = abs(energy_min - signal_hypothesis[second_closest[2]]) 
energy_lower   = abs(energy_min - signal_hypothesis[closest_to_interval[2]])

# create a plot of these
plt.title("Profile Log-Likelihood Scan")
plt.plot(signal_hypothesis, profile_ll[2,:], color = "orange", label = "Energy | " + rf"${round(energy_min, 1)}^{{+{round(energy_upper, 1)}}}_{{-{round(energy_lower, 1)}}}$")
plt.plot(signal_hypothesis, profile_ll[1, :], color = "green", label = "Multisite | " + rf"${round(multi_min, 1)}^{{+{round(multi_upper, 1)}}}_{{-{round(multi_lower, 1)}}}$")
plt.plot(signal_hypothesis, profile_ll[0, :], color = "black", label = "Combined | " + rf"${round(combined_min, 1)}^{{+{round(combined_upper, 1)}}}_{{-{round(combined_lower, 1)}}}$")
plt.axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$ frequentist")
plt.legend()
plt.xlabel("Signal Hypothesis")
plt.ylabel(r"$-2log(\mathcal{L})$")
plt.ylim((0, 3))
plt.savefig(f"../plots/asimov_study/real_mc/advanced/unconstrained_profileLL_{plot_name}.png")
plt.close()

# create a plot of the 3 fitted background + signal models vs data

# find the norms corresponding to the minimum LL value for each LL function
minimum_idx = np.argmin(profile_ll, axis = 1)

# find the normalisations fom energy, multisite and combined LL giving that minimum
norms_combined = norms[0, minimum_idx[0], :]
norms_multi    = norms[1, minimum_idx[1], :]
norms_energy   = norms[2, minimum_idx[2], :]

# create the model out of the scaled PDFs
combined_energy_result    = norms_combined[:, None] * energy_pdf_array
combined_multisite_result = norms_combined[:, None] * multisite_pdf_array
energy_result             = norms_energy[:, None]   * energy_pdf_array
multisite_result          = norms_multi[:, None]    * multisite_pdf_array

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 15))
font_s    = 14
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
plt.savefig(f"../plots/asimov_study/real_mc/advanced/fitted_background_model_{plot_name}.png")
plt.close()



