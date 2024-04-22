import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.optimize
import time
"""
Script is used to evaluate the bias and pull distributions of the fitted B8
signal. Bias is the percentage change between the fit result and the expected
value for a dataset, and the pull distribution allows us to judge whether the
fit result uncertainty is a good metric; if gaussian with unit std, we are good.

The bias and pull are calculated from performing the fit on many fake datasets.
Fake datasets are generated by sampling the binned PDFs for each isotope in the
background model (Bi212, Bi214, Tl210, Tl208, B8) N times, where N is drawn 
from a poisson distribution with mean = expectation for each isotope. 
"""
np.random.seed(42)

def calculate_uncertainty(profile_ll):
    """
    Function returns the 1 sigma confidence intervals on the minimum of profile LL
    fit.
    """
    
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

    return [combined_min, combined_upper, combined_lower], [multi_min, multi_upper, multi_lower], [energy_min, energy_upper, energy_lower]

def rescale_ll(ll):
        
    min_idx   = np.argmin(ll)
    diff_zero = 0 - np.ravel(ll)[min_idx]
    ll        = ll + diff_zero

    return ll 

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

def profile_likelihood_scan(fixed_backgrounds, initial_guess, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs):
    """
    Function performs a profile likelihood scan, calling the loglikelihood functions
    defined above via the scipy minimizer.
    """

    # the fixed number of B8 solar neutrinos assumed for each minimisation
    

    # arrays store the loglikelihood value for each assumed B8 number
    loglikelihood_array = np.zeros((3, len(signal_hypothesis)))
    normalisations      = np.zeros((3, len(signal_hypothesis), len(fixed_backgrounds) + 2))
    errors              = np.zeros((3, len(signal_hypothesis))) # single error for the Tl208 norm

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

        errors[0, counter] = np.sqrt(np.diag((combined_result.hess_inv.todense())))
        errors[1, counter] = np.sqrt(np.diag((multisite_result.hess_inv.todense())))
        errors[2, counter] = np.sqrt(np.diag((energy_result.hess_inv.todense())))

        counter += 1

    # return the best normalisations
    return loglikelihood_array, normalisations, errors

def generate_dataset():
    """
    Function creates a fluctuated Asimov dataset using inverse CDF sampling.
    
    OUTPUT: fake dataset as np.array of ints (binned) in multisite and
            energy.
    """

    # apply poisson fluctuation to the background model expectation
    fluctuated_number   = np.random.poisson(model_expectations)

    # keep track of which bins in the PDF are chosen for each isotope
    dataset_energy    = np.zeros((len(model_expectations), len(energy_bins) - 1)) 
    dataset_multisite = np.zeros((len(model_expectations), len(multisite_bins) - 1))

    # calculate the CDF over the columns --> keep the rows
    cdf_array_energy    = np.cumsum(pdf_array_energy, axis = 1)
    cdf_array_multisite = np.cumsum(pdf_array_multisite, axis = 1)
    
    # generate random numbers uniformly between 0 --> 1
    for isotope in range(len(fluctuated_number)):

        # generate number of uniform samples according to fluctuated number
        random_samples     = np.random.rand(fluctuated_number[isotope])

        # work out what bin idx each one falls into from cdf for each 
        bin_num_energy     = np.searchsorted(cdf_array_energy[isotope, :], random_samples)
        bin_num_multisite  = np.searchsorted(cdf_array_multisite[isotope, :], random_samples)

        # add 1 to every index specified in place
        np.add.at(dataset_energy[isotope, :], bin_num_energy, 1)      
        np.add.at(dataset_multisite[isotope, :], bin_num_multisite, 1)

    # create an output bar plot to check shit looks right
    # plt.figure()
    # for i in range(5):
    #     plt.step(energy_bins, dataset_energy[i, :].tolist() + [0], where = 'post', linewidth = 2, label = f"{isotope_names[i]}: {np.sum(dataset_energy[i, :])} Counts")
    #     # plt.step(multisite_bins, dataset_multisite[i, :].tolist() + [0], where = 'post', linewidth = 2, label = f"{isotope_names[i]}: {np.sum(dataset_multisite[i, :])} Counts")
    # plt.legend()
    # plt.savefig("../plots/asimov_study/real_mc/advanced/fake_dataset_energy.png")
    # plt.close()
        
    return dataset_energy, dataset_multisite

# define the background model expectations for each isotope
model_expectations  = np.array([101.2, 495.3, 0.93, 65.3, 38.5])
isotope_names       = ["B8", "Tl208", "Tl210", "Bi212", "Bi214"]
num_datasets        = 10000 # how many fake datasets to create? 
signal_hypothesis   = np.arange(0, 201, 1)

bias_bins           = np.arange(-1, +1, 0.05)
pull_bins           = np.arange(-10, +10, 0.25)
biases              = np.zeros((3, num_datasets))
pulls               = np.zeros((3, num_datasets))

# load the binned pdfs in energy and multisite
energy_bins         = np.arange(2.5, 5.05, 0.05)         # need to use the same binning as the pdfs were saved with
energy_mids         = energy_bins[:-1] + np.diff(energy_bins)[0] / 2
multisite_bins      = np.arange(-1.375, -1.325, 0.0005)
multisite_mids      = multisite_bins[:-1] + np.diff(multisite_bins)[0] / 2
pdf_array_energy    = np.load("./energy_pdf_array.npy")
pdf_array_multisite = np.load("./multisite_pdf_array.npy")

# run analysis for each fake dataset
start = time.time()
# for idataset in range(num_datasets):
#     if idataset % 100 == 0:
#         print(idataset, f"{time.time() - start} s.")
    
#     # generate fake dataset
#     dataset_energy, dataset_multisite = generate_dataset()

#     # perform profile log-likelihood analysis
#     profile_ll, norms, errors = profile_likelihood_scan(model_expectations[2:].tolist(), model_expectations[1], dataset_energy, dataset_multisite, pdf_array_energy, pdf_array_multisite)
    
#     # scale LL so the minimum of each curve is at 0
#     profile_ll[0, :] = rescale_ll(profile_ll[0, :])
#     profile_ll[1, :] = rescale_ll(profile_ll[1, :])
#     profile_ll[2, :] = rescale_ll(profile_ll[2, :])

#     # find 1 sigma frequentist confidence interval and fitted minimums 
#     combined_error, multisite_error, energy_error = calculate_uncertainty(profile_ll)

#     # find the norms corresponding to the minimum LL value for each LL function
#     minimum_idx = np.argmin(profile_ll, axis = 1)

#     # find the normalisations fom energy, multisite and combined LL giving that minimum
#     norms_combined = norms[0, minimum_idx[0], :]
#     norms_multi    = norms[1, minimum_idx[1], :]
#     norms_energy   = norms[2, minimum_idx[2], :]

#     # get the errors on the fitted norm of the Tl208
#     error_combined = errors[0, minimum_idx[0]]
#     error_multi    = errors[1, minimum_idx[1]]
#     error_energy   = errors[2, minimum_idx[2]]

#     # create the model out of the scaled PDFs
#     combined_energy_result    = norms_combined[:, None] * pdf_array_energy
#     combined_multisite_result = norms_combined[:, None] * pdf_array_multisite
#     energy_result             = norms_energy[:, None]   * pdf_array_energy
#     multisite_result          = norms_multi[:, None]    * pdf_array_multisite

#     if idataset < 10:
#         # create a plot for a few of these
#         fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (18, 8))
#         axes[0].set_title("Profile Log-Likelihood Scan")
#         axes[0].axvline(np.sum(dataset_energy[0, :]), color = "red", linestyle = "dotted", label = f"True Signal Counts: {np.sum(dataset_energy[0, :])}")
#         axes[0].plot(signal_hypothesis, profile_ll[2,:], color = "orange", label = "Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$")
#         axes[0].plot(signal_hypothesis, profile_ll[1, :], color = "green", label = "Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$")
#         axes[0].plot(signal_hypothesis, profile_ll[0, :], color = "black", label = "Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$")
#         axes[0].axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$ frequentist")
#         axes[0].legend(frameon = False, loc = "upper left")
#         axes[0].set_xlabel("Signal Hypothesis")
#         axes[0].set_ylabel(r"$-2log(\mathcal{L})$")
#         axes[0].set_ylim((0, 3))

#         for i in range(5):

#             # plot the 'true' underlying sampled data for each isotope
#             axes[1].step(energy_bins, dataset_energy[i, :].tolist() + [0], where = 'post', linestyle = "--", alpha = 0.5, linewidth = 2, label = f"{isotope_names[i]}: {np.sum(dataset_energy[i, :])} Counts")
#             axes[2].step(multisite_bins, dataset_multisite[i, :].tolist() + [0], where = 'post', linestyle = "--", alpha = 0.5, linewidth = 2, label = f"{isotope_names[i]}: {np.sum(dataset_multisite[i, :])} Counts")
        
#         # plot the total datasets (what the detector would see in real life)
#         axes[1].errorbar(energy_mids, np.sum(dataset_energy, axis = 0), marker = "^", linestyle = "", yerr = np.sqrt(np.sum(dataset_energy, axis = 0)), capsize = 2, color = "black", label = f"Total Dataset: {np.sum(dataset_energy)}")
#         axes[2].errorbar(multisite_mids, np.sum(dataset_multisite, axis = 0), marker = "^", linestyle = "", yerr = np.sqrt(np.sum(dataset_multisite, axis = 0)), capsize = 2, color = "black", label = f"Total Dataset: {np.sum(dataset_multisite)}")
        
#         # plot the fitted model for energy, multisite and combined together
#         axes[1].step(energy_bins, np.sum(energy_result, axis = 0).tolist() + [0], color = "red", where = 'post', linewidth = 2, label = f"Energy Fit: {np.sum(energy_result):.3g} Counts")
#         axes[1].step(energy_bins, np.sum(combined_energy_result, axis = 0).tolist() + [0], color = "#8B0000", where = 'post', linewidth = 2, label = f"Combined Fit: {np.sum(combined_energy_result):.3g} Counts")
#         axes[2].step(multisite_bins, np.sum(multisite_result, axis = 0).tolist() + [0], color = "red", where = 'post', linewidth = 2, label = f"Multisite Fit: {np.sum(multisite_result):.3g} Counts")
#         axes[2].step(multisite_bins, np.sum(combined_multisite_result, axis = 0).tolist() + [0], color = "#8B0000", where = 'post', linewidth = 2, label = f"Combined Fit: {np.sum(combined_multisite_result):.3g} Counts")


#         axes[1].set_xlabel("Reconstructed Energy (MeV)")
#         axes[2].set_xlabel("Multisite Discriminant")
#         axes[1].set_ylabel("Counts")
#         axes[2].set_ylabel("Counts")
#         axes[1].set_title("Fake Dataset: Energy")
#         axes[2].set_title("Fake Dataset: Multisite")
#         axes[2].set_xlim((-1.355, -1.334))
#         axes[1].legend(frameon = False, loc = "upper left")
#         axes[2].legend(frameon = False, loc = "upper left")
#         axes[1].set_ylim((0, 75))
#         axes[2].set_ylim((0, 75))
#         fig.tight_layout()
#         plt.savefig(f"../plots/asimov_study/real_mc/advanced/fluctuated_datasets/scan_{idataset}.png")
#         plt.close()

#     # calculate the bias and pull using each fit method
#     biases[0, idataset] = (norms_energy[0]   - model_expectations[0]) / model_expectations[0]
#     biases[1, idataset] = (norms_multi[0]    - model_expectations[0]) / model_expectations[0]
#     biases[2, idataset] = (norms_combined[0] - model_expectations[0]) / model_expectations[0]

#     if norms_energy[0] <= model_expectations[0]:
#         pulls[0, idataset] = ( model_expectations[0] - norms_energy[0] ) / energy_error[1]
#     else:
#         pulls[0, idataset] = ( norms_energy[0] - model_expectations[0] ) / -energy_error[2]

#     if norms_multi[0] <= model_expectations[0]:
#         pulls[1, idataset] = ( model_expectations[0] - norms_multi[0] ) / multisite_error[1]
#     else:
#         pulls[1, idataset] = ( norms_multi[0] - model_expectations[0] ) / -multisite_error[2]
        
#     if norms_combined[0] <= model_expectations[0]:
#         pulls[2, idataset] = ( model_expectations[0] - norms_combined[0] ) / combined_error[1]
#     else:
#         pulls[2, idataset] = ( norms_combined[0] - model_expectations[0] ) / -combined_error[2]

# np.save("./bias.npy", biases)
# np.save("./pull.npy", pulls)
biases = np.load("./bias.npy")
print(len(biases[np.isinf(biases)]))
biases[np.isinf(biases)] = -999

pulls  = np.load("./pull.npy")
print(len(pulls[np.isinf(pulls)]))
pulls[np.isinf(pulls)] = -999

# create a mask to calculate the mean of the distributions without including these bad values
bias_mask  = (biases != -999)
pulls_mask = (pulls != -999)
print(bias_mask.shape)
bias_means = np.ma.masked_array(biases, ~bias_mask).mean(axis = 1)
print(bias_means)
bias_std   = np.ma.masked_array(biases, ~bias_mask).std(axis = 1)
pulls_means = np.ma.masked_array(pulls, ~pulls_mask).mean(axis = 1)
pulls_std   = np.ma.masked_array(pulls, ~pulls_mask).std(axis = 1)

# make a plot of bias
fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 6))

axes[0].hist(biases[0,:], bins = bias_bins, histtype = "step", color = "orange")
axes[0].hist(biases[1,:], bins = bias_bins, histtype = "step", color = "green")
axes[0].hist(biases[2,:], bins = bias_bins, histtype = "step", color = "black")
axes[0].plot([], [], color = "orange", label = rf"Energy | $\mu$: {bias_means[0]:.2g}, $\sigma$: {bias_std[0]:.2g}")
axes[0].plot([], [], color = "green", label = rf"Multisite | $\mu$: {bias_means[1]:.2g}, $\sigma$: {bias_std[1]:.2g}")
axes[0].plot([], [], color = "black", label = rf"Combined | $\mu$: {bias_means[2]:.2g}, $\sigma$: {bias_std[2]:.2g}")

axes[0].set_xlabel("Bias")
axes[0].set_ylabel("Counts")
axes[0].legend(frameon = False, loc = "upper left")

axes[1].hist(pulls[0,:], bins = pull_bins, histtype = "step", color = "orange")
axes[1].hist(pulls[1,:], bins = pull_bins, histtype = "step", color = "green")
axes[1].hist(pulls[2,:], bins = pull_bins, histtype = "step", color = "black")
axes[1].plot([], [], color = "orange", label = rf"Energy | $\mu$: {pulls_means[0]:.2g}, $\sigma$: {pulls_std[0]:.2g}")
axes[1].plot([], [], color = "green", label = rf"Multisite | $\mu$: {pulls_means[1]:.2g}, $\sigma$: {pulls_std[1]:.2g}")
axes[1].plot([], [], color = "black", label = rf"Combined | $\mu$: {pulls_means[2]:.2g}, $\sigma$: {pulls_std[2]:.2g}")

axes[1].set_xlabel("Pull")
axes[1].set_ylabel("Counts")
axes[1].legend(frameon = False, loc = "upper left")
fig.tight_layout()
plt.savefig("../plots/asimov_study/real_mc/advanced/bias_fake_data.png")
plt.close()