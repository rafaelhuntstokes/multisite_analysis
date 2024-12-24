import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.optimize
import ROOT
from flux_conversion import *
from prettytable import PrettyTable
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.font_manager as fm

path2 = '/data/snoplus3/hunt-stokes/nuPhysPoster/scripts/Times_New_Roman_Normal.ttf'
prop_font = fm.FontProperties(fname=path2, size = 12)

plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams.update({'axes.unicode_minus' : False})

# set position of axis labels
plt.rcParams["xaxis.labellocation"] = 'right'
plt.rcParams["yaxis.labellocation"] = 'top'

# set global parameters with rcParams -- copy whole block below before plotting code to set plotting template globally

# set height and width of big markings on axis x
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.6
# set height and width of small markings on axis x
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1.6
# set height and width of big markings on axis y
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 1.6
# set height and width of small markings on axis y
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 1.6
# set thickness of axes
plt.rcParams['axes.linewidth'] = 1.6
# set plot background color
plt.rcParams['figure.facecolor'] = 'white'
# set plot aspect ratio -- change according to needs
plt.rcParams['figure.figsize'] = (8.5, 6.5)
# set padding (between ticks and axis label)
plt.rcParams['xtick.major.pad'] = '12' ## change me if the axis labels overlap! ## 
plt.rcParams['ytick.major.pad'] = '12'
# set padding (between plot and title)
plt.rcParams['axes.titlepad'] = 12
# set markings on axis to show on the inside of the plot, can change if needed
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# set ticks on both sides of x and y axes
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

histogram_style = {
    'histtype': 'step', 
    'color': 'blue',
    'alpha': 0.7,
    'linewidth': 2
}

scatter_style = {
    'marker': 's',
    'color': 'black',
    's': 25
}

errorbar_style = {
    'linestyle': 'None',
    'color': 'black',
    'capsize': 1.5
}

line_plot_style = {
    'linestyle': '-',
    'color': 'black',
    'linewidth': 2.5
}
"""
Simple profile likelihood fitter using unconstrained Tl208 number, and fixed
normalisations for the backgrounds Bi214, Bi212 and Tl210.
"""
def evaluate_loglikelihood_vectorised_grid_search(normalisations, dataset_energy, dataset_multisite, pdfs_energy, pdfs_multisite):
    """
    This function is made to speed up the grid search method used instead of the
    optimiser, so the loops are only computed over the signal norms.

    It is expected the normalisations array is (N_tl208_norms, N_norms), with
    a row for each Tl208 signal normalisation for a given B8 number.

    The PDFs are of shape (N_norms, N_bins), as normal.
    """

    # get the sum of norms for each row (i.e. remove the cols)
    norm_sum             = np.sum(normalisations, axis = 1)

    # see log book entry 29th April 2024
    pdf_norm_energy      = normalisations[:, :, None] * pdfs_energy[None, :] # (Results in a 'row' for each Tl208 norm used)
    pdf_norm_multisite   = normalisations[:, :, None] * pdfs_multisite[None, :]

    # calculate the sum of the log terms for each Tl208 norms (over each row in each Tl208 block)
    inner_log_energy     = np.sum(pdf_norm_energy, axis = 1) # result is (N_tl208, Nbins)
    inner_log_multisite  = np.sum(pdf_norm_multisite, axis = 1)

    # element-wise logarithm - easy
    log_energy           = np.log(inner_log_energy)
    log_multisite        = np.log(inner_log_multisite) 

    # multiply by the counts in the dataset and sum each Tl208 norm
    scaled_log_energy    = np.sum(dataset_energy * log_energy, axis = 1) # result is (N_tl208, ) 1D
    scaled_log_multisite = np.sum(dataset_multisite * log_multisite, axis = 1)
    
    # convert to a -2log(L) value for each Tl208 norm
    loglikelihood_energy    = 2 * (norm_sum - scaled_log_energy)    # elementwise - easy
    loglikelihood_multisite = 2 * (norm_sum - scaled_log_multisite)
    loglikelihood_combined  = loglikelihood_energy + loglikelihood_multisite

    return loglikelihood_combined, loglikelihood_multisite, loglikelihood_energy

def vectorised_profile_likelihood(fixed_backgrounds):
    """
    Calls the vectorised version of the evaluate_loglikelihood function for bias / pull studies.
    """

    # arrays store the loglikelihood value for each assumed B8 number
    loglikelihood_array = np.zeros((3, len(signal_hypothesis)))
    
    # fit result norms
    fitted_norms = np.zeros((3, len(signal_hypothesis), len(fixed_backgrounds) + 2))

    # create an array containing the normalisations for each Tl208
    normalisations = np.zeros((len(background_hypothesis), len(fixed_backgrounds) + 2))

    # fill in all the Tl208 norms
    normalisations[:, 1] = background_hypothesis

    # fill in the fixed backgrounds
    normalisations[:, 2:] = fixed_backgrounds

    # single error for the Tl208 norm
    errors              = np.zeros((3, len(signal_hypothesis))) 
    
    # loop over each assumed signal number and perform minimisation
    counter = 0
    for sig in signal_hypothesis:

        # fill in the signal values
        normalisations[:, 0] = sig
                
        # evaluate the loglikelihood with each method for each Tl208 norm
        combined_LL, multisite_LL, energy_LL = evaluate_loglikelihood_vectorised_grid_search(normalisations, dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array)

        # find the minimum value of the LL for each method and store
        loglikelihood_array[0, counter] = np.min(combined_LL)  # combined
        loglikelihood_array[1, counter] = np.min(multisite_LL) # multisite
        loglikelihood_array[2, counter] = np.min(energy_LL)    # energy
        
        # find the corresponding normalisations to this minimum
        fitted_norms[0, counter, :] = [sig] + [background_hypothesis[np.argmin(combined_LL)]]  + fixed_backgrounds
        fitted_norms[1, counter, :] = [sig] + [background_hypothesis[np.argmin(multisite_LL)]] + fixed_backgrounds
        fitted_norms[2, counter, :] = [sig] + [background_hypothesis[np.argmin(energy_LL)]]    + fixed_backgrounds

        # error on the fitted norms --> setting to zero for now as 'not important'
        errors[0, counter] = 0
        errors[1, counter] = 0
        errors[2, counter] = 0

        counter += 1

    # return the best normalisations
    return loglikelihood_array, normalisations, errors

def evaluate_loglikelihood(variable_background, fixed_backgrounds, signal_number, dataset, pdfs):
    """
    Function evaluates the extended binned log-likelihood value for a given
    number of signal events (B8).

    variable_background: normalisation of Tl208, float
    fixed_backgrounds:   array of floats [Tl210, Bi212, Bi214]
    signal_number:       norm of B8, float
    dataset:             binned counts in either energy or multisite of the dataset, array of floats length nbins
    pdfs:                array of binned PDFs in either multisite or energy (same binning as dataset)
                         contains the pdf for EACH normalisation -- > dimensions (n_pdfs, nbins)
    """

    # create array of total normalisations --> [B8 norm, Tl208 Norm, Tl210 norm, Bi212 Norm, Bi214 norm]
    normalisations   = np.array([signal_number] + [variable_background] + fixed_backgrounds)
    
    # normalisations   = np.array([signal_number] + [variable_background] + fixed_backgrounds)
    # we allow the total normalisation to float
    norm_sum         = np.sum(normalisations)   # sum the normalisations as the first term in eqn

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
    # expected_backg_conv   = np.array([66.3] + expected_backg)
    # constraint       = (normalisations - expected_backg_conv)**2 / (2 * np.array(expected_uncertainty)**2)
    # convert to 2 * -log(L)
    # loglikelihood    = 2 * (norm_sum - scaled_log + ( (normalisations[1] - expected_backg[0])**2/(2*tl208_uncert**2) ) )
    loglikelihood    = 2 * (norm_sum - scaled_log)       # element wise operation
    # print(f"B8: {signal_number}, Tl208: {variable_background}, LL: {loglikelihood}")
    # print(f"\nnorm sum: {norm_sum} | log: {log.shape} | normalisations: {normalisations} | | dataset shape: {dataset.shape} | pdf_shape: {pdfs.shape} | pdf_norms_shape: {pdf_norm.shape}\n")
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
    # print(normalisations)
    energy_pdf_norm         = normalisations * energy_pdfs   # elementwise multiplication results in (N_pdfs, N_bins)
    multisite_pdf_norm      = normalisations * multisite_pdfs

    # sum each column to get the inside of each log term contributing to outer sum over bins
    energy_inner_log        = np.sum(energy_pdf_norm, axis = 0) # remove the rows so axis = 0
    multisite_inner_log     = np.sum(multisite_pdf_norm, axis = 0)
    # print(energy_inner_log)
    energy_log              = np.log(energy_inner_log)          # 1D array of where each element is log(Nj Pj) for every bin in PDF
    multisite_log           = np.log(multisite_inner_log)

    # multiply each element by the number of data events in each PDF bin and sum all the terms
    energy_scaled_log       = np.sum(energy_dataset * energy_log)
    multisite_scaled_log    = np.sum(multisite_dataset * multisite_log)

    # convert to 2 * -log(L)
    # energy_loglikelihood    = 2 * (norm_sum - energy_scaled_log + ( (normalisations[1] - expected_backg[0])**2/(2*tl208_uncert**2) ))
    # multisite_loglikelihood = 2 * (norm_sum - multisite_scaled_log + ( (normalisations[1] - expected_backg[0])**2/(2*tl208_uncert**2) ))
    energy_loglikelihood    = 2 * (norm_sum - energy_scaled_log)
    multisite_loglikelihood = 2 * (norm_sum - multisite_scaled_log)

    combined_loglikelihood  = energy_loglikelihood + multisite_loglikelihood

    return combined_loglikelihood

def profile_likelihood_scan(fixed_backgrounds, initial_guess, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs):
    """
    Function performs a profile likelihood scan, calling the loglikelihood functions
    defined above via the scipy minimizer.
    """    

    # arrays store the loglikelihood value for each assumed B8 number
    loglikelihood_array = np.zeros((3, len(signal_hypothesis)))
    normalisations      = np.zeros((3, len(signal_hypothesis), len(fixed_backgrounds) + 2))
    errors              = np.zeros((3, len(signal_hypothesis))) # single error for the Tl208 norm
    
    # loop over each assumed signal number and perform minimisation
    counter = 0
    for sig in signal_hypothesis:
        
        # save the LL profiles for this B8 norm for each method
        intermediate_values = np.zeros((3, len(background_hypothesis)))
        
        # grid search over the Tl208 normalisations
        for iback in range(len(background_hypothesis)):
            intermediate_values[0, iback] = evaluate_combined_loglikelihood(background_hypothesis[iback], fixed_backgrounds, sig, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs)
            intermediate_values[1, iback] = evaluate_loglikelihood(background_hypothesis[iback], fixed_backgrounds, sig, multisite_dataset, multisite_pdfs)
            intermediate_values[2, iback] = evaluate_loglikelihood(background_hypothesis[iback], fixed_backgrounds, sig, energy_dataset, energy_pdfs)

        # find the minimum value of the LL for each method and store
        loglikelihood_array[0, counter] = np.min(intermediate_values[0, :]) # combined
        loglikelihood_array[1, counter] = np.min(intermediate_values[1, :]) # multisite
        loglikelihood_array[2, counter] = np.min(intermediate_values[2, :]) # energy
        
        # find the corresponding normalisations to this minimum
        normalisations[0, counter, :] = [sig] + [background_hypothesis[np.argmin(intermediate_values[0, :])]] + fixed_backgrounds
        normalisations[1, counter, :] = [sig] + [background_hypothesis[np.argmin(intermediate_values[1, :])]] + fixed_backgrounds
        normalisations[2, counter, :] = [sig] + [background_hypothesis[np.argmin(intermediate_values[2, :])]] + fixed_backgrounds

        # error on the fitted norms --> setting to zero for now as 'not important'
        errors[0, counter] = 0
        errors[1, counter] = 0
        errors[2, counter] = 0

        counter += 1

    # return the best normalisations
    return loglikelihood_array, normalisations, errors

def likelihood_gridsearch(fixed_backgrounds, initial_guess, energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs):
    """
    Function performs a 2D grid search to exhaustively explore the LL space. The
    minimum Tl208 norm for each B8 norm is returned and plotted.
    
    This is to compare with that found using the profile likelihood method. 
    
    The overall minimum and confidence interval (in 2D) is plotted.
    """

    l_space_energy    = np.zeros((len(signal_hypothesis), len(background_hypothesis)))
    l_space_multisite = np.zeros((len(signal_hypothesis), len(background_hypothesis)))
    l_space_combined  = np.zeros((len(signal_hypothesis), len(background_hypothesis)))

    for isig in range(len(signal_hypothesis)):
        for iback in range(len(background_hypothesis)):

            # evaluate log-likelihood with each of the 3 methods
            l_space_energy[isig, iback]    = evaluate_loglikelihood(background_hypothesis[iback], fixed_backgrounds, signal_hypothesis[isig], energy_dataset, energy_pdfs)
            l_space_multisite[isig, iback] = evaluate_loglikelihood(background_hypothesis[iback], fixed_backgrounds, signal_hypothesis[isig], multisite_dataset, multisite_pdfs)
            l_space_combined[isig, iback]  = evaluate_combined_loglikelihood(background_hypothesis[iback], fixed_backgrounds, signal_hypothesis[isig], energy_dataset, multisite_dataset, energy_pdfs, multisite_pdfs)

    # rescale the LL
    l_space_energy = rescale_ll(l_space_energy)
    l_space_multisite = rescale_ll(l_space_multisite)
    l_space_combined  = rescale_ll(l_space_combined)

    np.save("./l_space_energy.npy", l_space_energy)
    np.save("./l_space_multisite.npy", l_space_multisite)
    np.save("./l_space_combined.npy", l_space_combined)
    # find minimum of each
    idx_energy = np.argwhere(l_space_energy == np.min(l_space_energy))
    idx_multisite = np.argwhere(l_space_multisite == np.min(l_space_multisite))
    idx_combined = np.argwhere(l_space_combined == np.min(l_space_combined))
    print(idx_energy, idx_multisite, idx_combined)

    # find the minimum for the Tl208 profiles
    # find the minimum LL for each value of Tl208 norm (the cols)
    tl208_energy_profile = np.min(l_space_energy, axis = 0)
    tl208_multi_profile  = np.min(l_space_multisite, axis = 0)
    tl208_combi_profile  = np.min(l_space_combined, axis = 0)

    np.save(f"./tl208_profile_energy_THESIS.npy", tl208_energy_profile)
    np.save(f"./tl208_profile_multi_THESIS.npy", tl208_multi_profile)
    np.save(f"./tl208_profile_combi_THESIS.npy", tl208_combi_profile)

    fig, axes = plt.subplots(nrows = 3, ncols = 1, figsize = (8,11))
    # fig.set_size_inches(((8.27, 11.69)))
    # im = axes.imshow(np.log(l_space_energy), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # # im = axes[0].imshow(l_space_energy, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # axes.scatter(background_hypothesis[idx_energy[:, 1]], signal_hypothesis[idx_energy[:, 0]], color = "red", label = f"Minimum B8: {signal_hypothesis[idx_energy[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_energy[0,1]]:.3g}")
    # axes.scatter(495.3, 66.3, color = "blue", label = "Background Model Assumption")
    # divider = make_axes_locatable(axes)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(im, cax = cax)
    # axes.set_title("Energy", fontsize = 12)
    # axes.set_xlabel("Tl208 Normalisation", fontsize = 12)
    # axes.set_ylabel("B8 Normalisation", fontsize = 12)
    # axes.set_aspect("equal")
    # axes.set_ylim((800, 0))
    # axes.set_xlim((0,800))
    # axes.legend(fontsize = 12)

    im = axes[0].imshow(np.log(l_space_energy), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]], aspect = 'auto')
    # im = axes[0].imshow(l_space_energy, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    axes[0].scatter(background_hypothesis[idx_energy[:, 1]], signal_hypothesis[idx_energy[:, 0]], color = "red", label = f"Minimum B8: {signal_hypothesis[idx_energy[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_energy[0,1]]:.3g}")
    axes[0].scatter(expected_backg[0], expected_signal, alpha = 0.5, color = "blue", label = f"Expected (B8, Tl208): ({expected_signal:.3g}, {expected_backg[0]:.3g})")
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax = cax)
    # axes[0].set_title("Energy")
    axes[0].set_xlabel("Tl208 Normalisation", fontproperties = prop_font)
    axes[0].set_ylabel("B8 Normalisation", fontproperties = prop_font)
    axes[0].set_aspect("equal", adjustable = 'box')

    maxLL = np.max(np.log(l_space_energy))
    dchi2 = 2 * (maxLL - np.log(l_space_energy))
    # axes[0].contour(dchi2, levels = [2.3], colors = 'white', extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # axes[0].set_ylim((800, 0))
    # axes[0].set_xlim((0,800))
    # axes[0].set_ylim((background_hypothesis[0], background_hypothesis[-1]))
    # axes[0].set_xlim((signal_hypothesis[0], signal_hypothesis[-1]))
    for label in axes[0].get_xticklabels():
            label.set_fontproperties(prop_font)

    for label in axes[0].get_yticklabels():
        label.set_fontproperties(prop_font)
    axes[0].legend(frameon = True, prop=fm.FontProperties(fname=path2, size = 14))
    im = axes[1].imshow(np.log(l_space_multisite), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]], aspect = 'auto')
    # im = axes[1].imshow(l_space_multisite, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    axes[1].scatter(background_hypothesis[idx_multisite[:, 1]], signal_hypothesis[idx_multisite[:, 0]], color = "red", label = f"Minimum B8: {signal_hypothesis[idx_multisite[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_multisite[0,1]]:.3g}")
    axes[1].scatter(expected_backg[0], expected_signal, color = "blue", alpha = 0.5, label = f"Expected (B8, Tl208): ({expected_signal:.3g}, {expected_backg[0]:.3g})")
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax = cax)
    # axes[1].set_title("Multisite")
    axes[1].set_aspect("equal", adjustable = 'box')
    axes[1].set_xlabel("Tl208 Normalisation", fontproperties = prop_font)
    axes[1].set_ylabel("B8 Normalisation", fontproperties = prop_font)
    axes[1].set_ylim(axes[0].get_ylim())
    axes[1].legend(frameon = True, prop=fm.FontProperties(fname=path2, size = 14))
    for label in axes[1].get_xticklabels():
            label.set_fontproperties(prop_font)

    for label in axes[1].get_yticklabels():
        label.set_fontproperties(prop_font)
    maxLL = np.max(np.log(l_space_multisite))
    dchi2 = 2 * (maxLL - np.log(l_space_multisite))
    # axes[1].contour(dchi2, levels = [2.3], colors = 'white', extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # im = axes.imshow(np.log(l_space_multisite), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # # im = axes[1].imshow(l_space_multisite, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # axes.scatter(background_hypothesis[idx_multisite[:, 1]], signal_hypothesis[idx_multisite[:, 0]], color = "red", label = f"Minimum B8: {signal_hypothesis[idx_multisite[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_multisite[0,1]]:.3g}")
    # axes.scatter(495.3, 66.3, color = "blue", label = "Background Model Assumption")
    # divider = make_axes_locatable(axes)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(im, cax = cax)
    # axes.set_title("Multisite")
    # axes.set_xlabel("Tl208 Normalisation")
    # axes.set_ylabel("B8 Normalisation")
    # axes.set_ylim((800, 0))
    # axes.set_xlim((0,800))
    # axes.legend(fontsize = 12)
    # axes.legend()

    im = axes[2].imshow(np.log(l_space_combined), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]], aspect = 'auto')
    # im = axes[2].imshow(l_space_combined, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    axes[2].scatter(background_hypothesis[idx_combined[:, 1]], signal_hypothesis[idx_combined[:, 0]], color = "red", label = f"Minimum B8: {signal_hypothesis[idx_combined[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_combined[0, 1]]:.3g}")
    axes[2].scatter(expected_backg[0], expected_signal, alpha = 0.5, color = "blue", label = f"Expected (B8, Tl208): ({expected_signal:.3g}, {expected_backg[0]:.3g})")
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax = cax)
    # axes[2].set_title("Combined")
    axes[2].set_xlabel("Tl208 Normalisation", fontproperties = prop_font)
    axes[2].set_ylabel("B8 Normalisation", fontproperties = prop_font)
    axes[2].set_aspect("equal", adjustable = 'box')
    axes[2].legend(frameon = True, prop=fm.FontProperties(fname=path2, size = 14))

    for label in axes[2].get_xticklabels():
            label.set_fontproperties(prop_font)

    for label in axes[2].get_yticklabels():
        label.set_fontproperties(prop_font)
    maxLL = np.max(np.log(l_space_combined))
    dchi2 = 2 * (maxLL - np.log(l_space_combined))
    # axes[2].contour(dchi2, levels = [2.3], colors = 'white', extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    fig.tight_layout()
    # plt.savefig(f"../plots/asimov_study/real_mc/advanced/diff_bins/heatmap_{energy_string}_{num_bins}.png")
    plt.savefig(f"../plots/asimov_study/real_mc/advanced/heatmap_{energy_string}_{plot_name}.pdf")
    plt.close()

    return l_space_energy

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
                # print(f"{location}/{irun}.root")
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
                itr = ientry.itr
                # posFOM = ientry.posFOM / ientry.posFOM_hits
                if r > fv_cut:
                    continue
                # if itr < 0.22 or itr > 0.3:
                #     continue
                # if posFOM < 13.7:
                #     continue
                # print("dlogL: ", ientry.dlogL)
                # event passed FV so add info to PDFs
                energies.append(ientry.energy)
                multisites.append(ientry.dlogL)
        print(f"Found {len(missing_runs)} missing runs:\n{missing_runs}")
        # print("min - max: ", np.min(multisites), np.max(multisites))
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

def obtain_dataset(energy_range):

    run_list = np.loadtxt("../../data_studies/runlists/quiet_period.txt", dtype = int)

    energy_vals = []
    multi_vals  = []
    count = 0
    for irun in run_list:
        # if energy_string == "2p5_3p75" or energy_string == "3p75_5p0" or energy_string == "3p5_5p0":
        try:
            file = ROOT.TFile.Open(f"../../data_studies/extracted_data/full_analysis3/reprocessed_dataset_7.0.15/{irun}.root")
        except:
            continue
        # else:
            # file = ROOT.TFile.Open(f"../../data_studies/extracted_data/full_analysis/processed_dataset/{irun}.root")
        ntuple = file.Get(f"{energy_range}")

        for ientry in ntuple:
            x = ientry.x
            y = ientry.y
            z = ientry.z
            r = np.sqrt(x**2 + y**2 + z**2)
            if r > FV_CUT:
                continue
            # itr = ientry.itr
            # if itr < 0.22 or itr > 0.3:
            #     continue
            # posFOM = ientry.posFOM / ientry.posFOM_hits
            # if posFOM < 13.7:
            #     continue
            energy_vals.append(ientry.energy)
            multi_vals.append(ientry.dlogL)
            # print(ientry.dlogL)

        count += 1
        # print(count)

    return energy_vals, multi_vals

def create_pdfs_and_datasets(analyse_real_data, energy_string):
    """
    Extract information and create the binned PDFS for energy shape and multisite,
    for each isotope of interest.
    """
    # we will create an Asimov dataset using the multisite and energy PDFs of every included normalisation
    pdf_runlist    = np.loadtxt("../runlists/full_test.txt", dtype = int)
    mc_path        = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
    if energy_string == "2p5_3p75" or energy_string == "3p75_5p0" or energy_string == "3p5_5p0":
        sig_mc_path = f"{mc_path}/full_analysis3_B8_solar_nue"
        backg_path  = "full_analysis3"
    else:
        sig_mc_path    = f"{mc_path}/full_analysis_B8_solar_nue"#normed_pdfs_4p5m_B8_solar_nue" # signal path is always the same
        backg_path = "full_analysis"#normed_pdfs_4p5m"
    # loop over the included normalisations

    # create array containing the signal and background pdfs for multisite and energy
    multisite_pdf_array = np.zeros((5, multisite_bins.size - 1))  # (number backgrounds + 1, number of bins)
    energy_pdf_array    = np.zeros((5, energy_bins.size - 1))    # (number backgrounds + 1, number of bins)

    # obtain the signal PDF
    print("Obtaining signal PDF.")
    sig_energy_pdf, sig_multi_pdf = obtain_pdf(sig_mc_path, FV_CUT, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue")
    print(sig_energy_pdf)
    # add the signal pdfs as the first row in the pdfs arrays
    multisite_pdf_array[0, :] = sig_multi_pdf
    energy_pdf_array[0, :]    = sig_energy_pdf

    # obtain background pdfs and add to subsequent rows
    idx_counter = 1
    for iname in range(len(backg_names)):
        
        print("Obtaining background pdf: ", backg_names[iname])

        # we have an expected rate and want to include this normalisation in the Asimov dataset
        backg_mc_path  = f"{mc_path}/{backg_path}_{backg_names[iname]}"
        print(backg_mc_path)
        # obtain the normalised PDFs for this background and add to the respective arrays
        backg_energy_pdf, backg_multi_pdf   = obtain_pdf(backg_mc_path, FV_CUT, multisite_bins, energy_bins, pdf_runlist, energy_string, backg_names[iname])
        multisite_pdf_array[idx_counter, :] = backg_multi_pdf
        energy_pdf_array[idx_counter, :]    = backg_energy_pdf
        print(multisite_pdf_array)

        # incremement idx counter ... not the same idx as the number of backg names!
        idx_counter += 1

    np.save(f"./energy_pdf_array_{energy_string}.npy", energy_pdf_array)
    np.save(f"./multisite_pdf_array_{energy_string}.npy", multisite_pdf_array)

    # create the dataset --> either by building an Asimov dataset, or using real data
    if analyse_real_data == False:
        normalisations_stretched = normalisations.copy() 
        normalisations_stretched = normalisations_stretched[:, None]

        print("Creating Asimov dataset from normalisation-scaled PDFs.")
        dataset_energy    = normalisations_stretched * energy_pdf_array       # multiply every bin by corresponding normalisation
        dataset_energy    = np.sum(dataset_energy, axis = 0)                  # remove the rows so axis = 0

        dataset_multisite = normalisations_stretched * multisite_pdf_array    # multiply every bin by corresponding normalisation
        dataset_multisite = np.sum(dataset_multisite, axis = 0)               # remove the rows so axis = 0
        np.save(f"./energy_dataset_asimov_{energy_string}.npy", dataset_energy)
        np.save(f"./multisite_dataset_asimov_{energy_string}.npy", dataset_multisite)
    else:
        # use real data
        energy_vals, multisite_vals = obtain_dataset(energy_string)
        dataset_energy, _    = np.histogram(energy_vals, bins = energy_bins)
        dataset_multisite, _ = np.histogram(multisite_vals, bins = multisite_bins)
        np.save(f"./energy_dataset_real_{energy_string}.npy", dataset_energy)
        np.save(f"./multisite_dataset_real_{energy_string}.npy", dataset_multisite)
    return dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array

def create_fit_model_subplot(ax, font_s, fit_result, data, fit_bins, data_mids, xlims, words, fit_error):
    """
    Create a subplot of the fitted model vs data. Only non-obvious input is words: 
    
    words = ["title", "xlabel"]
    """

    # create the 4 subplots showing relative contributions of each normalisation, the sum and the data
    color_cycle = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for i in range(len(normalisations)):

        # plot each normalisation individually as a bar graph
        col = next(color_cycle)

        if col == '#d62728':
            col = next(color_cycle)
        ax.step(fit_bins, fit_result[i, :].tolist() + [0], where = 'post', linewidth = 2, color=col, label = f"{labels[i]} Num: {round(np.sum(fit_result[i,:]), 2)}")

    # plot the sum of the model
    sum_fit_model = np.sum(fit_result, axis = 0)
    ax.step(fit_bins, sum_fit_model.tolist() + [0], where = 'post', color = "red", linewidth = 2, label = f"Total Model: {round(np.sum(sum_fit_model), 2)}")
    
    # find the chi2 / NDF
    # chi2 = (data - sum_fit_model)**2 / 2*np.sqrt(sum_fit_model)**2

    # remove zero bins
    # idx_include = np.where(sum_fit_model >= 1)
    # NDF  = len(data) - 1 # number of bins - number of pdfs (1)
    # print(NDF)
    # ax.plot([], [], linestyle = "", label = r"$\delta _{^{208}Tl} \pm$" + f" {fit_error}")
    ax.plot([], [], linestyle = "", label = f"Num in Dataset: {sum(data):.1f}")

    # ax.set_title(words[0], fontsize = font_s)
    ax.set_xlabel(words[1], fontsize = font_s, fontproperties = prop_font)
    ax.set_ylabel("Counts", fontsize = font_s, fontproperties = prop_font)
    # ax.set_ylim((0, 65))
    ax.set_xlim(xlims)

    # plot the dataset and error bars
    ax.errorbar(data_mids, data, yerr = np.sqrt(data), color = "black", marker = "^", capsize = 2, linestyle = "", label = "Data", markersize = 4)
    ax.legend(frameon = False, fontsize = font_s, prop=fm.FontProperties(fname=path2, size = 14))
    # ax.legend(frameon = False, fontsize = 11)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, labelsize=font_s)  # Hide x-axis ticks and labels
    
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop_font)
    # Plot residuals
    residuals   = data - sum_fit_model
    residual_ax = ax.inset_axes([0, -0.2, 1, 0.2])
    residual_ax.axhline(y = 0, color = "black", linewidth = 1, linestyle = "dotted")
    residual_ax.scatter(data_mids, residuals, marker='o', s = 5, color='black')
    residual_ax.set_xlabel(words[1], fontsize = font_s, fontproperties = prop_font)
    residual_ax.set_ylabel(r"Data - Model", fontsize = font_s, fontproperties = prop_font)
    residual_ax.set_xlim(xlims)
    residual_ax.set_ylim((-15, 15))
    residual_ax.tick_params(axis = "both", labelsize = font_s)
    yticks = residual_ax.get_yticks()
    yticks = np.arange(-15, 20, 5)
    yticks[-1] = 0  # Replace the last tick label with 0
    residual_ax.set_yticks(yticks)
    print(f"Total value of deltas: {np.sum(residuals)}")

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

def evaluate_bias():
    """
    Function uses poisson-fluctuated asimov dataset to evaluate the bias in the
    fitting framework itself. The bias is the difference between the fitted
    minimum of profile LL and the true asimov value.
    """

    # load the asimov dataset (assuming it exists already)
    asimov_dataset_energy    = np.load("./energy_dataset_asimov.npy")
    asimov_dataset_multisite = np.load("./multisite_dataset_asimov.npy")

    # run the analysis on N poisson fluctuated datasets
    bias_combined  = []
    bias_multisite = []
    bias_energy    = []
    pull_combined  = []
    pull_multisite = []
    pull_energy    = []
    start = time.time()
    for i in range(10):

        fluctuated_energy    = np.random.poisson(asimov_dataset_energy)
        fluctuated_multisite = np.random.poisson(asimov_dataset_multisite)

        old_profile_ll, old_norms, old_errors = profile_likelihood_scan(expected_backg[1:], expected_backg[0], fluctuated_energy, fluctuated_multisite, energy_pdf_array, multisite_pdf_array)

        profile_ll, norms, errors = vectorised_profile_likelihood(expected_backg[1:])
        print(i)
        # find the minimum of the ll functions and see how far from true value it is
        minimum_idx   = np.argmin(profile_ll, axis = 1)
        min_combined  = signal_hypothesis[minimum_idx[0]]
        min_multisite = signal_hypothesis[minimum_idx[1]]
        min_energy    = signal_hypothesis[minimum_idx[2]]

        bias_combined.append( (min_combined - expected_signal) / expected_signal )
        bias_multisite.append( (min_multisite - expected_signal) / expected_signal )
        bias_energy.append( (min_energy - expected_signal) / expected_signal )

        # evaluate the pull of the distribution
        combined_error, multisite_error, energy_error = calculate_uncertainty(profile_ll)

        if min_combined < expected_signal:
            if min_combined < 20:
                print(min_combined, combined_error[1])#
                plt.plot(signal_hypothesis, profile_ll[0, :])
                plt.savefig("../plots/asimov_study/real_mc/advanced/fked_result.png")
                plt.close()
            pull_combined.append(( min_combined - expected_signal ) / combined_error[1] )
        else:
            pull_combined.append(( min_combined - expected_signal ) / combined_error[2] )
        
        if min_multisite < expected_signal:
            # print(min_multisite)
            pull_multisite.append(( min_multisite - expected_signal ) / multisite_error[1] )
        else:
            pull_multisite.append(( min_multisite - expected_signal ) / multisite_error[2] )

        if min_energy < expected_signal:
            # print(min_energy)
            pull_energy.append(( min_energy - expected_signal ) / energy_error[1] )
        else:
            pull_energy.append(( min_energy - expected_signal ) / energy_error[2] )
    end = time.time()
    print(f"Took {end -start} s; {(end-start)/10} s per dataset.")
    np.save("./bias_combined.npy", bias_combined)
    np.save("./bias_multisite.npy", bias_multisite)
    np.save("./bias_energy.npy", bias_energy)
    np.save("./pull_combined.npy", pull_combined)
    np.save("./pull_multisite.npy", pull_multisite)
    np.save("./pull_energy.npy", pull_energy)
    bias_combined  = np.load("./bias_combined.npy") 
    bias_multisite = np.load("./bias_multisite.npy")
    bias_energy    = np.load("./bias_energy.npy")
    pull_combined  = np.load("./pull_combined.npy") 
    pull_multisite = np.load("./pull_multisite.npy")
    pull_energy    = np.load("./pull_energy.npy")

    # find the mean and std of each bias distribution
    mean_combined  = np.mean(bias_combined)
    std_combined   = np.std(bias_combined)
    mean_multisite = np.mean(bias_multisite)
    std_multisite  = np.std(bias_multisite)
    mean_energy    = np.mean(bias_energy)
    std_energy     = np.std(bias_energy)

    binning = np.arange(-100/ expected_signal, 100/ expected_signal, 2/ expected_signal)
    plt.figure(figsize = (12,8))
    plt.hist(bias_combined, bins = binning, color = "black", histtype = "step",  linewidth = 2)
    plt.hist(bias_multisite, bins = binning, color = "green", histtype = "step", linewidth = 2)
    plt.hist(bias_energy, bins = binning, color = "orange", histtype = "step", linewidth = 2)

    # for nice labels
    plt.plot([], [], color = "black", label = rf"Combined | $\mu$: {mean_combined:.3g}, $\sigma$: {std_combined:.3g}", linewidth = 2)
    plt.plot([], [], color = "green", label = rf"Multisite | $\mu$: {mean_multisite:.3g}, $\sigma$: {std_multisite:.3g}", linewidth = 2)
    plt.plot([], [], color = "orange", label = rf"Energy | $\mu$: {mean_energy:.3g}, $\sigma$: {std_energy:.3g}", linewidth = 2)

    plt.xlabel(r"$\frac{Fit - Expected}{Expected}$", fontsize = font_s)
    plt.xlim((-1, 1))
    plt.ylabel("Counts", fontsize = font_s)
    plt.legend(fontsize = font_s, frameon=False)
    plt.savefig("../plots/asimov_study/real_mc/advanced/bias.png")
    plt.close()

    plt.figure(figsize = (12, 8))
    binning = np.arange(-1, 1, 3/ expected_signal)
    # find the mean and std of each bias distribution
    mean_combined  = np.mean(pull_combined[(pull_combined >= -1 and pull_combined <= 1)])
    std_combined   = np.std(pull_combined[(pull_combined >= -1 and pull_combined <= 1)])
    mean_multisite = np.mean(pull_multisite[(pull_multisite >= -1 and pull_multisite <= 1)])
    std_multisite  = np.std(pull_multisite[(pull_multisite >= -1 and pull_multisite <= 1)])
    mean_energy    = np.mean(pull_energy[(pull_energy >= -1 and pull_energy <= 1)])
    std_energy     = np.std(pull_energy[(pull_energy >= -1 and pull_energy <= 1)])

    plt.plot([], [], color = "black", label = rf"Combined | $\mu$: {mean_combined:.3g}, $\sigma$: {std_combined:.3g}", linewidth = 2)
    plt.plot([], [], color = "green", label = rf"Multisite | $\mu$: {mean_multisite:.3g}, $\sigma$: {std_multisite:.3g}", linewidth = 2)
    plt.plot([], [], color = "orange", label = rf"Energy | $\mu$: {mean_energy:.3g}, $\sigma$: {std_energy:.3g}", linewidth = 2)

    plt.hist(pull_combined, bins = binning, color = "black", histtype = "step", linewidth = 2)
    plt.hist(pull_multisite, bins = binning, color = "green", histtype = "step", linewidth = 2)
    plt.hist(pull_energy, bins = binning, color = "orange",  histtype = "step", linewidth = 2)
    plt.legend(frameon = False, fontsize = font_s)
    plt.savefig("../plots/asimov_study/real_mc/advanced/pull.png")
    plt.close()

def convert_fitted_tl_to_thorium(fitted_tl208):
    """
    Function takes the fitted Tl208 normalisation and converts this into a Th223
    concentration in the scintillator.
    """
    
    efficiency_4p5m = 0.1509 # (2.5 --> 5.0 MeV + 4.5 m Efficiency) x fraction of events in full ROI present in whatever smaller ROI used
    mass_av         = 7.84e8                # g
    Na              = 6.02e23               # avogadro number per mole
    mr_th232        = 232                   # molar mass Th232, g per mole
    branching_ratio = 0.3594                # ratio of Bi212 --> Tl208
    livetime        = 145.7                 # days
    half_life       = 14e9 * 3.154e7        # Th232 t1/2, seconds 
    decay_constant  = np.log(2) / half_life # per second

    # convert the fitted number in 4.5 m --> total number in 6 m 
    rate_in_av = fitted_tl208 / efficiency_4p5m # number in livetime

    # convert the rate into correct unit
    rate_in_av = rate_in_av / ( livetime * 24 * 60 * 60 ) # rate per second
    print(rate_in_av)
    # convert into a concentration
    th232_concentration = ( rate_in_av * mr_th232 ) / ( mass_av * Na * branching_ratio * decay_constant)

    return th232_concentration

def energy_bin_plots():
    """
    Function called to plot the energy only fit result (fit, profile LL, heatmap).
    This is useful when studing the impact of different numbers of energy bins on
    the fit result.
    """

    fig, axes = plt.subplots(nrows = 1, ncols = 3)
    fig.set_size_inches((8.27, 11.69*0.3))

    # create the energy fit result subplot
    create_fit_model_subplot(axes[0], 5, energy_result, dataset_energy, energy_bins, mids_energy, (2.5, 5.0), [r"Energy Fit", "Reconstructed Energy (MeV)"], error_energy)

    # create the 2D heatmap
    idx_energy = np.argwhere(l_space_energy == np.min(l_space_energy))
    im = axes[1].imshow(np.log(l_space_energy), cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    # im = axes[0].imshow(l_space_energy, cmap = "BuGn", extent = [background_hypothesis[0], background_hypothesis[-1], signal_hypothesis[-1], signal_hypothesis[0]])
    axes[1].scatter(background_hypothesis[idx_energy[:, 1]], signal_hypothesis[idx_energy[:, 0]], color = "red", s = 4, label = f"Minimum B8: {signal_hypothesis[idx_energy[0,0]]:.3g}\nMinimum Tl208: {background_hypothesis[idx_energy[0,1]]:.3g}")
    axes[1].scatter(495.3*fixed_backg_weights[0], 66.3*signal_weight, color = "blue", label = "Background Model Assumption", s = 4)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax = cax)
    cbar.ax.tick_params(labelsize=5)
    cbar.ax.set_ylabel(r"$log(-2log(\mathcal{L}))$", fontsize = 5)
    # axes[1].set_title("Energy")
    axes[1].set_xlabel("Tl208 Normalisation", fontsize = 5)
    axes[1].set_ylabel("B8 Normalisation", fontsize = 5)
    axes[1].set_aspect("equal")
    axes[1].set_ylim((800, 0))
    axes[1].set_xlim((0,800))
    axes[1].legend(fontsize = 5)
    axes[1].tick_params(axis = "both", labelsize = 5)

    # create the profile log likelihood curve plots
    axes[2].plot(signal_hypothesis, profile_ll[2,:], color = "orange", label = "Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$")
    axes[2].tick_params(axis = "both", labelsize = 5)
    axes[2].set_xlabel(r"$^8$B Normalisation" ,fontsize = 5)
    axes[2].set_ylabel(r"$-2log(\mathcal{L})$", fontsize = 5)
    axes[2].set_ylim((0, 10))
    axes[2].axhline(1, linestyle = "dashed", color = "red")
    axes[2].legend(fontsize = 5)
    fig.tight_layout()
    plt.savefig(f"../plots/asimov_study/real_mc/advanced/diff_bins/combined_plot_{num_bins}.pdf")

# binning for the energy and multisite discriminant PDFs
energy_string     = "2p5_5p0"
mfactor = 1.25#2.5
uncert = 0.1
FV_CUT = 4500
# extra_name        = f"{FV_CUT}_unconstrained"
extra_name        = f"{FV_CUT}_constrained"
num_bins          = 1

energy_bins       = np.arange(2.5, 5.05, 0.1)
# energy_bins = np.linspace(2.5, 5.0, num_bins+1)

if energy_string == "2p5_3p75":
    multisite_bins = np.arange(-1.44, -1.38, 0.0005)
    fixed_backg_weights = [0.611, 0.554, 1, 1]
    signal_weight = 0.549
if energy_string == "3p5_5p0":
    if FV_CUT == 4500:
        multisite_bins = np.arange(-1.309, -1.26, 0.0005*mfactor)
    if FV_CUT == 5000:
        multisite_bins = np.arange(-1.180, -1.15, 0.0005*mfactor)
    signal_weight = 0.551
    fixed_backg_weights = [0.634, 0.620, 0, 0]
if energy_string == "3p75_5p0":
    # multisite_bins = np.arange(-1.375, -1.325, 0.0005)
    multisite_bins = np.arange(-1.04, -0.98, 0.0005)
    signal_weight = 0.451
    fixed_backg_weights = [0.389, 0.445, 0, 0]
if energy_string == "2p5_5p0":
    multisite_bins      = np.linspace(-1.355, -1.335, 50)# np.linspace(-2.12, -2.06, 100)#np.arange(-2.12, -2.06, 0.0005)#np.arange(-1.375, -1.325, 0.0005)
    fixed_backg_weights = [1, 1, 1, 1] # Tl208, Tl210, BiPo212, BiPo214
    signal_weight = 1
elif energy_string == "2p5_3p0":
    multisite_bins    = np.linspace(0.58, 0.63, 100)
    fixed_backg_weights = [0.135, 0.969, 0.986]
elif energy_string == "3p0_3p5":
    multisite_bins    = np.linspace(-1.60, -1.555, 100)
    fixed_backg_weights = [0.245, 0.031, 0.013]
elif energy_string == "3p5_4p0":
    multisite_bins    = np.linspace(-1.95, -1.90, 100)
    fixed_backg_weights = [0.345, 0.0, 0.001]
elif energy_string == "4p0_4p5":
    multisite_bins    = np.linspace(-1.074, -1.038, 100)
    fixed_backg_weights = [0.222, 0.0, 0.0]
elif energy_string == "4p5_5p0":
    multisite_bins    = np.linspace(1.780, 1.81, 25)
    fixed_backg_weights = [0.222, 0.0, 0.0]

backg_names       = ["Tl208", "Tl210", "BiPo212", "BiPo214"] 
labels            = ["B8"] + backg_names
mids_energy       = energy_bins[:-1] + np.diff(energy_bins)[0] / 2
mids_multi        = multisite_bins[:-1] + np.diff(multisite_bins)[0] / 2
signal_hypothesis = np.arange(0, 400, 1)
background_hypothesis = np.arange(974, 2374, 1)#np.arange(1, 500, 1)
analyse_real_data = False
generate_datasets = True
eval_bias         = True
expected_signal   = 91.8 * 365/145.7# signal_weight#66.3 * signal_weight#101.2

# expected_backg    = [169, 0 * fixed_backg_weights[1], 59.227 * fixed_backg_weights[2], 38.528 * fixed_backg_weights[3]] # Tl208, Tl210, BiPo212, BiPo214
# expected_backg    = [468.9 * fixed_backg_weights[0], 0 * fixed_backg_weights[1], 59.227 * fixed_backg_weights[2], 38.528 * fixed_backg_weights[3]] # Tl208, Tl210, BiPo212, BiPo214
# expected_backg    = [250.481 * fixed_backg_weights[0], 0 * fixed_backg_weights[1], 0 * fixed_backg_weights[2], 0 * fixed_backg_weights[3]] # Tl208, Tl210, BiPo212, BiPo214] 
expected_backg = [1174.61, 2.27, 148.38, 96.52]#[277.9, 0, 0, 0]
tl208_uncert   = 0.089 * expected_backg[0]#14.6#15.0#482.5 * fixed_backg_weights[0] * 0.089 # 8.9 % error

# 1 yr scaled livetime [166.084, 1174.61, 2.27, 148.38, 96.52]
# expected_signal = 166.084 * signal_weight
# expected_backg = np.array([1174.61, 2.27, 148.38, 96.52]) * np.array(fixed_backg_weights)
# expected_backg = expected_backg.tolist()
# expected_uncertainty = [10.15, 50.65, 1, 32.65, 19.25]
normalisations    = np.array([expected_signal] + expected_backg)
print("expected signal: ", expected_signal)
print("expected tl: ", expected_backg[0])
# can we load a pre-made pdf and data array or need to create from MC and ROOT files?
if generate_datasets == True:
    print("Generating PDFs and datasets.")
    dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array = create_pdfs_and_datasets(analyse_real_data, energy_string)
    print(dataset_multisite.shape)
    for i in range(multisite_pdf_array.shape[0]):
        print(f"Min: {np.min(multisite_pdf_array[i, :])} | Max: {np.max(multisite_pdf_array[i, :])}")
    
    if analyse_real_data == True:
        plot_name = f"real_{energy_string}_{extra_name}"
    else:
        plot_name = f"asimov_{energy_string}"
else:
    print("loading the already - computed PDF arrays...")
    energy_pdf_array    = np.load(f"./energy_pdf_array_{energy_string}.npy")
    multisite_pdf_array = np.load(f"./multisite_pdf_array_{energy_string}.npy")

    if analyse_real_data == False:
        dataset_energy    = np.load(f"./energy_dataset_asimov_{energy_string}.npy")
        dataset_multisite = np.load(f"./multisite_dataset_asimov_{energy_string}.npy")
        plot_name         = f"asimov_{energy_string}"
    else:
        print("loading real datasets")
        dataset_energy    = np.load(f"./energy_dataset_real_{energy_string}.npy")
        dataset_multisite = np.load(f"./multisite_dataset_real_{energy_string}.npy")
        plot_name         = f"real_{energy_string}_{extra_name}"

# blank out everything except from the Tl208 multisite PDF
# multisite_pdf_array[1:, :]  = 1e-12

# now, whatever data we have, perform the profile likelihood scan
profile_ll, norms, errors = profile_likelihood_scan(expected_backg[1:], expected_backg[0], dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array)

profile_ll[0, :] = rescale_ll(profile_ll[0, :])
profile_ll[1, :] = rescale_ll(profile_ll[1, :])
profile_ll[2, :] = rescale_ll(profile_ll[2, :])
# profile_ll = np.sqrt(profile_ll)

# save the profile LL
# np.save("./profileLL_constrained.npy", profile_ll)

# find 1 sigma frequentist confidence interval #
combined_error, multisite_error, energy_error = calculate_uncertainty(profile_ll)
# np.save("./error_constrained.npy", [combined_error, multisite_error, energy_error])
# create a plot of these
# plt.title("Profile Log-Likelihood Scan")
plt.plot(signal_hypothesis, profile_ll[2,:], color = "orange", label = "Energy | " + rf"${energy_error[0]}^{{+{energy_error[1]:.3g}}}_{{-{energy_error[2]:.3g}}}$")
plt.plot(signal_hypothesis, profile_ll[1, :], color = "green", label = "Multisite | " + rf"${multisite_error[0]:.3g}^{{+{multisite_error[1]:.3g}}}_{{-{multisite_error[2]:.3g}}}$")
plt.plot(signal_hypothesis, profile_ll[0, :], color = "black", label = "Combined | " + rf"${combined_error[0]:.3g}^{{+{combined_error[1]:.3g}}}_{{-{combined_error[2]:.3g}}}$")
plt.axhline(1.0, color = "red", linestyle = "dotted", label = r"1 $\sigma$ frequentist")
plt.axhline(4.0, color = "blue", linestyle = "dotted", label = r"2 $\sigma$ frequentist")
plt.axhline(9.0, color = "green", linestyle = "dotted", label = r"3 $\sigma$ frequentist")
# plt.axvline(expected_signal, color = "grey", marker = "^", markersize = 1, linestyle = "dashdot", label = rf"$^8$B Expectation: {expected_signal}")
plt.legend(frameon = False, prop=fm.FontProperties(fname=path2, size = 14))
plt.xlabel("Signal Hypothesis", fontproperties = prop_font)
plt.ylabel(r"$-2log(\mathcal{L})$", fontproperties = prop_font)
plt.ylim((0, 10))
plt.xlim((0,150))
ax = plt.gca()
for label in ax.get_xticklabels():
    label.set_fontproperties(prop_font)

for label in ax.get_yticklabels():
    label.set_fontproperties(prop_font)
# plt.savefig(f"../plots/asimov_study/real_mc/advanced/diff_bins/unconstrained_profileLL_{plot_name}_{num_bins}.pdf")
plt.savefig(f"../plots/asimov_study/real_mc/advanced/profileLL_{plot_name}.pdf")
plt.close()

# np.save("profile_ll_asimov_1yr_fullROI_2.npy", profile_ll)
# create a plot of the 3 fitted background + signal models vs data

# find the norms corresponding to the minimum LL value for each LL function
minimum_idx = np.argmin(profile_ll, axis = 1)

# find the normalisations fom energy, multisite and combined LL giving that minimum
norms_combined = norms[0, minimum_idx[0], :]
norms_multi    = norms[1, minimum_idx[1], :]
norms_energy   = norms[2, minimum_idx[2], :]
print(norms_energy)
# get the errors on the fitted norm of the Tl208
error_combined = errors[0, minimum_idx[0]]
error_multi    = errors[1, minimum_idx[1]]
error_energy   = errors[2, minimum_idx[2]]

# create the model out of the scaled PDFs
combined_energy_result    = norms_combined[:, None] * energy_pdf_array
combined_multisite_result = norms_combined[:, None] * multisite_pdf_array
energy_result             = norms_energy[:, None]   * energy_pdf_array
multisite_result          = norms_multi[:, None]    * multisite_pdf_array

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 15))
font_s    = 14
# create the 4 subplots showing relative contributions of each normalisation, the sum and the data
count = 0
for i in range(2):
    for j in range(2):
        ax = axes[i,j]

        if j == 0:
            data      = dataset_energy
            fit_bins  = energy_bins
            fit_error = error_combined
            data_mids = mids_energy
            xlims     = (2.5, 5.0)
             
            if i == 0:
                fit_result = combined_energy_result
                words      = [r"Energy Fit: Combined $\Delta log(\mathcal{L})$", "Reconstructed Energy (MeV)"]
            else:
                fit_result = energy_result
                words      = [r"Energy Fit", "Reconstructed Energy (MeV)"]
        else:
            data      = dataset_multisite
            fit_bins  = multisite_bins
            data_mids = mids_multi
            xlims     =(multisite_bins[0], multisite_bins[-1]) #(-2.09, -2.07)#(multisite_bins[0], multisite_bins[-1])  #(-1.29, multisite_bins[-1])#(multisite_bins[0], multisite_bins[-1]) 
            
            if i == 0:
                fit_result = combined_multisite_result
                fit_error = error_multi
                words      = [r"Multisite Fit: Combined $\Delta log(\mathcal{L})$", "Multisite Discriminant"]
            else:
                fit_result = multisite_result
                fit_error = error_energy
                words      = [r"Multisite Fit", "Multisite Discriminant"]

        print("Fit Result: ", fit_result)
        create_fit_model_subplot(ax, font_s, fit_result, data, fit_bins, data_mids, xlims, words, fit_error)

fig.tight_layout()
# plt.savefig(f"../plots/asimov_study/real_mc/advanced/diff_bins/fitted_background_model_{plot_name}_{num_bins}.pdf")
plt.savefig(f"../plots/asimov_study/real_mc/advanced/fitted_background_model_{plot_name}.pdf")
plt.close()

# convert the fitted B8 number to a flux
energy_flux, energy_positive_err, energy_negative_err          = extract_flux(norms_energy[0], energy_error[1], energy_error[2], livetime = 145.7, neutrino_detection_efficiency = 0.0641)
multisite_flux, multisite_positive_err, multisite_negative_err = extract_flux(norms_multi[0], multisite_error[1], multisite_error[2], livetime = 145.7, neutrino_detection_efficiency = 0.0641)
combined_flux, combined_positive_err, combined_negative_err    = extract_flux(norms_combined[0], combined_error[1], combined_error[2], livetime = 145.7, neutrino_detection_efficiency = 0.0641)

# extract the fitted Th322 g/g for each method
th232_conc_energy    = convert_fitted_tl_to_thorium(norms_energy[1]) 
th232_conc_multisite = convert_fitted_tl_to_thorium(norms_multi[1]) 
th232_conc_combined  = convert_fitted_tl_to_thorium(norms_combined[1])
bm_flux, bm_positive_err, bm_negative_err = extract_flux(30.7, combined_error[1], combined_error[2], livetime = 145.7, neutrino_detection_efficiency = 0.0641)
print(f"Th232 From Background Model Counts: {convert_fitted_tl_to_thorium(169.1):.3g} g/g")
print(f"B8 flux From Background Model Counts: {bm_flux:.3g} + {bm_positive_err:.3g} - {bm_negative_err:.3g}")
# background_model_LL_energy    = evaluate_loglikelihood(495.3, expected_backg[1:], 66.3, dataset_energy, energy_pdf_array)
# background_model_LL_multisite = evaluate_loglikelihood(495.3, expected_backg[1:], 66.3, dataset_multisite, multisite_pdf_array)
# background_model_LL_combined  = evaluate_loglikelihood(495.3, expected_backg[1:], 66.3, dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array)

flux_table = np.array([["Fit Method", "Fitted Flux", "+ve Err", "-ve Err", "Fitted Tl208", "Th232 Conc (g/g)"], ["Energy", f"{energy_flux:.3g}", f"+{energy_positive_err:.3g}", f"-{energy_negative_err:.3g}", f"{norms_energy[1]:.3g}", f"{th232_conc_energy:.3g}"], ["Multisite", f"{multisite_flux:.3g}", f"+{multisite_positive_err:.3g}", f"-{multisite_negative_err:.3g}", f"{norms_multi[1]:.3g}", f"{th232_conc_multisite:.3g}"], ["Combined", f"{combined_flux:.3g}", f"+{combined_positive_err:.3g}", f"-{combined_negative_err:.3g}", f"{norms_combined[1]:.3g}", f"{th232_conc_combined:.3g}"]])
tab = PrettyTable()
tab.field_names = flux_table[0, :]
for i in range(flux_table.shape[0] - 1):
    tab.add_row(flux_table[i+1])
print(tab)

l_space_energy = likelihood_gridsearch(expected_backg[1:], expected_backg[0], dataset_energy, dataset_multisite, energy_pdf_array, multisite_pdf_array)
# evaluate_bias()
# energy_bin_plots()