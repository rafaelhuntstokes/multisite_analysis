import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

"""
Script to plot out (nicely) the PDF information for the different energy bins
for the directionality and multisite.

The plots generated are: 
    1. Time Residual comparison for Tl208 and B8 (multisite PDFs)
    2. Directionality cos(theta_gamma) vs b8 tRes PDF
    3. Energy distributions of PDF'd events
    4. Position distributions for entire dataset
        a. rho2-Z
        b. X-Y
    5. Radius vs Energy histogram
"""
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

        

        return energies, multisites

def plot_pdf_analytics(pdf_energy_string):

    # use the string to define energy cuts
    if pdf_energy_string == "2.5_5.0":
        energy_low  = 2.5
        energy_high = 5.0
    if pdf_energy_string == "2.5_3.0":
        energy_low  = 2.5
        energy_high = 3.0
    if pdf_energy_string == "3.0_3.5":
        energy_low  = 3.0
        energy_high = 3.5
    if pdf_energy_string == "3.5_4.0":
        energy_low  = 3.5
        energy_high = 4.0
    if pdf_energy_string == "4.0_4.5":
        energy_low  = 4.0
        energy_high = 4.5
    if pdf_energy_string == "4.5_5.0":
        energy_low = 4.5
        energy_high = 5.0
    if pdf_energy_string == "2.5_3.75":
        energy_low = 2.5
        energy_high = 3.75
    if pdf_energy_string == "3.75_5.0":
        energy_low = 3.75
        energy_high = 5.0
    if pdf_energy_string == "3.5_5.0":
        energy_low = 3.5
        energy_high = 5.0
    
    # load the information
    pdf_B8_file    = ROOT.TFile.Open("../run_by_run_pdf/full_analysis2_B8_solar_nue/total.root")
    pdf_Tl208_file = ROOT.TFile.Open("../run_by_run_pdf/full_analysis2_Tl208/total.root")

    # extract the desired PDF histograms for re-plotting with matplotlib
    multisite_pdf_B8    = pdf_B8_file.Get(f"multi_{pdf_energy_string}")
    dir_pdf_B8          = pdf_B8_file.Get(f"directionality_{pdf_energy_string}")
    multisite_pdf_Tl208 = pdf_Tl208_file.Get(f"multi_{pdf_energy_string}")
    dir_pdf_Tl208       = pdf_Tl208_file.Get(f"directionality_{pdf_energy_string}")
    
    # get the bin edges and contents for each PDF
    multi_B8_bins      = [multisite_pdf_B8.GetBinLowEdge(i) for i in range(1, multisite_pdf_B8.GetNbinsX() + 2)]
    multi_B8_counts    = [multisite_pdf_B8.GetBinContent(i) / multisite_pdf_B8.Integral() for i in range(1, multisite_pdf_B8.GetNbinsX()+1)]
    multi_Tl208_bins   = [multisite_pdf_Tl208.GetBinLowEdge(i) for i in range(1, multisite_pdf_Tl208.GetNbinsX() + 2)]
    multi_Tl208_counts = [multisite_pdf_Tl208.GetBinContent(i) / multisite_pdf_Tl208.Integral() for i in range(1, multisite_pdf_Tl208.GetNbinsX()+1)]

    np.save(f"/home/hunt-stokes/Tl208_tRes_normed_4p5mFV_{pdf_energy_string}.npy", multi_Tl208_counts)
    np.save(f"/home/hunt-stokes/B8_tRes_normed_4p5mFV_{pdf_energy_string}.npy", multi_B8_counts)
    np.save(f"/home/hunt-stokes/binning_4p5mFV_{pdf_energy_string}.npy", multi_B8_bins)

    # create numpy arrays to save the directionality 2D PDF bins
    directionality_pdf = np.zeros((dir_pdf_B8.GetNbinsX(), dir_pdf_B8.GetNbinsY()))
    for i in range(1, dir_pdf_B8.GetNbinsX() + 1):
        for j in range(1, dir_pdf_B8.GetNbinsY() + 1):
            directionality_pdf[i-1, j-1] = dir_pdf_B8.GetBinContent(i,j)
    x_min = dir_pdf_B8.GetXaxis().GetXmin()
    x_max = dir_pdf_B8.GetXaxis().GetXmax()
    y_min = dir_pdf_B8.GetYaxis().GetXmin()
    y_max = dir_pdf_B8.GetYaxis().GetXmax()

    # get the event level info for plotting the position and energy distributions
    event_info_B8    = pdf_B8_file.Get("PDF_event_by_event_Information")
    event_info_Tl208 = pdf_Tl208_file.Get("PDF_event_by_event_Information")
    
    num_events_B8    = event_info_B8.GetEntries()
    num_events_Tl208 = event_info_Tl208.GetEntries()
    energy_B8        = np.zeros(num_events_B8)
    rho2_B8          = np.zeros(num_events_B8)
    z_B8             = np.zeros(num_events_B8)
    x_B8             = np.zeros(num_events_B8)
    y_B8             = np.zeros(num_events_B8)
    r_B8             = np.zeros(num_events_B8)
    energy_Tl208     = np.zeros(num_events_Tl208)
    rho2_Tl208       = np.zeros(num_events_Tl208)
    z_Tl208          = np.zeros(num_events_Tl208)
    x_Tl208          = np.zeros(num_events_Tl208)
    y_Tl208          = np.zeros(num_events_Tl208)
    r_Tl208          = np.zeros(num_events_Tl208)
    counter = 0
    for ientry in event_info_B8:
        e = ientry.energy
        if e >= energy_low and e < energy_high:
            energy_B8[counter] = e

            X = ientry.x
            Y = ientry.y
            Z = ientry.z
            rho2_B8[counter] = (X**2 + Y**2) / 6000**2
            z_B8[counter]    = Z
            x_B8[counter]    = X
            y_B8[counter]    = Y
            r_B8[counter]    = np.sqrt(X**2 + Y**2 + Z**2)
        counter += 1
    counter = 0
    for ientry in event_info_Tl208:
        e = ientry.energy
        if e >= energy_low and e < energy_high:
            energy_Tl208[counter] = e
            #itr = ientry.itr
            X = ientry.x
            Y = ientry.y
            Z = ientry.z
            rho2_Tl208[counter] = (X**2 + Y**2) / 6000**2
            z_Tl208[counter]    = Z
            x_Tl208[counter]    = X
            y_Tl208[counter]    = Y
            r_Tl208[counter]    = np.sqrt(X**2 + Y**2 + Z**2)
        counter += 1
    
    # remove all the zeros (if any) left in the arrays
    rho2_B8 = rho2_B8[np.nonzero(rho2_B8)]
    z_B8 = z_B8[np.nonzero(z_B8)]
    x_B8 = x_B8[np.nonzero(x_B8)]
    y_B8 = y_B8[np.nonzero(y_B8)]
    energy_B8 = energy_B8[np.nonzero(energy_B8)]
    r_B8 = r_B8[np.nonzero(r_B8)]

    rho2_Tl208 = rho2_Tl208[np.nonzero(rho2_Tl208)]
    z_Tl208 = z_Tl208[np.nonzero(z_Tl208)]
    x_Tl208 = x_Tl208[np.nonzero(x_Tl208)]
    y_Tl208 = y_Tl208[np.nonzero(y_Tl208)]
    energy_Tl208 = energy_Tl208[np.nonzero(energy_Tl208)]
    r_Tl208 = r_Tl208[np.nonzero(r_Tl208)]

    fig, axes = plt.subplots(nrows = 4, ncols = 2)
    fig.set_size_inches((8.27, 11.69))
    font_s = 8
    # multisite PDF
    axes[0,0].step(multi_B8_bins[:-1], multi_B8_counts,  where = 'post', label = r"$^8B$" + f" | Num. Events: {len(energy_B8)}")
    axes[0,0].step(multi_B8_bins[:-1], multi_Tl208_counts,  where = 'post', label = r"$^{208}Tl$"+ f" | Num. Events: {len(energy_Tl208)}")
    axes[0,0].legend(frameon=False, fontsize = font_s)
    axes[0,0].set_xlabel("Time Residual (ns)", fontsize = font_s)
    axes[0,0].set_ylabel("Normalised Counts per 1 ns Bin", fontsize = font_s)
    axes[0,0].set_xlim((-5, 100))

    # directionality PDF
    # axes[0,1].imshow(directionality_pdf.T, origin = 'lower', aspect = 'auto', extent = [x_min, x_max, y_min, y_max])
    # axes[0,1].set_xlabel("Time Residual (ns)", fontsize = font_s)
    # axes[0,1].set_ylabel(r"$cos(\theta _\gamma)$", fontsize = font_s)
    # axes[0,1].set_xlim((-10, 80))

    # energy distribution
    width = 0.05
    energy_bins = np.arange(2.5, 5.0 + width, width)
    axes[0,1].plot([], [], label = r"$^8B$" + f" | Num. Events: {len(energy_B8)}")
    axes[0,1].plot([], [], label = r"$^{208}Tl$" + f" | Num. Events: {len(energy_Tl208)}")
    axes[0,1].hist(energy_B8, bins = energy_bins, density = True, color = "C0", histtype = "step")
    axes[0,1].hist(energy_Tl208, bins = energy_bins, density = True, color = "C1", histtype = "step")
    axes[0,1].set_xlabel("Reconstructed Energy (MeV)", fontsize = font_s)
    axes[0,1].set_ylabel(f"Normalised Counter per {width} MeV", fontsize = font_s)
    axes[0,1].legend(frameon=False, fontsize = font_s)
    
    # energy vs radius
    r_bins = np.linspace(0, 6000, 50)
    axes[1,0].hist2d(r_B8, energy_B8, bins = [r_bins, energy_bins], cmin = 1e-4, density = False)
    axes[1,0].set_xlabel("Reconstructed Radius (mm)", fontsize = font_s)
    axes[1,0].set_ylabel("Reconstructed Energy (MeV)", fontsize = font_s)
    axes[1,0].plot([], [], linestyle = "", label = r"$^8B$")
    axes[1,0].legend(frameon=False, fontsize = font_s)

    axes[1,1].hist2d(r_Tl208, energy_Tl208, bins = [r_bins, energy_bins], cmin = 1e-4, density = False)
    axes[1,1].set_xlabel("Reconstructed Radius (mm)", fontsize = font_s)
    axes[1,1].set_ylabel("Reconstructed Energy (MeV)", fontsize = font_s)
    axes[1,1].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[1,1].legend(frameon = False, fontsize = font_s)

    # rho2_z position distribution
    z_bins    = np.linspace(-6000, 6000, 100)
    rho2_bins = np.linspace(0, 1, 100)
    axes[2,0].hist2d(rho2_B8, z_B8, bins = [rho2_bins, z_bins], density = False, cmin = 1e-4)
    axes[2,0].plot([], [], linestyle = "", label = r"$^8B$")
    axes[2,0].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontsize = font_s)
    axes[2,0].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[2,0].legend(frameon=False, fontsize = font_s)

    axes[2,1].hist2d(rho2_Tl208, z_Tl208, bins = [rho2_bins, z_bins], density = False, cmin = 1e-4)
    axes[2,1].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[2,1].set_xlabel(r"$\left(\frac{\rho}{\rho_{AV}}\right)^2$", fontsize = font_s)
    axes[2,1].set_ylabel("Reconstructed Z (mm)", fontsize = font_s)
    axes[2,1].legend(frameon=False, fontsize = font_s)

    # X-Y position plot
    axes[3,0].hist2d(x_B8, y_B8, bins = [z_bins, z_bins], density = False, cmin = 1e-10)
    axes[3,0].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[3,0].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[3,0].plot([], [], linestyle = "", label = r"$^8B$")
    axes[3,0].legend(frameon = False, fontsize = font_s)

    axes[3,1].hist2d(x_Tl208, y_Tl208, bins = [z_bins, z_bins], density = False, cmin = 1e-10)
    axes[3,1].set_xlabel("Reconstructed X (mm)", fontsize = font_s)
    axes[3,1].set_ylabel("Reconstructed Y (mm)", fontsize = font_s)
    axes[3,1].plot([], [], linestyle = "", label = r"$^{208}Tl$")
    axes[3,1].legend(frameon = False, fontsize = font_s)
    
    fig.tight_layout()
    # plt.subplots_adjust(top=0.93)
    # plt.suptitle(f"PDF Information: {pdf_energy_string} MeV | B8 Evs: {len(x_B8)} | Tl208 Evs: {len(x_Tl208)}", fontsize = 12)
    plt.savefig(f"../plots/pdf_analytics/{pdf_energy_string}.pdf")
    
def plot_pdf_energy_bins():
    """
    Create a simple plot showing the energy bins used in the analysis overlaid on the
    full MC event energy distributions.
    """
    # load the information
    pdf_B8_file    = ROOT.TFile.Open("../run_by_run_pdf/B8_solar_nue/total.root")
    pdf_Tl208_file = ROOT.TFile.Open("../run_by_run_pdf/Tl208/total.root")
    
    # get the event level info for plotting the position and energy distributions
    event_info_B8    = pdf_B8_file.Get("PDF_event_by_event_Information")
    event_info_Tl208 = pdf_Tl208_file.Get("PDF_event_by_event_Information")
    
    # arrays to save the energy spectrum
    energy_Tl208 = []
    energy_B8    = []
    
    for ientry in event_info_B8:
        # print("B8: ", ientry)
        energy_B8.append(ientry.energy)
    for ientry in event_info_Tl208:
        # print("Tl208: ", ientry)
        energy_Tl208.append(ientry.energy)
        
    # create the plots with energy bins overlaid
    energy_bin_lines = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    energy_binning = np.arange(2.5, 5.5, 0.5)
    
    plt.figure(figsize = (8,4))
    counts_b8, _, _     = plt.hist(energy_B8, bins = energy_binning, histtype = "step", density = False, label = r"$^8B$")
    counts_tl208, _, _ = plt.hist(energy_Tl208, bins = energy_binning, histtype = "step", density = False, label = r"$^{208}Tl$")
    print(counts_tl208)
    for line in energy_bin_lines:
        plt.axvline(line, color = "black")
    
    plt.legend()
    plt.xlim((2.5, 5.0))
    plt.xlabel("Reconstructed Energy (MeV)")
    plt.ylabel("Counts per 0.05 MeV")
    plt.savefig("../plots/pdf_analytics/pdf_energy_binning.pdf")

    tot_b8 = np.sum(counts_b8)
    tot_tl = np.sum(counts_tl208)
    print("Percentage of B8 events in each energy bin:\n")
    print(f"{counts_b8[0]/tot_b8}, {counts_b8[1]/tot_b8}, {counts_b8[2]/tot_b8}, {counts_b8[3]/tot_b8}")
    print("\n")
    print("Percentage of Tl events in each energy bin:\n")
    print(f"{counts_tl208[0]/tot_tl}, {counts_tl208[1]/tot_tl}, {counts_tl208[2]/tot_tl}, {counts_tl208[3]/tot_tl}")
    print("\n")
    
def plot_multisite_energy_pdfs():
    """
    Create a nice 2 plot subplot showing normalised multisite and energy PDFs
    for each isotope used in the energy + multisite fits.
    """
    pdf_runlist    = np.loadtxt("../runlists/full_test.txt", dtype = int)
    energy_bins    = np.arange(2.5, 5.05, 0.05)
    multisite_bins = np.arange(-1.36, -1.325, 0.0005)
    energy_string  = "2p5_5p0"
    path           = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies/run_by_run_test"
    
    # load the information
    b8_energies, b8_multisites           = obtain_pdf(f"{path}/full_analysis_B8_solar_nue", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    tl208_energies, tl208_multisites     = obtain_pdf(f"{path}/full_analysis_Tl208", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    tl210_energies, tl210_multisites     = obtain_pdf(f"{path}/full_analysis_Tl210", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    bipo214_energies, bipo214_multisites = obtain_pdf(f"{path}/full_analysis_BiPo214", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    bipo212_energies, bipo212_multisites = obtain_pdf(f"{path}/full_analysis_BiPo212", 4500.0, multisite_bins, energy_bins, pdf_runlist, energy_string, "B8_nue") 
    
    # create the plots of the pdfs
    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    fig.set_size_inches((8.27, 11.69*0.3))
    axes[0].hist(b8_energies,      bins = energy_bins, histtype = "step", density = True)
    axes[0].hist(tl208_energies,   bins = energy_bins, histtype = "step", density = True)
    axes[0].hist(tl210_energies,   bins = energy_bins, histtype = "step", density = True)
    axes[0].hist(bipo214_energies, bins = energy_bins, histtype = "step", density = True)
    axes[0].hist(bipo212_energies, bins = energy_bins, histtype = "step", density = True)
    axes[0].set_xlabel("Reconstructed Energy (MeV)", fontsize = 8)
    axes[0].plot([], [], color = "C0", label = r"$^8$B $\nu _e$" + f" | Num. Events: {len(b8_energies)}")
    axes[0].plot([], [], color = "C1", label = r"$^{208}$Tl" + f" | Num. Events: {len(tl208_energies)}")
    axes[0].plot([], [], color = "C2", label = r"$^{210}$Tl" + f" | Num. Events: {len(tl210_energies)}")
    axes[0].plot([], [], color = "C3", label = r"BiPo214" + f" | Num. Events: {len(bipo214_energies)}")
    axes[0].plot([], [], color = "C4", label = r"BiPo212" + f" | Num. Events: {len(bipo212_energies)}")
    axes[0].legend(fontsize = 6, frameon = False)
    axes[0].tick_params(axis='both', labelrotation=0, labelsize = 8)
    

    axes[1].hist(b8_multisites,      bins = multisite_bins, histtype = "step", density = True)
    axes[1].hist(tl208_multisites,   bins = multisite_bins, histtype = "step", density = True)
    axes[1].hist(tl210_multisites,   bins = multisite_bins, histtype = "step", density = True)
    axes[1].hist(bipo214_multisites, bins = multisite_bins, histtype = "step", density = True)
    axes[1].hist(bipo212_multisites, bins = multisite_bins, histtype = "step", density = True)
    axes[1].set_xlabel(r"$\Delta log(\mathcal{L})$", fontsize = 8)
    axes[1].plot([], [], color = "C0", label = r"$^8$B $\nu _e$" + f" | Num. Events: {len(b8_multisites)}")
    axes[1].plot([], [], color = "C1", label = r"$^{208}$Tl" + f" | Num. Events: {len(tl208_multisites)}")
    axes[1].plot([], [], color = "C2", label = r"$^{210}$Tl" + f" | Num. Events: {len(tl210_multisites)}")
    axes[1].plot([], [], color = "C3", label = r"BiPo214" + f" | Num. Events: {len(bipo214_multisites)}")
    axes[1].plot([], [], color = "C4", label = r"BiPo212" + f" | Num. Events: {len(bipo212_multisites)}")
    axes[1].legend(fontsize = 6, frameon = False)
    axes[1].tick_params(axis='both', labelrotation=0, labelsize = 8)

    fig.tight_layout()
    plt.savefig(f"../plots/pdf_analytics/energy_multisites_{energy_string}.pdf")

    print(len(b8_energies), len(b8_multisites))

def plot_multisite_solo_pdfs():
    """
    Plot the multisite PDFs for Tl208 and B8 scaled to the total counts in the
    dataset. This demonstrates that the multisite fits pretty well with Tl208
    alone, and gives a reason why the B8 fit is poor.
    """

    binning = np.arange(-1.375, -1.325, 0.0005)
    mids    = binning[:-1] + np.diff(binning)[0] / 2
    data    = np.load("multisite_dataset_real_2p5_3p0.npy")
    energies = ["2p5_3p0", "3p0_3p5", "3p5_4p0", "4p0_4p5"]
    pdfs    = np.load("multisite_pdf_array_2p5_3p0.npy")

    tl208_pdf = pdfs[1, :] * np.sum(data)
    b8_pdf    = pdfs[0, :] * np.sum(data)

    fig, axes = plt.subplots(nrows = 1, ncols = 2)
    for energy in energies:
        print(energy)
        pdfs    = np.load(f"multisite_pdf_array_{energy}.npy")
        tl208_pdf = pdfs[1, :] * np.sum(data)
        b8_pdf    = pdfs[0, :] * np.sum(data)
        # axes[1].errorbar(mids, data, yerr = np.sqrt(data), color = "black", capsize = 2, linestyle = "", marker = "^", markersize = 4, label = "Data")
        axes[1].step(binning[:-1], b8_pdf.tolist() + [0], where = "post", linewidth = 2, color = "C0", label = r"$^{8}$B")
        axes[1].legend(frameon = False, fontsize = 5)
        axes[1].tick_params(axis = "both", labelsize = 5)
        axes[1].set_xlabel("Multisite Discriminant", fontsize = 5)
        axes[1].set_ylabel("Counts", fontsize = 5)
        # axes[1].set_ylim((0, 90))

        # axes[0].errorbar(mids, data, yerr = np.sqrt(data), color = "black", capsize = 2, linestyle = "", marker = "^", markersize = 4, label = "Data")
        axes[0].step(binning[:-1], tl208_pdf.tolist() + [0], where = "post", linewidth = 2, color = "red", label = r"$^{208}$Tl")
        axes[0].legend(frameon = False, fontsize = 5)
        axes[0].tick_params(axis = "both", labelsize = 5)
        axes[0].set_xlabel("Multisite Discriminant", fontsize = 5)
        axes[0].set_ylabel("Counts", fontsize = 5)
        # axes[0].set_ylim((0, 90))
    fig.set_size_inches(((8.27, 11.69 *0.3)))
    fig.tight_layout()

    # plt.savefig("../plots/asimov_study/real_mc/advanced/diff_bins/tl208_b8_multisite_only.pdf")
    plt.savefig("test.pdf")

plot_multisite_solo_pdfs()
# plot_pdf_energy_bins()
    
# plot_pdf_analytics("2.5_5.0")
# plot_pdf_analytics("2.5_3.0")
# plot_pdf_analytics("3.0_3.5")
# plot_pdf_analytics("3.5_4.0")
# plot_pdf_analytics("4.0_4.5")
# plot_pdf_analytics("4.5_5.0")
# plot_pdf_analytics("2.5_3.75")
# plot_pdf_analytics("3.75_5.0")
# plot_pdf_analytics("3.5_5.0")
# plot_pdf_analytics("2.5_3.125")
# plot_pdf_analytics("3.125_3.75")
# plot_pdf_analytics("3.75_4.375")
# plot_pdf_analytics("4.375_5.0")

# plot_multisite_energy_pdfs()