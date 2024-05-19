import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

"""
Create plots of the bias and pull as a function of signal fraction.
"""
f_size    = 14
# fraction  = np.arange(0.1, 1.0, 0.1)
# fraction  = np.arange(0.02, 0.22, 0.02)
fraction = np.linspace(0.01, 1.00, 20)
print(len(fraction))
bias_bins = np.arange(-2, +2, 0.1)

bias  = np.zeros((3, len(fraction)))
error = np.zeros((3, len(fraction)))

fig, axes = plt.subplots(nrows = 4, ncols = 5)
fig.set_size_inches((8.27, 11.69))
k = 0
j = 0
for i in range(len(fraction)):
#     print(j,k)
#     # load the bias arrays for this fraction
    x = np.load(f"./bias_and_pull_arrays/bias_{round(fraction[i], 3)}.npy")
#     print(x.shape)
#     # compute the mean and std of each dataset
    bias[:, i]  = np.mean(x, axis = 1) # should be 1D length (3)
    error[:, i] = np.std(x, axis = 1)  # should be 1D length (3)
#     print(x[0,:])
    print(j,k)
    axes[j,k].hist(x[0, :], bins = bias_bins, histtype = "step", color = "black")
    axes[j,k].hist(x[1, :], bins = bias_bins, histtype = "step", color = "green")
    axes[j,k].hist(x[2, :], bins = bias_bins, histtype = "step", color = "orange")
    axes[j,k].plot([], [], color = "black", label = "Combined")
    axes[j,k].plot([], [], color = "green", label = "Multisite")
    axes[j,k].plot([], [], color = "orange", label = "Energy")
    axes[j,k].legend(frameon = False, fontsize = 4)
    axes[j,k].set_title(f"{round(fraction[i], 3)}", fontsize = 4)
    axes[j,k].set_xlabel("Bias", fontsize = 4)
    axes[j,k].tick_params(which = "both", labelsize = 4)

    k += 1
    if k == 5:
        j += 1
        k = 0

plt.tight_layout()
plt.savefig("../plots/asimov_study/real_mc/advanced/bias_vs_signal_hists.pdf")
plt.close()


fig = plt.figure()
fig.set_size_inches(((8.27, 11.69*0.3)))
plt.errorbar(fraction, bias[0, :], yerr = error[0, :], linestyle = "", markersize = 2, capsize = 2, color = "black", marker = "o", label = "Combined")
plt.errorbar(fraction, bias[1, :], yerr = error[1, :], linestyle = "", markersize = 2, capsize = 2, color = "green", marker = "o", label = "Multisite")
plt.errorbar(fraction, bias[2, :], yerr = error[2, :], linestyle = "", markersize = 2, capsize = 2, color = "orange", marker = "o", label = "Energy")
plt.axhline(0, color = "black", linestyle = "dashed", linewidth = 1)
plt.legend(frameon = False, fontsize = 5)
plt.xlabel("S/N", fontsize = 5)
plt.ylabel("Bias", fontsize = 5)
plt.ylim((-2, 2))
ax =plt.gca()
ax.tick_params(which = "both", labelsize = 5)
plt.savefig("../plots/asimov_study/real_mc/advanced/bias_vs_signal.pdf")
plt.close()

