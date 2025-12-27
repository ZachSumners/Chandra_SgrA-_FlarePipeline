import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.optimize import curve_fit

from matplotlib import rcParams


df = pd.read_csv('flare_properties_dec8.csv')
flare_ids = df['obs_id']
max_rate = df['rate_max']
mjd_starts = df['start_mjd']
mjd_ends = df['end_mjd']
duration = df['duration_s']
fluence = df['fluence_ct']
flux = df['flux_1e_minus12']
luminosity = df['luminosity_1e34']
luminosity_err = df['luminosity_err']
fluence_err = df['fluence_err']
energy = df['energy_1e37erg']
energy_err = df['energy_err']

max_rate_err = df['rate_max_err']
mean_rate_err = df['rate_mean_err']
mean_rate = df['rate_mean']
shape = df['Shape']
hardnessratio = df['hardness_ratio']
hardnessratio_err = df['hardness_ratio_err']
#median_rate = df['median_rate']
#median_rate_err = df['median_rate_err']

correlation_AB = df['fluence_ct'].corr(df['rate_max'])
print(correlation_AB)

fluence_err = mean_rate_err * duration

mask = (fluence > 100)
duration_maxrates = duration[mask]
max_rate_highlighted = max_rate[mask]
fluence_mod = fluence[mask]

mask2 = (fluence > 500)
duration_maxrates2 = duration[mask2]
max_rate_highlighted2 = max_rate[mask2]
fluence_strong = fluence[mask2]

plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.size'] = 5

plt.rcParams.update({'text.usetex' : True, 'font.family' : 'Computer Modern Roman'})

fontsize = 14 
markersize = 6
capsize = 2
capthick = 0.6
elinewidth= 0.6
edgewidth = 0.6




mask0_duration = (fluence < 100) & (mean_rate > 0.03)
mask1_duration = (fluence > 100) & (fluence < 400) & (mean_rate > 0.03)
mask2_duration = (fluence > 400) & (mean_rate > 0.03)
mask3_duration = (fluence < 100) & (mean_rate < 0.03)
mask4_duration = (fluence > 100) & (fluence < 400) & (mean_rate < 0.03)
mask5_duration = (fluence > 400) & (mean_rate < 0.03)



plt.figure(figsize=(5, 4))
ratio = max_rate*duration/fluence

plt.errorbar(duration[mask0_duration], ratio[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
plt.errorbar(duration[mask1_duration], ratio[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
plt.errorbar(duration[mask2_duration], ratio[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
plt.errorbar(duration[mask3_duration], ratio[mask3_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(duration[mask4_duration], ratio[mask4_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(duration[mask5_duration], ratio[mask5_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth)
plt.axvline(6000, color='black', linestyle='--', linewidth=1)
plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Ratio')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("FlareCut_Ratio.pdf", dpi=300, bbox_inches='tight')
plt.close()






def power_law(x, A, alpha):
    return A * x**alpha

logx = np.log10(fluence)
logy = np.log10(max_rate)
sigma_logy = max_rate_err / (max_rate * np.log(10))

# linear model in log space: log y = log A + alpha * log x
def linear_model(logx, logA, alpha):
    return logA + alpha * logx

popt, pcov = curve_fit(linear_model, logx, logy, p0 = [1e-5, 1], sigma=sigma_logy, absolute_sigma=True)
logA, alpha_plot = popt
A_plot = 10**logA

print('Err Fluence, MaxRate = ', np.sqrt(np.diag(pcov)))

print(f"A Fluence, MaxRate= {A_plot:.4e}")
print(f"alpha Fluence, MaxRate = {alpha_plot:.4f}")

logy_pred = linear_model(logx, logA, alpha_plot)
residuals = logy - logy_pred
sigma_res = np.std(residuals, ddof=1)  # 1σ scatter of data around the fit

print(f"Intrinsic scatter in log10-space (sigma_res) = {sigma_res:.3f} dex")

k_pred = 1   # 1 ~ 68%, 2 ~ 95%, 3 ~ "extreme outliers"

logx_fit = np.linspace(logx.min(), logx.max(), 200)
logy_fit = linear_model(logx_fit, logA, alpha_plot)
logy_upper = logy_fit + k_pred * sigma_res
logy_lower = logy_fit - k_pred * sigma_res

x_fit = 10**logx_fit
y_fit = 10**logy_fit
y_upper_line = 10**logy_upper
y_lower_line = 10**logy_lower

logy_pred = linear_model(np.log10(fluence), logA, alpha_plot)

# upper and lower log-space bounds
logy_upper = logy_pred + k_pred * sigma_res
logy_lower = logy_pred - k_pred * sigma_res
y_upper = 10**logy_upper
y_lower = 10**logy_lower

cut_mask = (max_rate < y_lower) | (max_rate > y_upper)

mask0_duration = (fluence < 100) & ~cut_mask
mask1_duration = (fluence > 100) & (fluence < 400) & ~cut_mask
mask2_duration = (fluence > 400) & ~cut_mask
mask3_duration = (fluence < 100) & cut_mask
mask4_duration = (fluence > 100) & (fluence < 400) & cut_mask
mask5_duration = (fluence > 400) & cut_mask
fig, axs = plt.subplots(1, 1, figsize=(5, 4))

# scatter data
plt.errorbar(fluence[mask0_duration], max_rate[mask0_duration], yerr=max_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
plt.errorbar(fluence[mask1_duration], max_rate[mask1_duration], yerr=max_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
plt.errorbar(fluence[mask2_duration], max_rate[mask2_duration], yerr=max_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
plt.errorbar(fluence[mask3_duration], max_rate[mask3_duration], yerr=max_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(fluence[mask4_duration], max_rate[mask4_duration], yerr=max_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(fluence[mask5_duration], max_rate[mask5_duration], yerr=max_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth)

x_fit = np.logspace(np.log10(min(fluence)), np.log10(max(fluence)), 200)
y_fit = A_plot * x_fit**alpha_plot
plt.plot(x_fit, y_fit, color='black', linestyle='dashed', label=r'PL Fit ($x^{0.9}$)')
plt.fill_between(x_fit, y_lower_line, y_upper_line,
                    alpha=0.25, color='grey',
                    label=rf'{k_pred:.0f}$\sigma$ prediction band')
plt.xlabel('Fluence (ct)')
plt.ylabel('Max Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.legend()

#x_fit = np.linspace(min(duration), max(duration), 200)
#y_fit = A * x_fit**alpha
#axs[1].errorbar(fluence[~cut_mask], max_rate[~cut_mask], yerr=max_rate_err[~cut_mask], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='royalblue', markersize=markersize, markeredgewidth=edgewidth, label='Non Flagged Flares')
#axs[1].errorbar(fluence[cut_mask], max_rate[cut_mask], yerr=max_rate_err[cut_mask], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='lightgreen', markersize=markersize, markeredgewidth=edgewidth, label='Flagged Flares')
#axs[1].axvline(5800, color='black', linestyle='--', linewidth=1, label='5800s (Longest non-splittable flare)')
#axs[1].plot(x_fit, y_fit, color='red', label=r'PL Fit ($Ax^{-\alpha}$)')
#axs[1].set_xscale('log')
#axs[1].set_yscale('log')
#axs[1].set_xlabel('Fluence (s)')
#axs[1].set_ylabel('Max Rate (ct/s)')
#axs[1].legend()
plt.savefig("Flare_Cutting.pdf", dpi=300, bbox_inches='tight')

print(np.corrcoef(duration[~cut_mask], fluence[~cut_mask]))











def power_law(x, A, alpha):
    return A * x**alpha

popt, pcov = curve_fit(power_law, duration[~cut_mask], fluence[~cut_mask], sigma=fluence_err[~cut_mask], absolute_sigma=True)
A, alpha = popt

print(f"Duration Fluence A = {A:.4e}")
print(f"Duration Fluence alpha = {alpha:.4f}")
print(f'Err Duration Fluence {np.sqrt(np.diag(pcov))}')

plt.clf()
plt.figure(figsize=(5,4))


x_fit = np.linspace(min(duration), max(duration), 200)
y_fit = A * x_fit**alpha
plt.errorbar(duration[~cut_mask], fluence[~cut_mask], yerr=fluence_err[~cut_mask], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='cornflowerblue', markersize=markersize, markeredgewidth=edgewidth, label='Non Flagged Flares')
plt.errorbar(duration[cut_mask], fluence[cut_mask], yerr=fluence_err[cut_mask], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='gold', markersize=markersize, markeredgewidth=edgewidth, label='Flagged Flares')
#axs[1].axvline(5800, color='black', linestyle='--', linewidth=1, label='5800s (Longest non-splittable flare)')
plt.plot(x_fit, y_fit, color='black', linestyle='dashed', label=r'PL Fit ($x^{1.5}$)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Duration (s)')
plt.ylabel('Fluence (ct)')
plt.legend()
plt.savefig("Duration_Fluence_fit.pdf", dpi=300, bbox_inches='tight')




plt.clf()
fig, axs = plt.subplots(3, 1, figsize=(5, 12), sharex=True)
axs[0].errorbar(duration[mask0_duration], fluence[mask0_duration], yerr=fluence_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
axs[0].errorbar(duration[mask1_duration], fluence[mask1_duration], yerr=fluence_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
axs[0].errorbar(duration[mask2_duration], fluence[mask2_duration], yerr=fluence_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
axs[0].errorbar(duration[mask3_duration], fluence[mask3_duration], yerr=fluence_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth)
axs[0].errorbar(duration[mask4_duration], fluence[mask4_duration], yerr=fluence_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
axs[0].errorbar(duration[mask5_duration], fluence[mask5_duration], yerr=fluence_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth)
axs[0].plot(x_fit, y_fit, color='black', linestyle='dashed', label=r'PL Fit ($x^{1.5}$)')
axs[0].set_xlabel('Duration (s)')
axs[0].set_ylabel('Fluence (ct)')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend()
axs[0].set_xlim(300, 16000)
#plt.savefig("Duration_Fluence.pdf", dpi=300, bbox_inches='tight')





#plt.clf()

mask0_duration = (fluence < 100) & ~cut_mask
mask1_duration = (fluence > 100) & (fluence < 400) & ~cut_mask
mask2_duration = (fluence > 400) & ~cut_mask
mask3_duration = (fluence < 100) & cut_mask
mask4_duration = (fluence > 100) & (fluence < 400) & cut_mask
mask5_duration = (fluence > 400) & cut_mask


def power_law(x, A, alpha):
    return A * x**alpha

popt, pcov = curve_fit(power_law, duration[~cut_mask], max_rate[~cut_mask], sigma=max_rate_err[~cut_mask], absolute_sigma=True)
A, alpha = popt

print(f"Duration MaxRate A = {A:.4e}")
print(f"Duration MaxRate alpha = {alpha:.4f}")

print(f'Duration MaxRate p ~ {np.corrcoef(duration[~cut_mask], max_rate[~cut_mask])}')
print(f'Err Duration MaxRate {np.sqrt(np.diag(pcov))}')
#plt.clf()
#plt.figure(figsize=(5,4))


#fig, axs = plt.subplots(1, 1, figsize=(5, 4))

x_fit = np.linspace(min(duration), max(duration), 200)
y_fit = A * x_fit**alpha
axs[1].plot(x_fit, y_fit, color='black', linestyle='dashed', label=r'PL Fit ($x^{1.1}$)')
axs[1].errorbar(duration[mask0_duration], max_rate[mask0_duration], yerr=max_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
axs[1].errorbar(duration[mask1_duration], max_rate[mask1_duration], yerr=max_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
axs[1].errorbar(duration[mask2_duration], max_rate[mask2_duration], yerr=max_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
axs[1].errorbar(duration[mask3_duration], max_rate[mask3_duration], yerr=max_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth)
axs[1].errorbar(duration[mask4_duration], max_rate[mask4_duration], yerr=max_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
axs[1].errorbar(duration[mask5_duration], max_rate[mask5_duration], yerr=max_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth)
axs[1].set_xlabel('Duration (s)')
axs[1].set_ylabel('Max Rate (ct/s)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend()
#plt.savefig("Duration_MaxRate.pdf", dpi=300, bbox_inches='tight')








#fig, axs = plt.subplots(1, 1, figsize=(5, 4))
axs[2].errorbar(duration[mask0_duration], mean_rate[mask0_duration], yerr=mean_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
axs[2].errorbar(duration[mask1_duration], mean_rate[mask1_duration], yerr=mean_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
axs[2].errorbar(duration[mask2_duration], mean_rate[mask2_duration], yerr=mean_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
axs[2].errorbar(duration[mask3_duration], mean_rate[mask3_duration], yerr=mean_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='khaki', markersize=markersize, markeredgewidth=edgewidth)
axs[2].errorbar(duration[mask4_duration], mean_rate[mask4_duration], yerr=mean_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
axs[2].errorbar(duration[mask5_duration], mean_rate[mask5_duration], yerr=mean_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='crimson', markersize=markersize, markeredgewidth=edgewidth)
axs[2].set_xlabel('Duration (s)')
axs[2].set_ylabel('Mean Rate (ct/s)')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].legend()

plt.subplots_adjust(hspace=0)
plt.savefig("ParameterPlots.pdf", dpi=300, bbox_inches='tight')











mask0_duration = (fluence < 100) 
mask1_duration = (fluence > 100) & (fluence < 400) 
mask2_duration = (fluence > 400) 




plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(max_rate.min()), np.log10(max_rate.max()), 10)
plt.hist([max_rate[mask0_duration], max_rate[mask1_duration], max_rate[mask2_duration]],
    bins=bins,
    edgecolor='black',
    color=['khaki', 'darkorange', 'crimson'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares']
)
plt.xlabel('Max Flare Rate (ct/s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("MaxRateHist.pdf", dpi=300, bbox_inches='tight')
plt.close()



plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(duration.min()), np.log10(duration.max()), 10)
plt.hist([duration[mask0_duration], duration[mask1_duration], duration[mask2_duration]],
    bins=bins,
    edgecolor='black',
    color=['khaki', 'darkorange', 'crimson'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares']
)
plt.xlabel('Duration (s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("DurationHist.pdf", dpi=300, bbox_inches='tight')
plt.close()







plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(duration.min()), np.log10(duration.max()), 10)
datasets = [duration[mask0_duration], duration[mask3_duration], duration[mask1_duration], duration[mask4_duration], duration[mask2_duration],  duration[mask5_duration]]

labels = ['Weak Flares', '', 'Moderate Flares', '', 'Strong Flares', '']
colors = ['khaki', 'khaki', 'darkorange', 'darkorange', 'crimson', 'crimson']
hatches = ['', '///', '', '///', '', '///']

fig, ax = plt.subplots(figsize=(5, 4))
bottom = np.zeros(len(bins) - 1)

for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    ax.bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1, 60)
ax.set_xlabel('Duration (s)', fontsize=fontsize)
ax.set_ylabel('Frequency', fontsize=fontsize)
ax.legend()
plt.savefig("DurationHist.pdf", dpi=300, bbox_inches='tight')



plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(duration.min()), np.log10(duration.max()), 10)
datasets = [np.concatenate((duration[mask0_duration], duration[mask3_duration])), np.concatenate((duration[mask1_duration], duration[mask4_duration])), np.concatenate((duration[mask2_duration],  duration[mask5_duration]))]

labels = ['Weak Flares', 'Moderate Flares', 'Strong Flares']
colors = ['khaki', 'darkorange', 'crimson']
hatches = ['', '', '']

fig, ax = plt.subplots(figsize=(5, 4))
bottom = np.zeros(len(bins) - 1)

for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    ax.bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1, 60)
ax.set_xlabel('Duration (s)', fontsize=fontsize)
ax.set_ylabel('Frequency', fontsize=fontsize)
ax.legend()
plt.savefig("DurationHist_NoHatch.pdf", dpi=300, bbox_inches='tight')





plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(max_rate.min()), np.log10(max_rate.max()), 10)
datasets = [max_rate[mask0_duration], max_rate[mask3_duration], max_rate[mask1_duration], max_rate[mask4_duration], max_rate[mask2_duration],  max_rate[mask5_duration]]

labels = ['Weak Flares', '', 'Moderate Flares', '', 'Strong Flares', '']
colors = ['khaki', 'khaki', 'darkorange', 'darkorange', 'crimson', 'crimson']
hatches = ['', '///', '', '///', '', '///']

fig, ax = plt.subplots(figsize=(5, 4))
bottom = np.zeros(len(bins) - 1)

for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    ax.bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Max Rate (ct/s)', fontsize=fontsize)
ax.set_ylabel('Frequency', fontsize=fontsize)
ax.set_ylim(0, 70)
ax.legend()
plt.savefig("MaxRateHist.pdf", dpi=300, bbox_inches='tight')




plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(max_rate.min()), np.log10(max_rate.max()), 10)
datasets = [np.concatenate((max_rate[mask0_duration].to_numpy(), max_rate[mask3_duration].to_numpy())), np.concatenate((max_rate[mask1_duration].to_numpy(), max_rate[mask4_duration].to_numpy())), np.concatenate((max_rate[mask2_duration].to_numpy(), max_rate[mask5_duration].to_numpy()))]

labels = ['Weak Flares', 'Moderate Flares', 'Strong Flares']
colors = ['khaki', 'darkorange', 'crimson']
hatches = ['', '', '']

fig, ax = plt.subplots(figsize=(5, 4))
bottom = np.zeros(len(bins) - 1)

for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    ax.bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Max Rate (ct/s)', fontsize=fontsize)
ax.set_ylabel('Frequency', fontsize=fontsize)
ax.set_ylim(0, 70)
ax.legend()
plt.savefig("MaxRateHist_NoHatch.pdf", dpi=300, bbox_inches='tight')











fig, axs = plt.subplots(4, 1, figsize=(5, 12), sharex=True)

bins = np.linspace((hardnessratio.min()), (hardnessratio.max()), 10)
datasets = [hardnessratio[mask0_duration], hardnessratio[mask3_duration]]
labels = ['Weak Flares', '']
colors = ['khaki', 'khaki']
hatches = ['', '///']
bottom = np.zeros(len(bins) - 1)
for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    axs[0].bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist
axs[0].legend()

bins = np.linspace((hardnessratio.min()), (hardnessratio.max()), 10)
datasets = [hardnessratio[mask1_duration], hardnessratio[mask4_duration]]
labels = ['Moderate Flares', '']
colors = ['darkorange', 'darkorange']
hatches = ['', '///']
bottom = np.zeros(len(bins) - 1)
for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    axs[1].bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist
axs[1].legend()

bins = np.linspace((hardnessratio.min()), (hardnessratio.max()), 10)
datasets = [hardnessratio[mask2_duration], hardnessratio[mask5_duration]]
labels = ['Strong Flares', '']
colors = ['crimson', 'crimson']
hatches = ['', '///']
bottom = np.zeros(len(bins) - 1)
for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    axs[2].bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist
axs[2].set_ylim(0, 2.5)
axs[2].legend()


bins = np.linspace((hardnessratio.min()), (hardnessratio.max()), 10)
datasets = [hardnessratio[mask0_duration], hardnessratio[mask3_duration], hardnessratio[mask1_duration], hardnessratio[mask4_duration], hardnessratio[mask2_duration],  hardnessratio[mask5_duration]]
labels = ['Weak Flares', '', 'Moderate Flares', '', 'Strong Flares', '']
colors = ['khaki', 'khaki', 'darkorange', 'darkorange', 'crimson', 'crimson']
hatches = ['', '///', '', '///', '', '///']
bottom = np.zeros(len(bins) - 1)
for data, label, color, hatch in zip(datasets, labels, colors, hatches):
    hist, _ = np.histogram(data, bins=bins)
    axs[3].bar(bins[:-1], hist, width=np.diff(bins), bottom=bottom, align='edge', color=color, edgecolor='black', label=label, hatch=hatch)
    bottom += hist
axs[3].set_ylim(0, 32)
axs[3].legend()

plt.subplots_adjust(hspace=0)
plt.xlabel('Hardness Ratio')
plt.savefig("HardnessRatioHist_Multiplot.pdf", dpi=300, bbox_inches='tight')


print(f'Mean Weak: {np.mean(np.concatenate((hardnessratio[mask0_duration], hardnessratio[mask3_duration])))} +/- {1/ np.sqrt(np.sum(1/ np.concatenate((hardnessratio_err[mask0_duration], hardnessratio_err[mask3_duration]))**2))}')
print(f'Mean Moderate: {np.mean(np.concatenate((hardnessratio[mask1_duration], hardnessratio[mask4_duration])))} +/- {1/ np.sqrt(np.sum(1/ np.concatenate((hardnessratio_err[mask1_duration], hardnessratio_err[mask4_duration]))**2))}')
print(f'Mean Strong: {np.mean(np.concatenate((hardnessratio[mask2_duration], hardnessratio[mask5_duration])))} +/- {1/ np.sqrt(np.sum(1/ np.concatenate((hardnessratio_err[mask2_duration], hardnessratio_err[mask5_duration]))**2))}')



plt.figure(figsize=(5, 4))
singlepeak = (shape == 1) & ~cut_mask
doublepeak = (shape == 2) & ~cut_mask
higherpeak = (shape == 3) & ~cut_mask
mask_singlepeak = (shape == 1) & cut_mask
mask_doublepeak = (shape == 2) & cut_mask
mask_higherpeak = (shape == 3) & cut_mask
plt.errorbar(duration[singlepeak], fluence[singlepeak], yerr=fluence_err[singlepeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='yellowgreen', markersize=markersize, label='Single Peaked Flares')
plt.errorbar(duration[doublepeak], fluence[doublepeak], yerr=fluence_err[doublepeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='pink', markersize=markersize, label='Double Peaked Flares')
plt.errorbar(duration[higherpeak], fluence[higherpeak], yerr=fluence_err[higherpeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='mediumaquamarine', markersize=markersize, label='Complex Substructure')
plt.errorbar(duration[mask_singlepeak], fluence[mask_singlepeak], yerr=fluence_err[mask_singlepeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='yellowgreen', markersize=markersize)
plt.errorbar(duration[mask_doublepeak], fluence[mask_doublepeak], yerr=fluence_err[mask_doublepeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='pink', markersize=markersize)
plt.errorbar(duration[mask_higherpeak], fluence[mask_higherpeak], yerr=fluence_err[mask_higherpeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='mediumaquamarine', markersize=markersize)
#plt.axhline(y=100, color='black', linestyle='-')
#plt.axhline(y=400, color='black', linestyle='-')
plt.fill_between(np.linspace(0, 20000), 0, 100 ,alpha=0.2, color='khaki')
plt.fill_between(np.linspace(0, 20000), 100, 400 ,alpha=0.2, color='darkorange')
plt.fill_between(np.linspace(0, 20000), 400, 4000 ,alpha=0.2, color='crimson')
plt.xlim(300, 15000)
plt.ylim(10, 4000)
plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Fluence (ct)')
plt.xscale('log')
plt.yscale('log')
plt.savefig("Duration_Fluence_Shapes.pdf", dpi=300, bbox_inches='tight')









pl = [2.07, 2.20, 2.3, 2.62]
pl_err = [0.11, 0.14, 0.13, 0.14]
avg_fluence = [1965, 615.4, 291.5, 142.621]
avg_fluence_err = [43.36, 24.63, 16.9, 11.96]
avg_lumin = np.array([28.74, 9.04, 6.31, 3.61])# * 10**34
avg_lumin_err = np.array([0.044, 0.063, 0.0366, 0.0477])# * 10**34


plt.figure(figsize=(5, 4))
plt.errorbar(avg_lumin * 1e34, pl, xerr=avg_lumin_err * 1e34, yerr=pl_err, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='-o', mec='black', mfc='steelblue', markersize=markersize)
plt.xlabel('Average Luminosity of Flares Used to Calculate Power Law')
plt.ylabel('Power Law Index')
plt.xscale('log')
plt.yscale('log')
plt.savefig("Fluence_PL.pdf", dpi=300, bbox_inches='tight')
plt.close()







plt.figure(figsize=(5, 4))
plt.errorbar(duration, luminosity*1e34, yerr=luminosity_err*1e34, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.xlabel('Duration (s)')
plt.ylabel('Luminosity (erg/s)')
plt.xscale('log')
plt.yscale('log')
plt.savefig("Duration_Luminosity.pdf", dpi=300, bbox_inches='tight')

plt.figure(figsize=(5, 4))
plt.errorbar(duration, energy*1e37, yerr=energy_err*1e37, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.xlabel('Duration (s)')
plt.ylabel('Energy (erg)')
plt.xscale('log')
plt.yscale('log')
plt.savefig("Duration_Energy.pdf", dpi=300, bbox_inches='tight')



plt.clf()
plt.figure(figsize=(5, 4))
shape_sep = np.array([6, 4, 3, 5, 4, 7, 5, 7, 3, 4, 8, 6, 4, 5, 3, 6, 4, 3, 8, 6, 10, 4, 4, 6, 6, 5, 6, 5, 5, 5, 3, 3, 4, 8, 5, 5, 6, 6, 6])
#print('shape sep', durations/(shape_sep*300))
plt.hist(shape_sep * 300, bins=np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])*300, edgecolor='black', color='steelblue')
#plt.hist(max_rate, bins=bins, edgecolor='black', label='Weak Flares')
#plt.hist(max_rate_highlighted, bins=bins, edgecolor='black', color='orange', label='Moderate Flares')
#plt.hist(max_rate_highlighted2, bins=bins, edgecolor='black', color='red', label='Strong Flares')

plt.xlabel('Approximate Resolved Peak Separations (s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
#plt.xscale('log')
#plt.yscale('log')xx
plt.savefig("ShapeSepHist.pdf", dpi=300, bbox_inches='tight')
plt.close()

#plt.clf()
#plt.figure(figsize=(5, 4))
#plt.errorbar(duration, energy/duration, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
#plt.xlabel('Duration (s)')
#plt.ylabel('Energy/Duration (erg/s)')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()

print(f'Duration MeanRate {np.corrcoef(duration[~cut_mask], mean_rate[~cut_mask])}')