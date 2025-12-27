import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess

from matplotlib import rcParams

df = pd.read_csv('flare_properties_dec2_final.csv')
flare_ids = df['obs_id']
max_rate = df['rate_max']
mjd_starts = df['start_mjd']
mjd_ends = df['end_mjd']
duration = df['duration_s']
fluence = df['fluence_ct']
flux = df['flux_1e_minus12']
luminosity = df['luminosity_1e34']
fluence_err = df['fluence_err']

max_rate_err = df['rate_max_err']
mean_rate_err = df['rate_mean_err']
mean_rate = df['rate_mean']
shape = df['Shape']
hardnessratio = df['hardness_ratio']
hardnessratio_err = df['hardness_ratio_err']
#median_rate = df['median_rate']
#median_rate_err = df['median_rate_err']

#fluence_err = mean_rate_err * duration


fontsize = 12

plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.titlesize'] = 10


plt.rcParams['xtick.major.size'] = 8  # length of major ticks
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.size'] = 5

plt.rcParams.update({'text.usetex' : True, 'font.family' : 'Computer Modern Roman'})

fontsize = 14  # or whatever you use





mask0_duration = (fluence < 100) & (mean_rate > 0.03)
mask1_duration = (fluence > 100) & (fluence < 400) & (mean_rate > 0.03)
mask2_duration = (fluence > 400) & (mean_rate > 0.03)
mask3_hist = (mean_rate < 0.03)



plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(max_rate.min()), np.log10(max_rate.max()), 10)
plt.hist([max_rate[mask0_duration], max_rate[mask1_duration], max_rate[mask2_duration], max_rate[mask3_hist]],
    bins=bins,
    edgecolor='black',
    color=['khaki', 'darkorange', 'crimson'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares', 'Flares mean $<$ 0.03ct/s']
)
#plt.hist(max_rate, bins=bins, edgecolor='black', label='Weak Flares')
#plt.hist(max_rate_highlighted, bins=bins, edgecolor='black', color='orange', label='Moderate Flares')
#plt.hist(max_rate_highlighted2, bins=bins, edgecolor='black', color='red', label='Strong Flares')

plt.xlabel('Max Flare Rate (ct/s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("MaxRateHist.pdf", dpi=300, bbox_inches='tight')
plt.close()






plt.clf()
plt.figure(figsize=(5, 4))
shape_sep = np.array([3, 5, 5, 5, 4, 5, 6, 6, 6])
durations = np.array([3336, 3555, 6405, 3290, 6437, 10483, 5538, 5986, 5358])

shape_sep = np.array([6, 4, 3, 5, 4, 7, 5, 7, 3, 4, 8, 6, 4, 5, 3, 6, 4, 3, 8, 6, 10, 4, 4, 6, 6, 5, 6, 5, 5, 5, 3, 3, 4, 8, 5, 5, 6, 6, 6])
#print('shape sep', durations/(shape_sep*300))
plt.hist(shape_sep * 300, bins=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], edgecolor='black', color='steelblue')
#plt.hist(max_rate, bins=bins, edgecolor='black', label='Weak Flares')
#plt.hist(max_rate_highlighted, bins=bins, edgecolor='black', color='orange', label='Moderate Flares')
#plt.hist(max_rate_highlighted2, bins=bins, edgecolor='black', color='red', label='Strong Flares')

plt.xlabel('Approximate Resolved Peak Separations (s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.savefig("ShapeSepHist.pdf", dpi=300, bbox_inches='tight')
plt.close()






plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(duration.min()), np.log10(duration.max()), 10)
plt.hist([duration[mask0_duration], duration[mask1_duration], duration[mask2_duration], duration[mask3_hist]],
    bins=bins,
    edgecolor='black',
    color=['steelblue', 'orange', 'red', 'gray'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares', 'Flares mean $<$ 0.03ct/s']
)
#plt.hist(duration[mask0_duration], bins=bins, edgecolor='black', color='steelblue', label='Weak Flares', stacked=True)
#plt.hist(duration[mask1_duration], bins=bins, edgecolor='black', color='orange', label='Moderate Flares', stacked=True)
#plt.hist(duration[mask2_duration], bins=bins, edgecolor='black', color='red', label='Strong Flares', stacked=True)
#plt.hist(duration[mask3_duration], bins=bins, edgecolor='black', color='steelblue', label='Weak Flares', hatch='///', stacked=True)
#plt.hist(duration[mask4_duration], bins=bins, edgecolor='black', color='orange', label='Moderate Flares', hatch='///', stacked=True)
#plt.hist(duration[mask5_duration], bins=bins, edgecolor='black', color='red', label='Strong Flares', hatch='///', stacked=True)
plt.xlabel('Flare Duration (s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("DurationHist.pdf", dpi=300, bbox_inches='tight')
plt.close()


plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(fluence.min()*2), np.log10(fluence.max()), 10)

plt.hist([fluence[mask0_duration], fluence[mask1_duration], fluence[mask2_duration], fluence[mask3_hist]],
    bins=bins,
    edgecolor='black',
    color=['steelblue', 'orange', 'red', 'gray'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares', 'Flares mean $<$ 0.03ct/s']
)
#plt.hist(fluence, bins=bins, edgecolor='black', label='Weak Flares')
#plt.hist(fluence[mask], bins=bins, edgecolor='black', color='orange', label='Moderate Flares')
#plt.hist(fluence[mask2], bins=bins, edgecolor='black', color='red', label='Strong Flares')
#plt.plot(np.arange(10, 4000), 67*np.arange(10, 4000)**-0.8)
plt.xlabel('Fluence (ct)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("FluenceHist.pdf", dpi=300, bbox_inches='tight')
plt.close()


fluence_mask = fluence > 100

#for flare in flare_ids[fluence_mask]:
#    print(f'FLARE {flare}')
#    flare_5digit = str(flare).zfill(5)
#    subprocess.call(f'fitsheader ./{flare}/repro/acisf{flare_5digit}_bary_evt2.fits | grep DETNAM', shell=True, cwd='/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub')
#    subprocess.call(f'fitsheader ./{flare}/repro/acisf{flare_5digit}_bary_evt2.fits | grep EXPTIME', shell=True, cwd='/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub')

#df_fluence_cut = pd.DataFrame(data = {'id': flare_ids[fluence_mask], 'start_mjd': mjd_starts[fluence_mask].values, 'end_mjd': mjd_ends[fluence_mask].values, 'fluence': fluence[fluence_mask].values})
#df_fluence_cut.to_csv('fluence_cut.csv')

markersize = 6
capsize = 2
capthick = 0.6
elinewidth= 0.6
edgewidth = 0.6


def linear_model(logx, logA, alpha):
    return logA + alpha * logx

logy_pred = linear_model(np.log10(fluence), -2.7136256638552467, 0.8554)

k_pred = 1
sigma_res = 0.217
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

plt.figure(figsize=(5, 4))

ratio = mean_rate/max_rate

plt.errorbar(duration[mask0_duration], ratio[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
plt.errorbar(duration[mask1_duration], ratio[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
plt.errorbar(duration[mask2_duration], ratio[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
plt.errorbar(duration[mask3_duration], ratio[mask3_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(duration[mask4_duration], ratio[mask4_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(duration[mask5_duration], ratio[mask5_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth)

plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Ratio')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("MaxRate_Fluence_Ratio.pdf", dpi=300, bbox_inches='tight')
plt.close()


'''
plt.figure(figsize=(5, 4))
bins = np.logspace(np.log10(ratio.min()), np.log10(ratio.max()), 10)
plt.hist(ratio, bins=bins, edgecolor='black', label='Weak Flares')
plt.hist(ratio[mask], bins=bins, edgecolor='black', color='orange', label='Moderate Flares')
plt.hist(ratio[mask2], bins=bins, edgecolor='black', color='red', label='Strong Flares')
plt.xlabel('Max Flare Rate (ct/s)', fontsize=fontsize)
plt.ylabel('Frequency', fontsize=fontsize)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("Ratio_Hist.pdf", dpi=300, bbox_inches='tight')
plt.close()

'''





mask0_duration = (fluence < 100) & (mean_rate > 0.03)
mask1_duration = (fluence > 100) & (fluence < 400) & (mean_rate > 0.03)
mask2_duration = (fluence > 400) & (mean_rate > 0.03)
mask3_duration = (fluence < 100) & (mean_rate < 0.03)
mask4_duration = (fluence > 100) & (fluence < 400) & (mean_rate < 0.03)
mask5_duration = (fluence > 400) & (mean_rate < 0.03)

plt.figure(figsize=(5, 4))
plt.errorbar(duration[mask0_duration], max_rate[mask0_duration], yerr=max_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(duration[mask1_duration], max_rate[mask1_duration], yerr=max_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(duration[mask2_duration], max_rate[mask2_duration], yerr=max_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(duration[mask3_duration], max_rate[mask3_duration], yerr=max_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(duration[mask4_duration], max_rate[mask4_duration], yerr=max_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='orange', markersize=markersize)
plt.errorbar(duration[mask5_duration], max_rate[mask5_duration], yerr=max_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='orange', markersize=markersize)
plt.axvline(6000, color='black', linestyle='--', linewidth=1)

plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Max Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Duration_MaxRate.pdf", dpi=300, bbox_inches='tight')
plt.close()






mask0_duration = (fluence < 100) & ~cut_mask
mask1_duration = (fluence > 100) & (fluence < 400) & ~cut_mask
mask2_duration = (fluence > 400) & ~cut_mask
mask3_duration = (fluence < 100) & cut_mask
mask4_duration = (fluence > 100) & (fluence < 400) & cut_mask
mask5_duration = (fluence > 400) & cut_mask

plt.figure(figsize=(5, 4))
plt.errorbar(mean_rate[mask0_duration], max_rate[mask0_duration], yerr=max_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(mean_rate[mask1_duration], max_rate[mask1_duration], yerr=max_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(mean_rate[mask2_duration], max_rate[mask2_duration], yerr=max_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(mean_rate[mask3_duration], max_rate[mask3_duration], yerr=max_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(mean_rate[mask4_duration], max_rate[mask4_duration], yerr=max_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='orange', markersize=markersize)
plt.errorbar(mean_rate[mask5_duration], max_rate[mask5_duration], yerr=max_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='orange', markersize=markersize)

plt.legend()
plt.xlabel('Mean Rate (ct/s)')
plt.ylabel('Max Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("MeanRate_MaxRate.pdf", dpi=300, bbox_inches='tight')
plt.close()





plt.figure(figsize=(5, 4))
plt.errorbar(fluence[mask0_duration], max_rate[mask0_duration], xerr=fluence_err[mask0_duration], yerr=max_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth, label='Weak Flares')
plt.errorbar(fluence[mask1_duration], max_rate[mask1_duration], xerr=fluence_err[mask1_duration], yerr=max_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, markeredgewidth=edgewidth, label='Moderate Flares')
plt.errorbar(fluence[mask2_duration], max_rate[mask2_duration], xerr=fluence_err[mask2_duration], yerr=max_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth, label='Strong Flares')
plt.errorbar(fluence[mask3_duration], max_rate[mask3_duration], xerr=fluence_err[mask3_duration], yerr=max_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(fluence[mask4_duration], max_rate[mask4_duration], xerr=fluence_err[mask4_duration], yerr=max_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='orange', markersize=markersize, markeredgewidth=edgewidth)
plt.errorbar(fluence[mask5_duration], max_rate[mask5_duration], xerr=fluence_err[mask5_duration], yerr=max_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, markeredgewidth=edgewidth)
plt.legend()
plt.xlabel('Fluence (ct)')
plt.ylabel('Max Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Fluence_MaxRate.pdf", dpi=300, bbox_inches='tight')
plt.close()






plt.figure(figsize=(5, 4))
plt.errorbar(duration[mask0_duration], mean_rate[mask0_duration], yerr=mean_rate_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(duration[mask1_duration], mean_rate[mask1_duration], yerr=mean_rate_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(duration[mask2_duration], mean_rate[mask2_duration], yerr=mean_rate_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(duration[mask3_duration], mean_rate[mask3_duration], yerr=mean_rate_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(duration[mask4_duration], mean_rate[mask4_duration], yerr=mean_rate_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize)
plt.errorbar(duration[mask5_duration], mean_rate[mask5_duration], yerr=mean_rate_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize)

plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Mean Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Duration_MeanRate.pdf", dpi=300, bbox_inches='tight')
plt.close()


plt.figure(figsize=(5, 4))
plt.errorbar(duration[mask0_duration], fluence[mask0_duration], yerr=fluence_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(duration[mask1_duration], fluence[mask1_duration], yerr=fluence_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(duration[mask2_duration], fluence[mask2_duration], yerr=fluence_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(duration[mask3_duration], fluence[mask3_duration], yerr=fluence_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(duration[mask4_duration], fluence[mask4_duration], yerr=fluence_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize)
plt.errorbar(duration[mask5_duration], fluence[mask5_duration], yerr=fluence_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize)
plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Fluence (ct)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Duration_Fluence.pdf", dpi=300, bbox_inches='tight')
plt.close()




plt.figure(figsize=(5, 4))
plt.errorbar(mean_rate[mask0_duration], fluence[mask0_duration], yerr=fluence_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(mean_rate[mask1_duration], fluence[mask1_duration], yerr=fluence_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(mean_rate[mask2_duration], fluence[mask2_duration], yerr=fluence_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(mean_rate[mask3_duration], fluence[mask3_duration], yerr=fluence_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(mean_rate[mask4_duration], fluence[mask4_duration], yerr=fluence_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize)
plt.errorbar(mean_rate[mask5_duration], fluence[mask5_duration], yerr=fluence_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize)
plt.legend()
plt.xlabel('Mean Rate (ct/s)')
plt.ylabel('Fluence (ct)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("MeanRate_Fluence.pdf", dpi=300, bbox_inches='tight')
plt.close()


plt.figure(figsize=(5, 4))
singlepeak = (shape == 1) & (mean_rate > 0.03)
doublepeak = (shape == 2) & (mean_rate > 0.03)
higherpeak = (shape == 3) & (mean_rate > 0.03)
mask_singlepeak = (shape == 1) & (mean_rate < 0.03)
mask_doublepeak = (shape == 2) & (mean_rate < 0.03)
mask_higherpeak = (shape == 3) & (mean_rate < 0.03)
plt.errorbar(duration[singlepeak], fluence[singlepeak], yerr=fluence_err[singlepeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Single Peaked Flares')
plt.errorbar(duration[doublepeak], fluence[doublepeak], yerr=fluence_err[doublepeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='green', markersize=markersize, label='Double Peaked Flares')
plt.errorbar(duration[higherpeak], fluence[higherpeak], yerr=fluence_err[higherpeak], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Complex Substructure')
plt.errorbar(duration[mask_singlepeak], fluence[mask_singlepeak], yerr=fluence_err[mask_singlepeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(duration[mask_doublepeak], fluence[mask_doublepeak], yerr=fluence_err[mask_doublepeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='green', markersize=markersize)
plt.errorbar(duration[mask_higherpeak], fluence[mask_higherpeak], yerr=fluence_err[mask_higherpeak], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize)

plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Fluence (ct)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Duration_Fluence_Shapes.pdf", dpi=300, bbox_inches='tight')
plt.close()


'''
plt.figure(figsize=(5, 4))
plt.errorbar(duration, median_rate, yerr=median_rate_err, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(duration[mask], median_rate[mask], yerr=median_rate_err[mask], ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(duration[mask2], median_rate[mask2], yerr=median_rate_err[mask2], ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.legend()
plt.xlabel('Duration (s)')
plt.ylabel('Median Rate (ct/s)')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Duration_MedianRate.pdf", dpi=300, bbox_inches='tight')
plt.close()
'''



hist, bin_edges = np.histogram(hardnessratio[hardnessratio < 3], bins=25)

import numpy as np
from scipy.optimize import curve_fit

def gaussian(x, amp, mean, sigma, offset):
    """
    1D Gaussian with constant offset.

    f(x) = amp * exp(-(x - mean)^2 / (2 * sigma^2)) + offset
    """
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2) + offset


def fit_gaussian(x, y, p0=None):
    """
    Fit a Gaussian to data (x, y).

    Parameters
    ----------
    x, y : array-like
        Data points.
    p0 : tuple or list, optional
        Initial guess (amp, mean, sigma, offset).
        If None, a reasonable guess is estimated from the data.

    Returns
    -------
    popt : ndarray
        Best-fit parameters [amp, mean, sigma, offset].
    perr : ndarray
        1σ uncertainties on the parameters.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if p0 is None:
        amp_guess = y.max() - y.min()
        mean_guess = np.sum(x * y) / np.sum(y)
        sigma_guess = 0.25 * (x.max() - x.min())
        offset_guess = y.min()
        p0 = (amp_guess, mean_guess, sigma_guess, offset_guess)

    popt, pcov = curve_fit(gaussian, x, y, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr

# Fit



mask0_duration = (fluence < 100) #& (mean_rate > 0.03)
mask1_duration = (fluence > 100) & (fluence < 400) #& (mean_rate > 0.03)
mask2_duration = (fluence > 400) #& (mean_rate > 0.03)
mask3_duration = (fluence < 100) #& (mean_rate < 0.03)
mask4_duration = (fluence > 100) & (fluence < 400) #& (mean_rate < 0.03)
mask5_duration = (fluence > 400) #& (mean_rate < 0.03)

fig, axs = plt.subplots(4, 1, figsize=(5, 10), sharex=True)
axs[0].hist(hardnessratio[mask0_duration], edgecolor='black', bins=np.linspace(0.5, 3, 10), color='steelblue', label='Weak Flares')
axs[0].legend()
axs[1].hist(hardnessratio[mask1_duration], edgecolor='black', bins=np.linspace(0.5, 3, 10), color='orange', label='Moderate Flares')
axs[1].legend()
axs[2].hist(hardnessratio[mask2_duration], edgecolor='black', bins=np.linspace(0.5, 3, 10), color='red', label='Strong Flares')
axs[2].legend()

axs[3].hist([hardnessratio[mask0_duration], hardnessratio[mask1_duration], hardnessratio[mask2_duration]],
    bins=np.linspace(0.5, 3, 10),
    edgecolor='black',
    color=['steelblue', 'orange', 'red'], 
    stacked = True,
    label=['Weak Flares', 'Moderate Flares', 'Strong Flares']
)
axs[3].legend()

hist, bin_edges = np.histogram(hardnessratio, bins=np.linspace(0.5, 3, 10))
popt, perr = fit_gaussian(bin_edges[:-1], hist)
amp, mean, sigma, offset = popt

axs[3].plot(bin_edges, gaussian(bin_edges, amp, mean, sigma, offset))

plt.subplots_adjust(hspace=0)
plt.xlabel('Hardness Ratio')
plt.savefig("HardnessRatioHist_Multiplot.pdf", dpi=300, bbox_inches='tight')



mask0_duration = (fluence < 100) & ~cut_mask
mask1_duration = (fluence > 100) & (fluence < 400) & ~cut_mask
mask2_duration = (fluence > 400) & ~cut_mask
mask3_duration = (fluence < 100) & cut_mask
mask4_duration = (fluence > 100) & (fluence < 400) & cut_mask
mask5_duration = (fluence > 400) & cut_mask

plt.figure(figsize=(5, 4))
plt.errorbar(fluence[mask0_duration], hardnessratio[mask0_duration], yerr=hardnessratio_err[mask0_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize, label='Weak Flares')
plt.errorbar(fluence[mask1_duration], hardnessratio[mask1_duration], yerr=hardnessratio_err[mask1_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize, label='Moderate Flares')
plt.errorbar(fluence[mask2_duration], hardnessratio[mask2_duration], yerr=hardnessratio_err[mask2_duration], ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize, label='Strong Flares')
plt.errorbar(fluence[mask3_duration], hardnessratio[mask3_duration], yerr=hardnessratio_err[mask3_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='steelblue', markersize=markersize)
plt.errorbar(fluence[mask4_duration], hardnessratio[mask4_duration], yerr=hardnessratio_err[mask4_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='darkorange', markersize=markersize)
plt.errorbar(fluence[mask5_duration], hardnessratio[mask5_duration], yerr=hardnessratio_err[mask5_duration], marker='X', ecolor='black', elinewidth=elinewidth, capsize=capsize, capthick=capthick, fmt='o', mec='black', mfc='red', markersize=markersize)
plt.legend()
plt.xlabel('Fluence (ct)')
plt.ylabel('Hardness Ratio')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.savefig("Fluence_HardnessRatio.pdf", dpi=300, bbox_inches='tight')
plt.close()




def plot_dNdF(x, nbins=10, color='k', label=None):
    bins = np.logspace(np.log10(min(x)), np.log10(max(x)), nbins)
    N, edges = np.histogram(x, bins=bins)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths  = edges[1:] - edges[:-1]
    dNdx = N / widths
    err = np.sqrt(N) / widths

    # 6. step histogram
    plt.step(edges[:-1], dNdx, where='post', color=color, lw=2, label=label)
    plt.errorbar(centers, dNdx, yerr=err, fmt='none', ecolor=color, capsize=3, lw=1)

    return centers, dNdx, err, edges

def model(x, N, g, c):
    return N * x**g * np.exp(-x / c)

def model2(x, N, g):
    return N * x**g

def reduced_chi2(x, y, sigma, model, params, k):
    y_model = model(x, *params)
    chi2 = np.sum(((y - y_model) / sigma)**2)
    nu = len(y) - k
    return chi2, chi2/nu

def fitting(x, label, label_letter):
    plt.figure(figsize=(6,4))

    centers, dNdx, err, edges = plot_dNdF(x, nbins=12, color='black', label='All Flares')

    p0 = [14, -1, 160]
    bounds = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])

    err[err == 0] = 1.0
    popt, pcov = curve_fit(model, centers, dNdx, p0=p0, bounds=bounds, sigma=err, absolute_sigma=True, maxfev=20000)  # uncomment if you have y_err)

    chi2, chi2_red = reduced_chi2(centers, dNdx, err, model, popt, k=3)
    #print("chi^2 =", chi2)
    #print("chi^2 / nu =", chi2_red)

    y_model = model(centers, *popt)

    plt.step(edges[:-1], y_model, where='post', color='red', lw=2, label='Model fit')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(label)
    plt.ylabel(f"dN/d{label_letter}")
    plt.legend()
    plt.title(f'$\chi^2$ {chi2:.3f}, Reduced $\chi^2$ {chi2_red:.3f}')
    plt.savefig(f'Param_Dist_{label}', dpi=300, bbox_inches='tight')

fitting(fluence, 'Fluence', 'F')
fitting(max_rate, 'Max Rate', 'R')
fitting(duration, 'Duration', 'T')
fitting(luminosity, 'Luminosity', 'L')
fitting(hardnessratio, 'Hardness Ratio', 'H')





pl = [2.07, 2.20, 2.3, 2.62]
pl_err = [0.11, 0.14, 0.13, 0.14]
avg_fluence = [1965, 615.4, 291.5, 142.621]
avg_fluence_err = [43.36, 24.63, 16.9, 11.96]
avg_lumin = np.array([28.74, 9.04, 6.31, 3.61])# * 10**34
avg_lumin_err = np.array([0.044, 0.063, 0.0366, 0.0477])# * 10**34


plt.figure(figsize=(5, 4))
plt.errorbar(avg_lumin, pl, xerr=avg_lumin_err, yerr=pl_err, ecolor='black', elinewidth=1, capsize=3, capthick=1, fmt='-o', mec='black', mfc='steelblue', markersize=markersize)
plt.xlabel('Average Luminosity of Flares Used to Calculate Power Law')
plt.ylabel('Power Law Index')
plt.xscale('log')
plt.yscale('log')
plt.savefig("Fluence_PL.pdf", dpi=300, bbox_inches='tight')
plt.close()