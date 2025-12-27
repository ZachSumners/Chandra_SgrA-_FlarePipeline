import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('flare_properties.csv')
mean_rate = df['rate_mean']
mean_rate_err = df['rate_mean_err']
max_rate = df['rate_max']
max_rate_err = df['rate_max_err']
duration = df['duration_s']
fluence = df['fluence_ct']
flare_ids = df['obs_id']
flare_numbers = df['flare_number']

df_neilsen = pd.read_csv('NeilsenTable.csv')

neilsen_flares = df_neilsen[df_neilsen["Flare"].notna()]
neilsen_duration = neilsen_flares['Duration']
neilsen_max_rate = neilsen_flares['Peak Rate']

neilsen_max_rate_upper_err = neilsen_flares['Peak Rate Upper Err']
neilsen_max_rate_lower_err = neilsen_flares['Peak Rate Lower Err']


mask = (df['obs_id'].isin(df_neilsen['ObsID'].unique()))

def classify_flare(obs_id, flare_number):
    key = (obs_id, flare_number)

    single_peak = {
        (10556, 1), (10556, 3), (10556, 4), (13838, 1), (13839, 1),
        (13842, 2), (13854, 1), (13854, 2), (13854, 3), (13854, 4), (13852, 1), (13852, 2),
        (13857, 1), (14463, 1), (14466, 1), (15041, 2), (15042, 1), (15045, 2),
        (16963, 1), (16508, 1), (18731, 1), (20446, 1), (21454, 1), (22230, 1), (23666, 1),
        (23665, 1), (23739, 1), (28231, 1), (3393, 3), (4684, 1), (6363, 1), (1561, 1)
    }

    double_peak = {
        (10556, 2), (13843, 1), (13847, 1), (13851, 1), (13849, 1), (14462, 2), (14468, 1),
        (14944, 1), (15043, 1), (15042, 2), (15570, 1), (16966, 1), (18732, 1), (20041, 1),
        (20346, 1), (21454, 2), (21456, 1), (25297, 1), (26760, 1), (3393, 1), (3393, 2),
        (5953, 1), (22595, 2), (22595, 1), (14392, 1), (13839, 2)
    }

    higher_substructure = {
        (13842, 1), (14432, 2), (14392, 2), (16218, 1), (20751, 1), (22594, 1), (22595, 3), (22592, 1), (1561, 2), (11843, 1)
    }

    ambiguous = {
        (13016, 1), (13840, 1), (13842, 3), (13845, 1), (14427, 1), (14427, 2),
        (14462, 1), (14439, 1), (14465, 1), (14465, 2), (14468, 2), (15041, 1), (15045, 1),
        (16217, 1), (18055, 1), (19703, 1), (19703, 2),
        (22937, 1), (26759, 1), (3392, 1), (3392, 2), (3392, 3),
        (3663, 1), (5952, 1), (9173, 1), (18055, 2), (14432, 1)
    }
    
    if key in single_peak:
        return "single_peak"
    elif key in double_peak:
        return "double_peak"
    elif key in higher_substructure:
        return "higher_substructure"
    elif key in ambiguous:
        return "ambiguous"
    else:
        return "unclassified"

plotted = {
    'single_peak': False,
    'double_peak': False,
    'higher_substructure': False,
    'ambiguous': False,
    'unclassified': False
}

for i, flare in enumerate(flare_ids):
    flaretype = classify_flare(flare, flare_numbers[i])

    label = None
    if not plotted[flaretype]:
        label = flaretype.replace('_', ' ').title()
        plotted[flaretype] = True
    

    if flaretype=='single_peak':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='r', ecolor="black", alpha=0.8, label=label)
        #single.append(flare)
        plt.errorbar(duration[i], max_rate[i], yerr=max_rate_err[i], fmt='o', ecolor='black', c='r', alpha=0.8, label=label)
    elif flaretype=='double_peak':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='b', ecolor="black", alpha=0.8, label=label)
        #double.append(flare)
        plt.errorbar(duration[i], max_rate[i], yerr=max_rate_err[i], fmt='o', ecolor='black', c='b', alpha=0.8, label=label)
    elif flaretype=='higher_substructure':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='g', ecolor="black", alpha=0.8, label=label)
        #higher.append(flare)
        plt.errorbar(duration[i], max_rate[i], yerr=max_rate_err[i], fmt='o', ecolor='black', c='g', alpha=0.8, label=label)
    #elif flaretype=='ambiguous':
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='orange', ecolor="black", alpha=0.8, label=label)
        #ambiguous.append(flare)
        #plt.errorbar(duration[i], fluence[i], yerr=(duration*mean_rate_err), s=12, c='orange', alpha=0.8, label=label)
    #else:
        #plt.errorbar(duration[i], max_rate[i], yerr=mean_rate_err[i], fmt='o', c='gray', ecolor="black", alpha=0.8, label=str(flare))
        #unclassified.append(flare)
        #plt.errorbar(duration[i], fluence[i], yerr=(duration*mean_rate_err), s=12, c='gray', alpha=0.8, label=label)
   

#plt.errorbar(, max_rate[i], yerr=mean_rate_err[i], fmt='o', ecolor="black", alpha=0.8, label='Sumners Ford XVP')
    #plt.errorbar(neilsen_duration, neilsen_max_rate, yerr=[np.abs(neilsen_max_rate_lower_err), np.abs(neilsen_max_rate_lower_err)], fmt='o', marker='v', ecolor='black', alpha=0.8, label='Neilsen et al. (2013)')

plt.legend()
#plt.axhline(y=np.mean(mean_rate), color="black", linestyle="--", linewidth=1.0, label="y = 0.5")
#plt.text(500, 0.12, f"{np.mean(mean_rate):.2f} (ct/s)", ha="center", va="bottom", fontsize=8, color="black")
plt.yscale('log')
plt.xscale('log')
plt.title('Flare Duration and Max Rate')
plt.xlabel('Duration (s)')
plt.ylabel('Max Flare Rate (ct/s)')
plt.grid()
plt.show()
