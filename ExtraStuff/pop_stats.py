import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.optimize import curve_fit

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


print(df['hardness_ratio'].describe())
print(df['hardness_ratio'].median())