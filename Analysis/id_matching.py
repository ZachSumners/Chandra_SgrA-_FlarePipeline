import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

df = pd.read_csv('exposure_date_and_flare_counts.csv')

exposure = df['exposure_ks']
flares = df['num_flares']

# ── 1) Drop any rows where obs_date is null ─────────────────────────────
df = df.dropna(subset=["obs_date"])

# ── 2) Let pandas parse the ISO-strings into datetimes ──────────────────
#    We specify the exact format so that any malformed string becomes NaT.
df["dt"] = pd.to_datetime(
    df["obs_date"],
    format="%Y-%m-%dT%H:%M:%S",
    errors="coerce"
)

# ── 3) Drop any rows that failed parsing (dt == NaT) ────────────────────
df = df.dropna(subset=["dt"])

# ── 4) Convert that datetime64[ns] array into Astropy Time → get MJD ───
#    Astropy can directly consume a numpy array of datetime64[ns] with format="datetime64"
t = Time(df["dt"].values, format="datetime64", scale="utc")
df["quiescent_mjd"] = t.mjd

januaryonemjd = np.array([
    49718.5, 50083.5, 50449.5, 50814.5, 51179.5,
    51544.5, 51910.5, 52275.5, 52640.5, 53005.5,
    53371.5, 53736.5, 54101.5, 54466.5, 54832.5,
    55197.5, 55562.5, 55927.5, 56293.5, 56658.5,
    57023.5, 57388.5, 57754.5, 58119.5, 58484.5,
    58849.5, 59215.5, 59580.5, 59945.5, 60310.5,
    60676.5, 61041.5
])
year = np.arange(1994, 2026, 1)
flare_year_rate = np.zeros(len(januaryonemjd))

for i in range(len(januaryonemjd)-1):
    mask = (df['quiescent_mjd'] > januaryonemjd[i]) & (df['quiescent_mjd'] < januaryonemjd[i+1]) # produces a boolean Series: [False, False, True, True]
    exposure_lengths = df['exposure_ks'][mask]
    flares = df['num_flares'][mask]

    total_flares = np.sum(flares.values)
    total_time = np.sum(exposure_lengths.values)

    if total_time != 0:
        flare_year_rate[i] = total_flares/total_time
    else:
        flare_year_rate[i] = np.nan
    print(total_time, total_flares, flare_year_rate[i]*86.4, year[i])


flare_year_rate_nonans = flare_year_rate[~np.isnan(flare_year_rate)]
print(np.nanmean(flare_year_rate))

q1, q2, q3 = np.quantile(flare_year_rate_nonans, [0.25, 0.50, 0.75])
# or, equivalently:
# q1, q2, q3 = np.percentile(x, [25, 50, 75])

print("Q1 =", q1*86.4)   # 25th percentile
print("Q2 =", q2*86.4)   # 50th percentile (median)
print("Q3 =", q3*86.4)   # 75th percentile

# Interquartile range
iqr = q3 - q1
print("IQR =", iqr*86.4)

year = np.arange(1994, 2026, 1)
print(year.shape, flare_year_rate.shape)
plt.scatter(year[:-1]+1, flare_year_rate[:-1]*86.4, s=10)
plt.xlabel('Year')
plt.ylabel('Flares/day of observing')
plt.title('Flaring Rate over Chandra Operation')
plt.grid()
#plt.ylim(0, 5)
plt.show()