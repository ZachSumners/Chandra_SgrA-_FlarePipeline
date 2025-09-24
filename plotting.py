import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.time import Time

df = pd.read_csv('flare_properties.csv')

#df_q = pd.read_csv('quiescent_rates.csv')

mean_rate = df['rate_mean']
mean_rate_err = df['rate_mean_err']
max_rate = df['rate_max']
max_rate_err = df['rate_max_err']
duration = df['duration_s']
fluence = df['fluence_ct']

q1, q2, q3 = np.quantile(duration, [0.25, 0.50, 0.75])
# or, equivalently:
# q1, q2, q3 = np.percentile(x, [25, 50, 75])

print("Q1 =", q1)   # 25th percentile
print("Q2 =", q2)   # 50th percentile (median)
print("Q3 =", q3)   # 75th percentile

#quiescent = df_q['quiescent_rate']
#quiescent_err = df_q['quiescent_rate_err']
#obs_date = df_q['obs_date']

# ── Read your CSV as before ─────────────────────────────────────────────
#df_q = pd.read_csv("quiescent_rates.csv")

# ── 1) Drop any rows where obs_date is null ─────────────────────────────
#df_q = df_q.dropna(subset=["obs_date"])

# ── 2) Let pandas parse the ISO-strings into datetimes ──────────────────
#    We specify the exact format so that any malformed string becomes NaT.
#df_q["dt"] = pd.to_datetime(
#    df_q["obs_date"],
#    format="%Y-%m-%dT%H:%M:%S",
#    errors="coerce"
#)

# ── 3) Drop any rows that failed parsing (dt == NaT) ────────────────────
#df_q = df_q.dropna(subset=["dt"])

# ── 4) Convert that datetime64[ns] array into Astropy Time → get MJD ───
#    Astropy can directly consume a numpy array of datetime64[ns] with format="datetime64"
#t = Time(df_q["dt"].values, format="datetime64", scale="utc")
#df_q["quiescent_mjd"] = t.mjd

# ── 5) (Optional) If you still want to inspect the first few rows: ──────

#plt.errorbar(df_q['quiescent_mjd'], quiescent/1000, yerr=quiescent_err/1000, fmt='o', ecolor="black", alpha=0.8)
#plt.xlabel('MJD')
#plt.ylabel('Quiescent Count Rate (ct/s)')
#plt.ylim(0, 0.02)
#plt.grid()
#plt.title('Quiescent Rate over Time')
#plt.show()


print(df['duration_s'], df['rate_max'], df['obs_id'])

plt.errorbar(duration, max_rate, yerr=max_rate_err, fmt='o', ecolor="black", alpha=0.8)
plt.axhline(y=np.mean(max_rate), color="black", linestyle="--", linewidth=1.0, label="y = 0.5")
#plt.text(500, 0.12, f"{np.mean(max_rate):.2f} (ct/s)", ha="center", va="bottom", fontsize=8, color="black")
plt.yscale('log')
plt.xscale('log')
plt.title('Flare Duration and Max Rate')
plt.xlabel('Duration (s)')
plt.ylabel('Max Count Rate (ct/s)')
plt.grid()
plt.show()


plt.errorbar(duration, fluence, yerr=(duration*mean_rate_err), fmt='o', ecolor='black', alpha=0.8)
plt.yscale('log')
plt.xscale('log')
plt.title('Flare Duration and Fluence')
plt.xlabel('Duration (s)')
plt.ylabel('Fluence (ct)')
plt.grid()
plt.show()



plt.errorbar(fluence, max_rate, yerr=max_rate_err, xerr=(duration*mean_rate_err), fmt='o', ecolor='black', alpha=0.8)
plt.yscale('log')
plt.xscale('log')
plt.title('Flare Fluence and Max Rate')
plt.xlabel('Fluence (ct)')
plt.ylabel('Max Count Rate (ct/s)')
plt.grid()
plt.show()

#plt.hist(mean_rate, bins=50, edgecolor="black")           # 30 equally-spaced bins
#plt.xlabel("Mean Count Rate (ct/s)")
#plt.title("Mean Flare Strength")
#plt.tight_layout()
#plt.show()

mask = fluence > 200
candidates = df[['obs_id', 'fluence_ct']][mask]

print(candidates, len(candidates))