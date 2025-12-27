import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.time import Time
from datetime import datetime

df = pd.read_csv('quiescent_rates.csv')
df_flares = pd.read_csv('flare_properties.csv')

dates = df['obs_date']
exposure = df['exposure_ks']
duration = df_flares['duration_s']
flare_start = df_flares['start_mjd']

#[1998, 1999, ... , 2024, 2025]
total_exposure = np.zeros(28)
flare_exposure = np.zeros(28)
num_flares = np.zeros(28)

for i, date in enumerate(flare_start):
    t = Time(float(date), format='mjd')
    decimal_year = t.decimalyear 
    flare_exposure[int(decimal_year) - 1998] += duration[i]
    num_flares[int(decimal_year) - 1998] += 1

for j, date in enumerate(dates):
    dt = datetime.fromisoformat(date)

    # Compute the fraction of the year elapsed
    start_of_year = datetime(dt.year, 1, 1)
    start_next_year = datetime(dt.year + 1, 1, 1)
    year_length = (start_next_year - start_of_year).total_seconds()
    elapsed = (dt - start_of_year).total_seconds()

    decimal_year = dt.year + elapsed / year_length

    total_exposure[int(decimal_year) - 1998] += exposure[j]*1000

print(total_exposure, flare_exposure, num_flares)

days_obs = total_exposure/86400
print(days_obs)
rate = num_flares/days_obs

years = np.arange(1998, 2026)

plt.scatter(years, rate)
plt.plot(years, days_obs, label='Days Observed')
plt.grid()
plt.ylabel('Flares/day')
plt.xlabel('Year')
plt.legend()
plt.show()
    