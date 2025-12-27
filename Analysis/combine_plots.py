import os
import glob
from PIL import Image
import pandas as pd

# Load flare IDs (as strings to match filenames)
df = pd.read_csv('/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/flare_properties.csv')
flare_ids = set(df['obs_id'].astype(str).values)
print(len(flare_ids))

# Step 1: Find all PNG files in standard folders
png_paths = sorted(glob.glob("/Users/zachsumners/Desktop/Research/Chandra/Chandra_SgrA-_FlarePipelineNoGithub/*/repro/*lc_2-8keV.png"))
#print(png_paths)
# Step 2: Include special 1561 subfolders
#png_paths += [
#    '1561/1561_000/repro/01561_300lc_2-8keV.png'
#]

# Step 3: Filter only those where the filename starts with a known flare ID
def has_flare_id(path):
    base = os.path.basename(path)
    obsid = base.split('_')[0].lstrip('0')  # remove leading zeros for match
    return obsid in flare_ids

filtered_paths = [p for p in png_paths if has_flare_id(p)]
print(filtered_paths)
# Step 4: Open images and save to PDF
images = [Image.open(p).convert('RGB') for p in filtered_paths]

if images:
    images[0].save("only_flare_plots.pdf", save_all=True, append_images=images[1:])
    print(f"Saved {len(images)} flare plots to all_flare_plots.pdf")
else:
    print("No matching flare plots found.")