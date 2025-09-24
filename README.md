# Chandra_SgrA-_FlarePipeline

Chandra_SgrA-_FlarePipeline is a Python package for identifying flares in Chandra data -- primarily targeting SgrA* -- using a Bayesian Blocks algorithm.

## Installation

You can make a local copy of the code by cloning this repository. The 'main' branch contains the most up to date version of the code.
If you would like to run a version of the code where no user inputs to the terminal are required, you can use the 'no-intervention' branch.

You will also need to install CIAO and obtain the associated CALDB files. Instructions for installing the latest release of CIAO can be found 
[here](https://cxc.cfa.harvard.edu/ciao/threads/ciao_install_tool/). Optionally, if you would like to enable remote CALDB access (so you don't need to always re-download the latest files), follow the instructions [here](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_remote_access.html).

Advanced pileup correction requires the installation of MARX found [here](https://space.mit.edu/asc/marx/index.html).

## Setup

Before running the pipeline, you will need to go to the Chandra Data Archive [CDA](https://cxc.harvard.edu/cda/), download your observation data folder(s) and put them in the **same directory** as the flare pipeline scripts. 

```ChandraLightCurvePipeline.py``` is the main script that will execute the rest of the pipeline. You can run it by entering ```python3 {path/to/Chandra_SgrA-_FlarePipeline}/ChandraLightCurvePipeline.py``` in the terminal.

## User Inputs

The ```ChandraLightCurvePipeline.py``` user inputs include: 
* observationID: (string) Sets the Chandra-assigned observation ID of the relevant data folder
* barycentric: (Boolean) Sets whether to perform a barycenter correction to the Chandra observation
* wcsCorrect: (Boolean) Sets whether to perform a WCS correction to improve the astrometry of the observation
* reprocess: (Boolean) Sets whether to reprocess the data using the latest Chandra calibration files
* search: (Boolean) Sets whether to manually identify all point sources in the image, allowing the user to select the source of interest (e.g., SgrA*)
* erange: (list) Sets the lower and upper energy range bounds for lightcurve extraction
* tbin: (integer) Sets the time bin size in seconds for lightcurve extraction
* src_coords: (list) Sets the coordinates in degrees of the source of interest (e.g., SgrA*) to use in the source region selection
* bkg_coords: (list) Sets the size of the inner and outer radii of the background selection region in pixels

## Outputs

The pipeline will output a 'repro' sub-folder into the folder containing the observeration data. 'repro'contains various CIAO output files 
as well as lightcurve files. The pipeline runs a Bayesian Blocks algorithm to detect flares in the light curve, and outputs the results of the flare detection to a 'Results' sub-folder.

## Additional Notes

For more information on the data calibration methods and flare detection algorithm, see [this paper](https://ui.adsabs.harvard.edu/abs/2019ApJ...884..148B/abstract) from Elie Bouffard,
and the associated code can be found [here](https://github.com/Elie23/X-ray-flare-simulator/tree/master).
