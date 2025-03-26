# Chandra_SgrA-_FlarePipeline

Chandra_SgrA-_FlarePipeline is a Python package for identifying flares in Chandra data -- primarily targeting SgrA* -- using a Bayesian Blocks algorithm.

## Installation

You can make a local copy of the code by cloning this repository. The 'main' branch contains the most up to date version of the code.
If you would like to run a version of the code where no user inputs to the terminal are required, you can use the 'no-intervention' branch.

## Setup

Before running the pipeline, you will need to go to the Chandra Data Archive [CDA](https://cxc.harvard.edu/cda/) 
and download your observation data folder(s) and put them in the **same directory** as the flare pipeline scripts. 

```ChandraLightCurvePipeline.py``` is the main script that will execute the rest of the pipeline. It takes several user inputs.

## User Inputs and Outputs

The user inputs include: 
* observationID: (string) Sets the Chandra-assigned observation ID of the relevant data folder
* barycentric: (Boolean) Sets whether to perform a barycenter correction to the Chandra observation
* 
