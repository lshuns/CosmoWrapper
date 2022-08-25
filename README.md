# CosmoWrapper

Scripts to wrap CosmoPipe to calibrate the redshifts and compute cosmological
constraints for the KiDS-1000 cosmic shear data:
- Define fiducial gold samples.
- Calculate the Nz.
- Run CosmoPipe:
  - Compute the COSEBI data vector.
  - Compute the covariance matrix.
  - Run MCMC chains.
  - Run least-squares minimiser (optional).
- Collect output data.

## Requirements

You need a local copy of the private repository `CosmoFisherForecast`, please
contact Benjamin Joachimi.

Optionally, you access to the non-public BitBucket
repository `marika_asgari/KiDS1000_chains` to run the least-squares minimiser,
please contact Marika Asgari. The wrapper script installs the code
automatically, which requires to configure an SSH access tokens for
`bitbucket.com`.

## Installation

Create conda environment for python 3.8 and R (recommended):

> conda env create --file environment.yml
> conda activate cosmowrapper
> bash INSTALL.sh

Create working directory at `PATH`:

> bash setup_workspace.sh [PATH]

This will link all scripts from `src` into `PATH` and copy the
`cosmowrapper.param` file.

## Configuration

Modify the copied `cosmowrapper.param` file in the working directory to your
needs. Then copy the required input files into the `$INPUTDIR`:
- `$SPECCAT_ALL`: Spectroscopic calibration data (van den Busch et al. 2022) as
  CSV file.
- `$PHOTCAT_ALL`: KiDS-1000 cosmic shear data, can be in FITS or LDAC format.
- `$MASKFILE`: KiDS-1000 mask file.
- `$SOMCOVFILE`: Covariance matrix of the delta-z nuisance parameters.

Additionally, some external code is required, specified through their
installation path:
- `$DIR_LDAC`: Folder containing the compiled LDAC tools.
- `$COSMOFISHER`: Path to CosmoFisherForecast root directory.

## Execution

Run the main wrapper script `cosmowrapper.sh`. The code is split into several
steps, which can be run separately (recommended). These steps are numbered:
 1) Preparing the input data.
 2) Reduce the size of the shear catalogue.
 3) Download SOM code and train the SOM.
 4) Define the Gold sample and compute the n(z).
 5) Combine the SOM columns into a single shear catalogue.
 6) Split catalogue in KiDS northern and southern field.
 7) Install CosmoPipe and its dependencies.
 8) Configure the CosmoPipe for the run.
 9) Collect and prepare the catalogues for CosmoPipe.
10) Compute data vector, covariance matrix and MCMC sampler.
11) Install and run minimiser (optional, requires access to `KiDS1000_chains`).
12) Download Planck chain and collect MCMC and minimiser output in `ChainDir`.
    These files are in `MultiNest` format and are named:  
    `GoldSet_*_multinest_*.txt` (MCMC)  
    `GoldSet_*_maxlike_*.txt` (minimiser)
