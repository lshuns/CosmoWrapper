# CosmoWrapper

Scripts to wrap the CosmoPipe for use in Wright et al 2020b


Script takes the KiDS+VIKING-450 dataset and computes: 
- 6 sets of gold samples:
  -> Fiducial gold class, restricting analysis to photometric sources which are represented by any specz 
  -> 3x gold classes computed without various spectroscopic surveys (DEEP2, zCOSMOS, VVDS)
  -> 1x redshift quality gold class (calibrating only with 'certain' confidence redshifts)
  -> 1x survey gold class, restricting analysis to sources represented by spectra from at least 3 specz samples
- Calculates the Nz for each of these samples
- Runs the CosmoPipe
- Plots the results 

## INSTALL

Create conda environment:
> conda env create --file environment.yml
> conda activate cosmowrapper
> bash INSTALL.sh

Create working directory at `PATH`:
> bash setup_workspace.sh [PATH]
This will link all scripts from `src` into `PATH` and copy the
`vandenBusch21.param` file.
