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

To run the code: 
  - Run the INSTALL.sh script to install the required R packages 
  - copy/link the R, shell, and python scripts into your desired running directory
  - Copy the specz compilation and (per-patch) KiDS shear catalogues into CosmoWrapper_Inputs/
  - Execute the Wright2019b.sh script 


