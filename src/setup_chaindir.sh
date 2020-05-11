
#Remove any previous chain directories
if [ -d ChainDir ]
then 
  rm -fr ChainDir
fi 
#Construct the chain directory
mkdir ChainDir
cd ChainDir
#Get the planck contour data
wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip -O base-plikHM-TTTEEE-lowl-lowE_R3.00.zip 
unzip base-plikHM-TTTEEE-lowl-lowE_R3.00.zip
cat base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE_?.txt > base_plikHM_TTTEEE_lowl_lowE.txt
ln -s base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE.paramnames . 
#Link the completed goldsets
ln -sf ../COSMOPIPE/GoldSet_NONE_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.txt GoldSet_NONE_shift.txt 
ln -sf ../COSMOPIPE/GoldSet_NONE_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.paramnames GoldSet_NONE_shift.paramnames 
ln -sf ../COSMOPIPE/GoldSet_NONE/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.txt GoldSet_NONE.txt  
ln -sf ../COSMOPIPE/GoldSet_NONE/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.paramnames GoldSet_NONE.paramnames  
ln -sf ../COSMOPIPE/GoldSet_noDEEP2_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.txt GoldSet_noDEEP2.txt 
ln -sf ../COSMOPIPE/GoldSet_noDEEP2_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.paramnames GoldSet_noDEEP2.paramnames 
ln -sf ../COSMOPIPE/GoldSet_Fid_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.txt GoldSet_Fiducial.txt 
ln -sf ../COSMOPIPE/GoldSet_Fid_shift/work_KV450/MCMC/output/sci_UNBLINDED/KV450__HEADER.paramnames GoldSet_Fiducial.paramnames 
