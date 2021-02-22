#Remove any previous chain directories
if [ -d ChainDir ]
then 
  rm -fr ChainDir
fi 
#Construct the chain directory
mkdir ChainDir
cd ChainDir
#Get the planck contour data
if [ ! -f base-plikHM-TTTEEE-lowl-lowE_R3.00.zip ]
then
  wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip -O base-plikHM-TTTEEE-lowl-lowE_R3.00.zip 
  unzip base-plikHM-TTTEEE-lowl-lowE_R3.00.zip
fi
cat base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE_?.txt > base_plikHM_TTTEEE_lowl_lowE.txt
ln -s base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE.paramnames . 
#Link the completed goldsets
if [ "$BLINDS" == "NONE" ]
then
  BLINDS="UNBLINDED"
fi
for GoldSet in ${GOLDLIST}
do
  ln -sf ../COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/MCMC/output/sci_${BLINDS}/${SURVEY}__HEADER.txt GoldSet_${GoldSet}.txt  
  ln -sf ../COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/MCMC/output/sci_${BLINDS}/${SURVEY}__HEADER.paramnames GoldSet_${GoldSet}.paramnames  
done
