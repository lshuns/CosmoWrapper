#
#
# Script for running the fiducial results for the Wright 2020 (b) paper. 
# Uses KV450 Catalogues to run Cosmological Chains with different GoldClasses
#

set -e 

#Source the default parameters file /*fold*/ {{{
source Wright2020b.param
#/*fend*/}}}

#Define the root directory  /*fold*/ {{{
ROOT=`pwd`
#/*fend*/}}}

#Check for the input directory and files  /*fold*/ {{{
if [ -d ${INPUTDIR} ]
then 
  if [ ! -f ${INPUTDIR}/${SPECCAT_ALL} ]
  then 
    echo "ERROR: The Input Directory does not contain the spectroscopic catalogue!"
    exit 1
  elif [ "`ls ${INPUTDIR}/${PHOTCAT_ALL/ALL/*} | wc -l`" == "0" ]
  then 
    echo "ERROR: The Input Directory does not contain the photometric catalogue!"
    exit 1
  elif [ ! -f ${INPUTDIR}/${MASKFILE} ]
  then
    echo "ERROR: The Input Directory does not contain the survey footprint mask!"
    exit 1
  fi 
else 
  echo "ERROR: The Input Directory (containing the spectroscopic & photometric catalogues, and survey footprint mask) does not exist!"
  exit 1
fi 
#/*fend*/}}}

#Make the output directory  /*fold*/ {{{
mkdir -p ${OUTPUTDIR}
#/*fend*/}}}

#Spec cat adapt file name  /*fold*/ {{{
SPECCAT_ALL_ADAPT=${SPECCAT_ALL//.csv/_adapt.fits}
#/*fend*/}}}

if [ "${1}" == "1" -o "${1}" == "" ]
then
  #Construct the Spectrosopic Adapt Catalogue  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} ]
  then 
    echo "Constructing Spectroscopic Adapt Catalogue {"
    ${P_RSCRIPT} construct_adapt_catalogue.R ${INPUTDIR}/${SPECCAT_ALL} ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT}
    echo "} - Done" 
  else 
    echo "Spectroscopic Adapt Catalogue Already Exists! Skipping!"
  fi 
  #/*fend*/}}}
fi

#Calculate Patch Numbers  /*fold*/ {{{
PHOTCAT_PATCHES=""
NPATCH=0
for PATCH in ${PATCHLIST}
do 
  let NPATCH=$NPATCH+1
  #Combine the patch filenames 
  PHOTCAT_PATCHES="${PHOTCAT_PATCHES} ${INPUTDIR}/${PHOTCAT_ALL//${PATCHALL}/${PATCH}}"
done 
#/*fend*/}}}

if [ "${1}" == "2" -o "${1}" == "" ]
then
  #Construct the Combined Photometry Catalogue from Patches  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL} ]
  then 
    echo -e "Constructing Combined Photometry Catalogue from Patches\n"
    #Paste the patch catalogues together
    ${DIR_LDAC}/ldacpaste${THELI} -i ${PHOTCAT_PATCHES} -o ${OUTPUTDIR}/${PHOTCAT_ALL}
    echo "\n- Done"
  else
    echo "Combined Photometry Catalogue Already Exists! Skipping!"
  fi
  #/*fend*/}}}
fi

#Phot cat DIR columns 
PHOTCAT_ALL_DCOL=${PHOTCAT_ALL//.cat/_DIRcols.cat}

if [ "${1}" == "3" -o "${1}" == "" ]
then
  #Select DIRcol subset (within LDAC)  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} ]
  then 
    echo -e "Constructing the DIR Column Photometry Catalogue\n"
    list=`${DIR_LDAC}/ldacdesc${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL} -t OBJECTS | 
      grep "Key name" | awk -F. '{print $NF}' | 
      grep -v "MAG_GAAP_\|recal\|autocal\|MAG_AUTO\|SeqNr\|THELI_\|_B\|ID"`
    ${DIR_LDAC}/ldacdelkey${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL} -k ${list} -o ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL}
    echo "\n- Done"
  else
    echo "DIR Column Photometry Catalogue Already Exists! Skipping!"
  fi
  #/*fend*/}}}
fi

if [ "${1}" == "4" -o "${1}" == "" ]
then
  if [ ! -f SOM_DIR.R ]
  then
    #Install the SOM_DIR code {{{
    cd ${ROOT}
    mkdir -p INSTALL
    cd ${ROOT}/INSTALL
    if [ ! -d SOM_DIR ]
    then 
      git clone https://github.com/AngusWright/SOM_DIR.git
    fi 
    bash SOM_DIR/INSTALL.sh 
    ln -s SOM_DIR/R/SOM_DIR.sh . 
    #}}}
  fi 
  #Construct the Fiducial SOM  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${SOMFILE} ]
  then 
    echo "Constructing the Fiducial SOM" 
    ${P_RSCRIPT} SOM_DIR.R \
      -r ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} -t ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
      --toroidal --topo hexagonal --som.dim 101 101 -np -fn Inf \
      -sc 32 --only.som \
      -o ${OUTPUTDIR} -of ${SOMFILE//_SOMdata/} \
      --zr.label Z_B --zt.label z_spec_B \
      -k MAG_GAAP_u-MAG_GAAP_g \
      MAG_GAAP_u-MAG_GAAP_r MAG_GAAP_g-MAG_GAAP_r \
      MAG_GAAP_u-MAG_GAAP_i MAG_GAAP_g-MAG_GAAP_i \
      MAG_GAAP_r-MAG_GAAP_i MAG_GAAP_u-MAG_GAAP_Z \
      MAG_GAAP_g-MAG_GAAP_Z MAG_GAAP_r-MAG_GAAP_Z \
      MAG_GAAP_i-MAG_GAAP_Z MAG_GAAP_u-MAG_GAAP_Y \
      MAG_GAAP_g-MAG_GAAP_Y MAG_GAAP_r-MAG_GAAP_Y \
      MAG_GAAP_i-MAG_GAAP_Y MAG_GAAP_Z-MAG_GAAP_Y \
      MAG_GAAP_u-MAG_GAAP_J MAG_GAAP_g-MAG_GAAP_J \
      MAG_GAAP_r-MAG_GAAP_J MAG_GAAP_i-MAG_GAAP_J \
      MAG_GAAP_Z-MAG_GAAP_J MAG_GAAP_Y-MAG_GAAP_J \
      MAG_GAAP_u-MAG_GAAP_H MAG_GAAP_g-MAG_GAAP_H \
      MAG_GAAP_r-MAG_GAAP_H MAG_GAAP_i-MAG_GAAP_H \
      MAG_GAAP_Z-MAG_GAAP_H MAG_GAAP_Y-MAG_GAAP_H \
      MAG_GAAP_J-MAG_GAAP_H MAG_GAAP_u-MAG_GAAP_Ks \
      MAG_GAAP_g-MAG_GAAP_Ks MAG_GAAP_r-MAG_GAAP_Ks \
      MAG_GAAP_i-MAG_GAAP_Ks MAG_GAAP_Z-MAG_GAAP_Ks \
      MAG_GAAP_Y-MAG_GAAP_Ks MAG_GAAP_J-MAG_GAAP_Ks \
      MAG_GAAP_H-MAG_GAAP_Ks MAG_AUTO 
  else 
    echo "Fiducial SOM Already Exists! Skipping!" 
  fi 
  #/*fend*/}}}
fi

if [ "${1}" == "4a" -o "${1}" == "" ]
then
  #Construct the Tomographic Bin Catalogues /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT//_adapt/_adapt_TOMO1} ]
  then 
    echo "Generating Tomographic Bin Catalogues" 
    ${P_RSCRIPT} construct_tomo_bins.R ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL}
    ${P_RSCRIPT} construct_tomo_bins.R ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT}
  fi 
  #/*fend*/}}}
  #Optimise the cluster size for each tomographic bin /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${SOMFILE//_SOMdata/_TOMO1} ]
  then 
    echo "Optimising Cluster Sizes" 
    for TOMO in `seq 5`
    do 
      ${P_RSCRIPT} SOM_DIR.R \
        -r ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//_DIRcols.cat/_DIRcols_TOMO${TOMO}.fits} -t ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT//_adapt/_adapt_TOMO${TOMO}} \
        --toroidal --topo hexagonal --som.dim 101 101 -np -fn Inf \
        -sc 4 --optimise --refr.flag -cr recal_weight \
        --old.som ${OUTPUTDIR}/${SOMFILE} \
        -o ${OUTPUTDIR} -of ${SOMFILE//_SOMdata/_TOMO${TOMO}} \
        --zr.label Z_B --zt.label z_spec_B \
        -k MAG_GAAP_u-MAG_GAAP_g \
        MAG_GAAP_u-MAG_GAAP_r MAG_GAAP_g-MAG_GAAP_r \
        MAG_GAAP_u-MAG_GAAP_i MAG_GAAP_g-MAG_GAAP_i \
        MAG_GAAP_r-MAG_GAAP_i MAG_GAAP_u-MAG_GAAP_Z \
        MAG_GAAP_g-MAG_GAAP_Z MAG_GAAP_r-MAG_GAAP_Z \
        MAG_GAAP_i-MAG_GAAP_Z MAG_GAAP_u-MAG_GAAP_Y \
        MAG_GAAP_g-MAG_GAAP_Y MAG_GAAP_r-MAG_GAAP_Y \
        MAG_GAAP_i-MAG_GAAP_Y MAG_GAAP_Z-MAG_GAAP_Y \
        MAG_GAAP_u-MAG_GAAP_J MAG_GAAP_g-MAG_GAAP_J \
        MAG_GAAP_r-MAG_GAAP_J MAG_GAAP_i-MAG_GAAP_J \
        MAG_GAAP_Z-MAG_GAAP_J MAG_GAAP_Y-MAG_GAAP_J \
        MAG_GAAP_u-MAG_GAAP_H MAG_GAAP_g-MAG_GAAP_H \
        MAG_GAAP_r-MAG_GAAP_H MAG_GAAP_i-MAG_GAAP_H \
        MAG_GAAP_Z-MAG_GAAP_H MAG_GAAP_Y-MAG_GAAP_H \
        MAG_GAAP_J-MAG_GAAP_H MAG_GAAP_u-MAG_GAAP_Ks \
        MAG_GAAP_g-MAG_GAAP_Ks MAG_GAAP_r-MAG_GAAP_Ks \
        MAG_GAAP_i-MAG_GAAP_Ks MAG_GAAP_Z-MAG_GAAP_Ks \
        MAG_GAAP_Y-MAG_GAAP_Ks MAG_GAAP_J-MAG_GAAP_Ks \
        MAG_GAAP_H-MAG_GAAP_Ks MAG_AUTO >> ${OUTPUTDIR}/cluster_optimisation.log
    done
  else 
    echo "Cluster Optimisation already done! Skipping!" 
  fi 
  #/*fend*/}}}
fi

if [ "${1}" == "5" -o "${1}" == "" ]
then
  ##Construct the Gold Classes  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} ]
  then 
    echo "Constructing the Goldclass subsets" 
    time ${P_RSCRIPT} construct_dr4_goldclasses.R \
      -p ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} \
      -s ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
      --som ${OUTPUTDIR}/${SOMFILE} \
      --blinds ${BLINDS} \
      --outputpath ${OUTPUTDIR}/ \
      --nzformat .asc \
      >> ${OUTPUTDIR}/construct_dr4_goldclasses.log
  else 
    echo "Gold Catalogues Already Exists! Skipping!" 
  fi 
  #/*fend*/}}}
fi

PHOTCAT_ALL_GOLD=${PHOTCAT_ALL//.cat/_goldclasses.cat}
if [ "${1}" == "6" -o "${1}" == "" ]
then
  #Merge the new GoldClasses back with the original catalogue  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} ] 
  then 
    echo "Constructing merged catalogue" 
    GOLDFLAGLIST=""
    for goldset in ${GOLDLIST}
    do 
      if [ "`echo ${goldset} | grep -c NONE`" == "0" ] 
      then 
        GOLDFLAGLIST=`echo $GOLDFLAGLIST Flag_SOM_${goldset//_shift/}`
      fi 
    done
    python merge_ldac_and_fits.py ${OUTPUTDIR}/${PHOTCAT_ALL} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} "${BLINDS}" "${GOLDFLAGLIST}" 
  else 
    echo "Merged catalogue already exists! Skipping!" 
  fi 
  #/*fend*/}}}
fi

if [ "${1}" == "7" -o "${1}" == "" ]
then
  #Recreate all patchwise catalogues  /*fold*/ {{{
  for PATCHNUM in `seq ${NPATCH}`
  do 
    #Get the patch definitions 
    PATCH=`echo $PATCHLIST | awk -v n=${PATCHNUM} '{print $n}'`
    TMPLO=`echo $PATCHLO   | awk -v n=${PATCHNUM} '{print $n}'`
    TMPHI=`echo $PATCHHI   | awk -v n=${PATCHNUM} '{print $n}'`
    #Filter the ALLPATCH catalogue on RA
    if [ $TMPLO -lt $TMPHI ]
    then 
      #Normal RA limits: pick lo < RA <= hi 
      ${DIR_LDAC}/ldacfilter${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
    else 
      #RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
      ${DIR_LDAC}/ldacfilter${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})OR(ALPHA_J2000<=${TMPHI}));"
    fi 
  done
  #/*fend*/}}}
fi

if [ "${1}" == "8" -o "${1}" == "" ]
then
  #Run the CosmoPipe Installation  /*fold*/ {{{
  cd ${ROOT}
  mkdir -p INSTALL
  cd ${ROOT}/INSTALL
  if [ ! -d COSMOPIPE ]
  then 
    git clone https://github.com/AngusWright/COSMOPIPE.git
  else 
    cd ${ROOT}/INSTALL/COSMOPIPE
    git pull 
  fi 
  cd ${ROOT}
  
  bash INSTALL/CosmoPipe/COSMOLOGY_MASTER_INSTALL.sh \
    --noconfig \
    --packroot ${ROOT}/INSTALL/CosmoPipe/ \
    --runroot ${ROOT}/COSMOPIPE/ \
    --storagepath GoldSet_@@GOLDSET@@/${STORAGEDIR}/ \
    --runtime GoldSet_@@GOLDSET@@/RUNTIME/ \
    --configpath GoldSet_@@GOLDSET@@/RUNTIME/config/ \
    --scriptpath GoldSet_@@GOLDSET@@/RUNTIME/scripts/ \
    --nzfileid ${PHOTCAT_ALL_DCOL//.cat/_@@GOLDSET@@_blindNONE_TOMO} \
    --nzfilesuffix _Nz.asc \
    --shearsubset Flag_SOM_@@GOLDSET@@_NONE \
    --patchpath ${ROOT}/COSMOPIPE/GoldSet_@@GOLDSET@@/PatchData/ \
    --filesuffix _v2_good_goldclasses \
    --cosmopipecfname COSMOPIPE_CF_@@GOLDSET@@ 
  #/*fend*/}}}
fi

if [ "${1}" == "9" -o "${1}" == "" ]
then
  #Update the configure files  /*fold*/ {{{
  for GoldSet in ${GOLDLIST} 
  do
    cd ${ROOT}/COSMOPIPE
    #Setup the Gold Set configurations  /*fold*/ {{{
    sed "s/@@GOLDSET@@/${GoldSet}/g" configure.sh > configure_${GoldSet}.sh 
    #/*fend*/}}}
    
    GoldSetLink=${GoldSet//_shift/}
    GoldSetTail=${GoldSet//$GoldSetLink/}
    if [ "${GoldSetLink}" == "NONE" -o "${GoldSet}" == "NONE_catastrophic" ]
    then 
      #Remove the shear subsetting (I.e. no gold selection!) /*fold*/ {{{
      echo "Removing shear subset for ${GoldSet}"
      sed -i.bak "s/Flag_SOM_${GoldSet}_NONE/NONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSetTail}" == "_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_${GoldSet}_NONE/Flag_SOM_${GoldSetLink}_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/${GoldSet}_blindNONE/${GoldSetLink}_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "multispec3_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_multispec3_shift_NONE/Flag_SOM_multispec3_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/multispec3_shift_blindNONE/multispec3_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "noVVDS_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_noVVDS_shift_NONE/Flag_SOM_noVVDS_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/noVVDS_shift_blindNONE/noVVDS_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "nozCOSMOS_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_nozCOSMOS_shift_NONE/Flag_SOM_nozCOSMOS_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/nozCOSMOS_shift_blindNONE/nozCOSMOS_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "noDEEP2_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_noDEEP2_shift_NONE/Flag_SOM_noDEEP2_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/noDEEP2_shift_blindNONE/noDEEP2_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "Fid_shift" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_Fid_shift_NONE/Flag_SOM_Fid_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/Fid_shift_blindNONE/Fid_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    elif [ "${GoldSet}" == "Fid_oldmcorr" ]
    then 
      #Remove the '_oldmcorr' from the subsetting  /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_Fid_oldmcorr_NONE/Flag_SOM_Fid_NONE/g" configure_${GoldSet}.sh
      sed -i.bak "s/Fid_oldmcorr_blindNONE/Fid_blindNONE/g" configure_${GoldSet}.sh
      #/*fend*/}}}
    fi 
    
    #For each selection, link the shear data  /*fold*/ {{{
    mkdir -p GoldSet_${GoldSet}/
    cd ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/
    #Remove previous links  /*fold*/ {{{
    if [ -d PatchData ] 
    then 
      rm -fr PatchData 
    fi 
    #/*fend*/}}}
    
    mkdir -p PatchData/
    cd PatchData/
    ln -s ${ROOT}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/*} .
    #/*fend*/}}}
    
    cd ${ROOT}/COSMOPIPE
    #If we are running 'shift' chains, add the dz shifts to the config  /*fold*/ {{{
    if [ "$GoldSet" == "NONE_shift" ]
    then 
      sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.047 0.025 0.032 -0.004 -0.013'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.008 0.010  0.008  0.008'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSet" == "Fid_shift" ] 
    then 
      sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.000 0.002 0.013 0.011 -0.006'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.012 0.012 0.008  0.010'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSet" == "nozCOSMOS_shift" ] 
    then 
      sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.005 0.005 0.032 0.030 0.002'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.026 0.016 0.014 0.010  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSet" == "noVVDS_shift" ] 
    then 
      sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.001 0.001 0.024 0.014 -0.007'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.012 0.014 0.010  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSet" == "noDEEP2_shift" ] 
    then 
      sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='-0.001 0.002 -0.002 -0.009 -0.015'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.012 0.012 0.010  0.010'  #/g" configure_${GoldSet}.sh
    fi 
    #/*fend*/}}}
    
    #Correct the mbias values for each Gold Set  /*fold*/ {{{
    GoldSetLink=${GoldSet//_shift/}
    if [ "$GoldSetLink" == "Fid" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0145 -0.0176 -0.0125 0.0045 0.0122'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSetLink" == "noDEEP2" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0137 -0.0162 -0.0112 0.0054 0.0130'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSetLink" == "noVVDS" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0143 -0.0172 -0.0116 0.0047 0.0125'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSetLink" == "nozCOSMOS" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0143 -0.0159 -0.0106 0.0053 0.0135'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSetLink" == "speczquality4" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0141 -0.0163 -0.0121 0.0043 0.0115'  #/g" configure_${GoldSet}.sh
    elif [ "$GoldSetLink" == "multispec3" ] 
    then 
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.0158 -0.0203 -0.0173 -0.0033 -0.0012'  #/g" configure_${GoldSet}.sh
    fi 
    #/*fend*/}}}
    
  done
  cd ${ROOT}/
  #/*fend*/}}}
fi 

if [ "${1}" == "10" -o "${1}" == "" ]
then
  #Run the Gold Samples  /*fold*/ {{{
  cd ${ROOT}/${INPUTDIR}/
  #Download and decompress the shear data products  /*fold*/ {{{
  if [ ! -f KV450_COSMIC_SHEAR_DATA_RELEASE.tar.gz ]
  then 
    wget http://kids.strw.leidenuniv.nl/cs2018/KV450_COSMIC_SHEAR_DATA_RELEASE.tar.gz 
    tar -xf KV450_COSMIC_SHEAR_DATA_RELEASE.tar.gz
  fi 
  #/*fend*/}}}
  
  for GoldSet in ${GOLDLIST}
  do
    cd ${ROOT}/COSMOPIPE/
    #Run the configure script  /*fold*/ {{{
    #Remove any previous run scripts  /*fold*/ {{{
    if [ -f run_COSMOLOGY_PIPELINE.sh ]
    then 
      rm -f run_COSMOLOGY_PIPELINE.sh
    fi 
    #/*fend*/}}}
    
    bash configure_${GoldSet}.sh 
    mv -f run_COSMOLOGY_PIPELINE.sh run_COSMOLOGY_PIPELINE_${GoldSet}.sh 
    #/*fend*/}}}
    
    #Prepare the Redshift distributions for the run  /*fold*/ {{{
    mkdir -p ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/
    cd ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/
    #Remove any previous Nz  /*fold*/ {{{
    for tomo in `seq 5`
    do
      if [ -f ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}${tomo}_Nz.asc ]
      then 
        rm ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}${tomo}_Nz.asc
      fi 
    done
    #/*fend*/}}}
    
    if [ "${GoldSet}" == "NONE" -o "${GoldSet}" == "NONE_shift" -o "${GoldSet}" == "NONE_catastrophic" ]
    then 
      #Copy the KV450 DIR Redshift distributions to the NONE runs  /*fold*/ {{{
      ln -s ${ROOT}/${INPUTDIR}/KV450_COSMIC_SHEAR_DATA_RELEASE/REDSHIFT_DISTRIBUTIONS/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.1t0.3.asc \
        ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}1_Nz.asc
      ln -s ${ROOT}/${INPUTDIR}/KV450_COSMIC_SHEAR_DATA_RELEASE/REDSHIFT_DISTRIBUTIONS/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.3t0.5.asc \
        ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}2_Nz.asc
      ln -s ${ROOT}/${INPUTDIR}/KV450_COSMIC_SHEAR_DATA_RELEASE/REDSHIFT_DISTRIBUTIONS/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.5t0.7.asc \
        ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}3_Nz.asc
      ln -s ${ROOT}/${INPUTDIR}/KV450_COSMIC_SHEAR_DATA_RELEASE/REDSHIFT_DISTRIBUTIONS/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.7t0.9.asc \
        ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}4_Nz.asc
      ln -s ${ROOT}/${INPUTDIR}/KV450_COSMIC_SHEAR_DATA_RELEASE/REDSHIFT_DISTRIBUTIONS/Nz_DIR/Nz_DIR_Mean/Nz_DIR_z0.9t1.2.asc \
        ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}5_Nz.asc
      #/*fend*/}}}
    else 
      #Copy the new SOM Redshift distributions to the Gold runs  /*fold*/ {{{
      GoldSetLink=${GoldSet//_shift/}
      GoldSetLink=${GoldSetLink//_oldmcorr/}
      rm -f *${GoldSetLink}*_Nz.asc
      ln -s ${ROOT}/${OUTPUTDIR}/*${GoldSetLink}*_Nz.asc . 
      #/*fend*/}}}
    fi 
    #/*fend*/}}}
    
    #Copy the KV450 Survey Mask  /*fold*/ {{{
    if [ ! -f ${MASKFILE} ]
    then 
      ln -s ${ROOT}/${INPUTDIR}/${MASKFILE} ${MASKFILE}
    fi 
    cd ${ROOT}/COSMOPIPE/
    #/*fend*/}}}
    
    #IF running a catastrophic run, lower the number of live points  /*fold*/ {{{
    if [ "${GoldSet}" == "NONE_catastrophic" ]
    then 
      sed -i.bak "s/NS_n_live_points 1000 /NS_n_live_points 400 /g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/run_MCMC.sh 
    fi 
    #/*fend*/}}}
    
    #Check if we can launch another run  /*fold*/ {{{
    while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge ${MAXRUNS} ]
    do
      #If this is the first loop of the wait, then print what is running  /*fold*/ {{{
      if [ "${prompt}" != "${GoldSet}" ]
      then
        echo "Paused before starting ${GoldSet}: Maximum simultaneous runs reached (`date`)"
        prompt=${GoldSet}
      fi
      sleep 600
      #/*fend*/}}}
    done
    echo "Launching GoldSet ${GoldSet}: `ps au | grep -v 'bash -c ' | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge ${MAXRUNS} (`date`)" 
    #/*fend*/}}}
    
    #Run the main script  /*fold*/ {{{
    screen -S CosmoWrapper_Goldset_${GoldSet}_$$.sh -d -m bash -c "nice bash run_COSMOLOGY_PIPELINE_${GoldSet}.sh > run_COSMOLOGY_PIPELINE_${GoldSet}.log 2>&1"
    sleep 1
    #/*fend*/}}}
  done
  #/*fend*/}}}

  #Check if we can continue  /*fold*/ {{{
  while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge ${MAXRUNS} ]
  do
    #If this is the first loop of the wait, then print what is running  /*fold*/ {{{
    if [ "${prompt}" != "${GoldSet}" ]
    then
      echo "Waiting before post-processing (`date`)"
      prompt=${GoldSet}
    fi
    sleep 600
    #/*fend*/}}}
  done
  echo "Running Post-processing (`date`)" 
  #/*fend*/}}}
fi 
    
#Construct the output figures and tables
if [ "${1}" == "11" -o "${1}" == "" ]
then 
  #Construct the Chain Directories 
  bash setup_chaindir.sh
  #Make the contour plot 
  ${P_RSCRIPT} plot_Om_S8.R 
  #Make the stackplot 
  ${P_RSCRIPT} stackplot.R > S8_constraints.txt
  #Construct the marginal constraints 
  ${P_RSCRIPT} compare_marginals.R > marginal_constraints.txt 
fi 


