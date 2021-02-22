#!/usr/bin/env bash
#
# Script for running the fiducial results for the Wright 2020 (b) paper. 
# Uses KV450 Catalogues to run Cosmological Chains with different GoldClasses
#

set -e 

#Source the default parameters file /*fold*/ {{{
source vandenBusch2020b.param
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
      -sc 64 --only.som \
      -o ${OUTPUTDIR} -of ${SOMFILE//_SOMdata/} \
      --zr.label Z_B --zt.label Zbest \
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
        -sc 64 --optimise --refr.flag -cr recal_weight \
        --old.som ${OUTPUTDIR}/${SOMFILE} \
        -o ${OUTPUTDIR} -of ${SOMFILE//_SOMdata/_TOMO${TOMO}} \
        --zr.label Z_B --zt.label Zbest \
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
      if [ "${goldset}" != "ORIG" ]
      then
        if [ "`echo ${goldset} | grep -c NONE`" == "0" ] 
        then 
          GOLDFLAGLIST=`echo $GOLDFLAGLIST Flag_SOM_${goldset//_shift/}`
        fi 
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
      ${P_LDACFILTER} -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
    else 
      #RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
      ${P_LDACFILTER} -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
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
  if [ ! -d CosmoPipe ]
  then 
    git clone https://github.com/AngusWright/CosmoPipe.git
  else 
    cd ${ROOT}/INSTALL/CosmoPipe
    git stash 
    git pull
  fi
  cd ${ROOT}/INSTALL/CosmoPipe
  git checkout ldacfilter.py
  cd ${ROOT}
  
  bash INSTALL/CosmoPipe/COSMOPIPE_MASTER_INSTALL.sh \
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

  # fixed the file paths in the scripts
  sed -i.bak 's/reweight_${RECALGRID}/${FILEBODY}/g' ${ROOT}/COSMOPIPE/configure.sh
  sed -i.bak "s/^RECALGRID=3x4x4/FILEBODY=V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2  #/g" ${ROOT}/COSMOPIPE/configure.sh
  sed -i.bak 's/\\\@RECALGRID\\\@#${RECALGRID}/\\\@FILEBODY\\\@#${FILEBODY}/g' ${ROOT}/COSMOPIPE/configure.sh
  
  sed -i.bak "s/reweight_\@RECALGRID\@/\@FILEBODY\@/g" ${ROOT}/INSTALL/CosmoPipe/config/kv450_cf_likelihood_public/kv450_cf_likelihood_public.data
  
  sed -i.bak 's/reweight_${RECALGRID}/${FILEBODY}/g' ${ROOT}/INSTALL/CosmoPipe/scripts/configure_raw.sh
  sed -i.bak 's/\\\@RECALGRID\\\@#${RECALGRID}/\\\@FILEBODY\\\@#${FILEBODY}/g' ${ROOT}/INSTALL/CosmoPipe/scripts/configure_raw.sh

  sed -i.bak "s/reweight_\@RECALGRID\@/\@FILEBODY\@/g" \
    ${ROOT}/INSTALL/CosmoPipe/scripts/{calculate_2ptcorr.sh,create_cov_mat_corr_func.py,post_process_MCMC.sh,prepare_xis.sh,rebinned_corr_func_data_vec.py,run_COSMOLOGY_PIPELINE_raw.sh,run_covariance.sh}

  # this is a bit more complicated: the naming convetion has changed in K1000 and there is no recalgrid anymore
  sed -i.bak "s/\@RECALGRID\@/\@FILEBODY\@/g" ${ROOT}/INSTALL/CosmoPipe/scripts/{calculate_2ptcorr.sh,run_COSMOLOGY_PIPELINE_raw.sh}
  sed -i.bak "s/^RECALGRID=\@RECALGRID\@/FILEBODY=\@FILEBODY\@/g" ${ROOT}/INSTALL/CosmoPipe/scripts/{calculate_2ptcorr.sh,configure_raw.sh}
  sed -i.bak 's/\@SURVEY\@_%s_reweight_%s\@FILESUFFIX\@.cat/\@SURVEY\@_%s_%s\@FILESUFFIX\@.cat/g' ${ROOT}/INSTALL/CosmoPipe/scripts/ave_e_vs_ZB_rot.py

  #Update the configure files  /*fold*/ {{{
  for GoldSet in ${GOLDLIST} 
  do
    cd ${ROOT}/COSMOPIPE
    #Setup the Gold Set configurations  /*fold*/ {{{
    sed "s/@@GOLDSET@@/${GoldSet}/g" configure.sh > configure_${GoldSet}.sh 
    #/*fend*/}}}
    
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
    if [ "${GoldSet}" == "ORIG" ]
    then
      ################## LINK THE DR4 Photometric data Here ############################
      ${DIR_LDAC}/ldacdelkey${THELI} \
        -i /net/fohlen12/home/awright/KiDS/DR4/DIR/Iteration2/GoldClasses/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses.cat \
        -t OBJECTS -k PSFeabs \
        -o ${PHOTCAT_ALL_GOLD}
      ###Separate the N & S
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
          ${P_LDACFILTER} -i ${PHOTCAT_ALL_GOLD} -o ${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
            -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
        else 
          #RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
          ${P_LDACFILTER} -i ${PHOTCAT_ALL_GOLD} -o ${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
            -c "((ALPHA_J2000>${TMPLO})OR(ALPHA_J2000<=${TMPHI}));"
        fi 
      done
      #/*fend*/}}}
    else
      ln -s ${ROOT}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/*} .
    fi
    #/*fend*/}}}
    
    cd ${ROOT}/COSMOPIPE
    # sed -i.bak "s/^/  #/g" configure_${GoldSet}.sh
    sed -i.bak "s|^COSMOFISHER=/path/to/CosmoFisherForecast/|COSMOFISHER=/net/home/fohlen12/awright/src/CosmoFisherForecast/  #|g" configure_${GoldSet}.sh
    sed -i.bak "s/^DATE=/DATE=$(date +%Y-%m-%d)  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^ALLPATCH=/ALLPATCH=${PATCHALL}  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^PATCHLIST=/PATCHLIST=\"${PATCHLIST}\"  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^SURVEY=/SURVEY=K1000  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^SURVEYAREA=/SURVEYAREA=2.79864e+06  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.000 0.002 0.013 0.011 -0.006'  #/g" configure_${GoldSet}.sh
    if [ "${GoldSet}" == "noDz" ]
    then
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.0 0.0 0.0 0.0 0.0'  #/g" configure_${GoldSet}.sh
    else
      sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.011 0.012 0.008 0.010'  #/g" configure_${GoldSet}.sh
    fi
    sed -i.bak "s/^BLIND=/BLIND=${BLINDS}  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/_blindNONE_TOMO/_blind${BLINDS}_TOMO  #/g" configure_${GoldSet}.sh
    sed -i.bak "s|^STORAGEPATH=|STORAGEPATH=GoldSet_${GoldSet}/${STORAGEDIR}//  #|g" configure_${GoldSet}.sh
    if [ "${GoldSet}" == "XIcuts" ]
    then
      sed -i.bak "s/^XIPLUSLIMS=/XIPLUSLIMS='0.5 300.0'  #/g" configure_${GoldSet}.sh
      sed -i.bak "s/^XIMINUSLIMS=/XIMINUSLIMS='4.0 300.0'  #/g" configure_${GoldSet}.sh
    fi
    # catalogue specific
    sed -i.bak "s/^FILESUFFIX=/FILESUFFIX=_goldclasses  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^MASKFILE=/MASKFILE=KiDS_K1000_healpix.fits  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^E1VAR=/E1VAR='autocal_e1_${BLINDS}'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^E2VAR=/E2VAR='autocal_e2_${BLINDS}'  #/g" configure_${GoldSet}.sh  
    sed -i.bak "s/^XPIX=/XPIX='Xpos'  #/g" configure_${GoldSet}.sh  
    sed -i.bak "s/^YPIX=/YPIX='Ypos'  #/g" configure_${GoldSet}.sh  
    sed -i.bak "s/^GAAPFLAG=/GAAPFLAG='FLAG_GAAP_ugriZYJHKs'  #/g" configure_${GoldSet}.sh  
    if [ "${GoldSet}" == "ORIG" ] || [ "${GoldSet}" == "XIcuts" ] || [ "${GoldSet}" == "noDz" ]
    then
      ######### USE THE Original Fiducial ###########
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_Fid_${BLINDS}  #/g" configure_${GoldSet}.sh
    else
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_${GoldSet}_${BLINDS}  #/g" configure_${GoldSet}.sh
    fi
    sed -i.bak "s/^WEIGHTNAME=/WEIGHTNAME=recal_weight_${BLINDS}  #/g" configure_${GoldSet}.sh
    #/*fend*/}}}
    
    #Correct the mbias values for each Gold Set  /*fold*/ {{{
    sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.009 -0.011 -0.015 0.002 0.007'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^MBIASERRORS=/MBIASERRORS='0.019 0.020 0.017 0.012 0.010'  #/g" configure_${GoldSet}.sh
    #/*fend*/}}}
    
  done
  cd ${ROOT}/
  #/*fend*/}}}
fi 

if [ "${1}" == "10" -o "${1}" == "" ]
then
  #Run the Gold Samples  /*fold*/ {{{
  cd ${ROOT}/${INPUTDIR}/
  
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
   
    # fix CosmoFisherForecast
    cd ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/CosmoFisherForecast/psvareos/proc
    sed -i.bak "s/nmaxline_surveywindow 10000/nmaxline_surveywindow 20000/g" thps_func.c
    make > ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/rebuild_CosmoFisherForecast.log 2>&1
    cd ${ROOT}/COSMOPIPE/
    
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
    
    #Copy the new SOM Redshift distributions to the Gold runs  /*fold*/ {{{
    if [ "${GoldSet}" == "ORIG" ]
    then
      rm -f *_Nz.asc
      for fpath in /net/home/fohlen12/hendrik/KiDS/data_store/KiDS-1000/SOM_N_of_Z/K1000_NS*blind${BLINDS}*_Nz.asc
      do
        fbase=$(basename $fpath)
        # remove the commented header in the n(z) file from the data store
        tail -n +2 $fpath > $(echo "$fbase" | sed "s/_Fid_/_${GoldSet}_/g")
      done
    elif [ "${GoldSet}" == "XIcuts" ] || [ "${GoldSet}" == "noDz" ]
    then
      rm -f *_Nz.asc
      for fpath in ${ROOT}/${OUTPUTDIR}/*Fid*_Nz.asc
      do
        fbase=$(basename $fpath)
        ln -s ${fpath} $(echo "$fbase" | sed "s/_Fid_/_${GoldSet}_/g")
      done
    else
      GoldSetLink=${GoldSet//_shift/}
      GoldSetLink=${GoldSetLink//_oldmcorr/}
      rm -f *${GoldSetLink}*_Nz.asc
      ln -s ${ROOT}/${OUTPUTDIR}/*${GoldSetLink}*_Nz.asc . 
    fi
    #/*fend*/}}}
    #/*fend*/}}}
    
    #Copy the KV450 Survey Mask  /*fold*/ {{{
    if [ ! -f ${MASKFILE} ]
    then 
      ln -s ${ROOT}/${INPUTDIR}/${MASKFILE} ${MASKFILE}
    fi 
    cd ${ROOT}/COSMOPIPE/
    #/*fend*/}}}
  done
  #/*fend*/}}}
fi
    
if [ "${1}" == "11" -o "${1}" == "" ]
then
  #Run the Gold Samples  /*fold*/ {{{
  for GoldSet in ${GOLDLIST}
  do
    cd ${ROOT}/COSMOPIPE/
    #lower the number of live points for preliminary chains /*fold*/ {{{
    sed -i.bak "s/NS_n_live_points 1000 /NS_n_live_points 200 /g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/run_MCMC.sh 
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
  while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge 1 ]
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
if [ "${1}" == "12" -o "${1}" == "" ]
then 
  export BLINDS
  export STORAGEDIR
  export GOLDLIST
  export SURVEY="K1000"
  #Construct the Chain Directories 
  bash setup_chaindir.sh
  #Make the contour plot 
  ${P_RSCRIPT} plot_Om_S8.R 
  #Make the stackplot 
  ${P_RSCRIPT} stackplot.R > S8_constraints.txt
  #Construct the marginal constraints 
  ${P_RSCRIPT} compare_marginals.R > marginal_constraints.txt 
fi 


