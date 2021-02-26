#!/usr/bin/env bash
#
# Script for running the fiducial results for the Wright 2020 (b) paper. 
# Uses KV450 Catalogues to run Cosmological Chains with different GoldClasses
#

set -e 

#Source the default parameters file /*fold*/ {{{
source vandenBusch2021.param
export OMP_NUM_THREADS=${MAXTHREADS}
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
    bash SOM_DIR/INSTALL.sh 2&>1 install_SOM_DIR.log
    cd ${ROOT}
    ln -s INSTALL/SOM_DIR/R/SOM_DIR.R . 
    #}}}
  fi 
  #Construct the Fiducial SOM  /*fold*/ {{{
  if [ ! -f ${OUTPUTDIR}/${SOMFILE} ]
  then 
    echo "Constructing the Fiducial SOM" 
    ${P_RSCRIPT} SOM_DIR.R \
      -r ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} -t ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
      --toroidal --topo hexagonal --som.dim 101 101 -np -fn Inf \
      -sc ${MAXTHREADS} --only.som \
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
        -sc ${MAXTHREADS} --optimise --refr.flag -cr recal_weight \
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
      python3 ${P_LDACFILTER} \
        -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
    else 
      #RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
      python3 ${P_LDACFILTER} \
        -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
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
  git checkout K1000
  cd ${ROOT}

  RESTORE_PATH=${PATH}
  if [ "P_GCC" != "" ]
  then
    export PATH=${P_GCC}:${RESTORE_PATH}
  fi
  bash INSTALL/CosmoPipe/COSMOPIPE_MASTER_INSTALL.sh \
    --noconfig \
    --packroot ${ROOT}/INSTALL/CosmoPipe/ \
    --runroot ${ROOT}/COSMOPIPE/ \
    --cosmofisher ~awright/KiDS/src/CosmoKiDS/ \
    --storagepath GoldSet_@@GOLDSET@@/${STORAGEDIR}/ \
    --runtime GoldSet_@@GOLDSET@@/RUNTIME/ \
    --configpath GoldSet_@@GOLDSET@@/RUNTIME/config/ \
    --scriptpath GoldSet_@@GOLDSET@@/RUNTIME/scripts/ \
    --nzfileid ${PHOTCAT_ALL_DCOL//.cat/_@@GOLDSET@@_blind${BLINDS}_TOMO} \
    --nzfilesuffix _Nz.asc \
    --patchpath ${ROOT}/COSMOPIPE/GoldSet_@@GOLDSET@@/PatchData/ \
    --filebody ${FILEBODY} \
    --filesuffix _goldclasses \
    --surveyarea ${SURVEYAREA} \
    --blind ${BLINDS} \
    --blinding UNBLINDED
  export PATH=${RESTORE_PATH}
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
    cp /net/home/fohlen12/awright/KiDS/Cosmo/PatchData/SOM_cov_multiplied.asc ./
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
          python3 ${P_LDACFILTER} \
            -i ${PHOTCAT_ALL_GOLD} -o ${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
            -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
        else 
          #RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
          python3 ${P_LDACFILTER} \
            -i ${PHOTCAT_ALL_GOLD} -o ${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
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
    sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.000 0.002 0.013 0.011 -0.006'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.011 0.012 0.008 0.010'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.009 -0.011 -0.015 0.002 0.007'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^MBIASERRORS=/MBIASERRORS='0.019 0.020 0.017 0.012 0.010'  #/g" configure_${GoldSet}.sh
    if [ "${GoldSet}" == "ORIG" ]
    then
      ######### USE THE Original Fiducial ###########
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_Fid_${BLINDS}  #/g" configure_${GoldSet}.sh
    else
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_${GoldSet}_${BLINDS}  #/g" configure_${GoldSet}.sh
    fi
    #/*fend*/}}}
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
    #sed -i.bak "s/NS_n_live_points 1000 /NS_n_live_points 200 /g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/run_MCMC.sh 
    sed -i.bak "s/live_points = 1000/live_points = ${NLIVE}/g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/COSEBIs_chain.ini
    sed -i.bak "s/tolerance = 0.01/tolerance = ${TOLERANCE}/g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/COSEBIs_chain.ini
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
  #Construct the Chain Directories 
  cd ${ROOT}
  #Construct the chain directory
  mkdir -p ChainDir
  cd ChainDir
  #Get the planck contour data
  if [ ! -f base-plikHM-TTTEEE-lowl-lowE_R3.00.zip ]
  then
    wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip -O base-plikHM-TTTEEE-lowl-lowE_R3.00.zip 
    unzip base-plikHM-TTTEEE-lowl-lowE_R3.00.zip
    cat base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE_?.txt > base_plikHM_TTTEEE_lowl_lowE.txt
    ln -s base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE.paramnames . 
  fi
  #Link the completed goldsets
  for GoldSet in ${GOLDLIST}
  do
    cp -fv ../COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/MCMC/output/K1000_UNBLINDED/cosebis/chain/output_multinest_${BLINDS}.txt \
          GoldSet_${GoldSet}_multinest_${BLINDS}.txt
  done
fi 


