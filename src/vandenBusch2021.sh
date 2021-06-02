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
    ${P_RSCRIPT} get_SOM_cellcounts.R \
      ${OUTPUTDIR}/${SOMFILE} \
      ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT//.fits/_allgoldclass.fits}
    ${P_RSCRIPT} get_SOM_cellcounts.R \
      ${OUTPUTDIR}/${SOMFILE} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits}
  else 
    echo "Gold Catalogues Already Exists! Skipping!" 
    echo "(${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits})"
  fi 
  #/*fend*/}}}
fi

#### this is only for the m-bias calibration ####
if [ "${1}" == "5a" -o "${1}" == "" ]
then
  outputdirMBIAS="output_for_mbias"
  ##Construct the Gold Classes  /*fold*/ {{{
  if [ ! -d ${outputdirMBIAS}/ ]
  then
    mkdir ${outputdirMBIAS}
    for cosmosfile in /net/home/fohlen12/awright/KiDS/DR4/DIR/Iteration1/COSMOS/COSMOSadaptdepth_ugri.V0.5.9A_ZYJHK.V2.0_photoz_LF_mask_SG_recalweight{,_goodgals_3x4x4}.cat
    do
      cosmosbase=$(basename ${cosmosfile})
      adaptfile=${outputdirMBIAS}/${cosmosbase//.cat/_adapt.fits}
      echo "==> processing ${cosmosbase}"
      echo "Constructing COSMOS Adapt Catalogue {"
      ${P_RSCRIPT} construct_adapt_catalogue.R ${cosmosfile} ${adaptfile}
      echo "} - Done" 
      echo "Constructing the COSMOS Goldclass subsets" 
      time ${P_RSCRIPT} construct_dr4_goldclasses.R \
        -p ${adaptfile} \
        -s ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
        --som ${OUTPUTDIR}/${SOMFILE} \
        --blinds "NONE" \
        --outputpath ${outputdirMBIAS}/ \
        --nzformat .asc \
        >> ${outputdirMBIAS}/construct_cosmos_goldclasses.log
      done
  else 
    echo "Output Folder Already Exist! Skipping!" 
    echo "(${outputdirMBIAS})"
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
      if [ "${goldset}" != "ORIG" ] && [ "${goldset}" != "ORIGmbias" ]
      then
        if [ "`echo ${goldset} | grep -c NONE`" == "0" ] 
        then 
          GoldSetLink=${goldset//_nodz/}
          GoldSetLink=${GoldSetLink//_trunc3/}
          GoldSetLink=${GoldSetLink//_trunc/}
          GoldSetLink=${GoldSetLink//_HartleyBest/}
          GoldSetLink=${GoldSetLink//_HartleyWorst/}
          GoldSetLink=${GoldSetLink//_Hartley/}
          GOLDFLAGLIST=`echo $GOLDFLAGLIST Flag_SOM_${GoldSetLink}`
        fi 
      fi
    done
    python merge_ldac_and_fits.py ${OUTPUTDIR}/${PHOTCAT_ALL} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} "${BLINDS}" "${GOLDFLAGLIST}" 
  else 
    echo "Merged catalogue already exists! Skipping!" 
    echo "(${OUTPUTDIR}/${PHOTCAT_ALL_GOLD})"
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

    GoldSetLink=${GoldSet//_nodz/}
    GoldSetLink=${GoldSetLink//_trunc3/}
    GoldSetLink=${GoldSetLink//_trunc/}
    GoldSetLink=${GoldSetLink//_HartleyBest/}
    GoldSetLink=${GoldSetLink//_HartleyWorst/}
    GoldSetLink=${GoldSetLink//_Hartley/}
    GoldSetTail=${GoldSet//$GoldSetLink/}
    if [ "${GoldSetTail}" != "" ]
    then 
      #Remove the '_shift' from the subsetting & Nz /*fold*/ {{{
      sed -i.bak "s/Flag_SOM_${GoldSet}/Flag_SOM_${GoldSetLink}/g" configure_${GoldSet}.sh
      sed -i.bak "s/${GoldSet}_blind${BLINDS}/${GoldSetLink}_blind${BLINDS}/g" configure_${GoldSet}.sh
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
    cp /net/home/fohlen12/awright/KiDS/Cosmo/PatchData/SOM_cov_multiplied.asc ./
    if [ "${GoldSet}" == "ORIG" ] || [ "${GoldSet}" == "ORIGmbias" ]
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
    if [ "${GoldSetLink}" == "Fid" ]
    then
      # bin1: -0.0016485259165871778 0.008904576570799095
      # bin2: -0.003807666962710251 0.006404266275882894
      # bin3: -0.006643348533635084 0.006419427909272947
      # bin4: 0.009886030380149706 0.006303345896370764
      # bin5: 0.012568191803766029 0.00805985325907987
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.002 -0.004 -0.007  0.010  0.013'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "Fid_noDEVILS" ]
    then
      # bin1: -0.0006296010406884851 0.008959283310110677
      # bin2: -0.0035735855359605093 0.006388020567536477
      # bin3: -0.006306263945970698 0.006469089845731768
      # bin4: 0.009919566143734643 0.006314314949824238
      # bin5: 0.012115237814203416 0.008050799624364158
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.001 -0.004 -0.006  0.010  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "plusPAUS" ]
    then
      # bin1: -0.0023269334232442337 0.008881225542425614
      # bin2: -0.0029790066285479925 0.006441456034541765
      # bin3: -0.0070332933389989985 0.006332818762900888
      # bin4: 0.010782697374635396 0.006191342761070044
      # bin5: 0.01278004986505398 0.008077782673125838
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.002 -0.003 -0.007  0.011  0.013'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "plusPAUSCOS15" ]
    then
      # bin1: -0.002792870134373022 0.008812701227445935
      # bin2: -0.0016483355948774612 0.006352654044074562
      # bin3: -0.009192201411269077 0.006355127294115379
      # bin4: 0.011582944825210256 0.006148865119996802
      # bin5: 0.012778279589247288 0.008105332732460821
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.003 -0.002 -0.009  0.012  0.013'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "nQ4" ]
    then
      # bin1: -0.0026642073311282087 0.008818247240945178
      # bin2: -0.0032242871494365362 0.006441279974643757
      # bin3: -0.006050157239036261 0.006465505709428217
      # bin4: 0.00972924294508508 0.006306496528352218
      # bin5: 0.012057486009379663 0.008172350342013841
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.003 -0.003 -0.006  0.010  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "ORIGmbias" ]
    then
      # bin1: -0.00995963912967607 0.008829221143656404
      # bin2: -0.00948511955203911 0.006330450857610293
      # bin3: -0.010991440899624444 0.006405850945378147
      # bin4: 0.007935261534419119 0.0064177568310609875
      # bin5: 0.012048666602049432 0.008072118337971662
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.010 -0.009 -0.011  0.008  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "onlyPAUS" ]
    then
      # bin1: -0.00458597126644202 0.008741059158244661
      # bin2: -0.003643567173407195 0.006439610824645876
      # bin3: -0.008529529637254405 0.0064225852457679075
      # bin4: 0.007650934867507904 0.00630644123981858
      # bin5: 0.008388372515304082 0.008299943911760298
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.005 -0.004 -0.009  0.008  0.008'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "onlyCOS15" ]
    then
      # bin1: -0.005770889954998207 0.008793911790837761
      # bin2: -0.0020081134463004425 0.006317168825689079
      # bin3: -0.009806280648000813 0.006403457221758943
      # bin4: 0.010404674721148812 0.00615987688108042
      # bin5: 0.012236512101840432 0.008077742500524721
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.006 -0.002 -0.010  0.010  0.012'  #/g" configure_${GoldSet}.sh
    elif [ "${GoldSetLink}" == "onlyPAUSCOS15" ]
    then
      # bin1: -0.006582387963870873 0.008872398780046082
      # bin2: -0.002017237730717207 0.006312960389114971
      # bin3: -0.009500122495017251 0.006392180202636211
      # bin4: 0.010547957903941751 0.0061493420460075375
      # bin5: 0.012232415677480046 0.008077742500524721
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.007 -0.002 -0.010  0.011  0.012'  #/g" configure_${GoldSet}.sh
    else
      # K1000 fiducial (Asgari et al. 2021)
      # bin1: -0.009 0.019
      # bin2: -0.011 0.020
      # bin3: -0.015 0.017
      # bin4: 0.002 0.012
      # bin5: 0.007 0.010
      sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.009 -0.011 -0.015  0.002  0.007'  #/g" configure_${GoldSet}.sh
    fi
    sed -i.bak "s/^MBIASERRORS=/MBIASERRORS=' 0.019  0.020  0.017  0.012  0.010'  #/g" configure_${GoldSet}.sh
    if [ "${GoldSet}" == "ORIG" ] || [ "${GoldSet}" == "ORIGmbias" ]
    then
      ######### USE THE Original Fiducial ###########
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_Fid_${BLINDS}  #/g" configure_${GoldSet}.sh
    else
      sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_${GoldSetLink}_${BLINDS}  #/g" configure_${GoldSet}.sh
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
    GoldSetLink=${GoldSet//_nodz/}
    GoldSetLink=${GoldSetLink//_trunc3/}
    GoldSetLink=${GoldSetLink//_trunc/}
    GoldSetLink=${GoldSetLink//_HartleyBest/}
    GoldSetLink=${GoldSetLink//_HartleyWorst/}
    GoldSetLink=${GoldSetLink//_Hartley/}
    GoldSetTail=${GoldSet//$GoldSetLink/}
    if [ "${GoldSet}" == "ORIG" ] || [ "${GoldSet}" == "ORIGmbias" ]
    then
      rm -f *_Nz.asc
      for fpath in /net/home/fohlen12/hendrik/KiDS/data_store/KiDS-1000/SOM_N_of_Z/K1000_NS*blind${BLINDS}*_Nz.asc
      do
        fbase=$(basename $fpath)
        # remove the commented header in the n(z) file from the data store
        tail -n +2 $fpath > $(echo "$fbase" | sed "s/_Fid_/_${GoldSet}_/g")
      done
    elif [ "${GoldSetTail}" == "_trunc" ]
    then
      rm -f *_${GoldSetLink}_blind*_Nz.asc
      for fpath in ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc
      do
        # clip n(z) files at 95-th percentile
        python ${ROOT}/nz_trim_percentile.py $fpath $(basename $fpath) 95.0
      done
    elif [ "${GoldSetTail}" == "_trunc3" ]
    then
      rm -f *_${GoldSetLink}_blind*_Nz.asc
      for fpath in ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc
      do
        # clip n(z) files at 97-th percentile
        python ${ROOT}/nz_trim_percentile.py $fpath $(basename $fpath) 97.0
      done
    elif [ "${GoldSetTail}" == "_nodz" ]
    then
      python ${ROOT}/nz_shift.py \
        -i ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc \
        -s 0.000 0.002 0.013 0.011 -0.006 \
        -o $(pwd)
    # offset the K1000 values with the data from Fig. 6 in Hartley+20
    # the mapping from the 4 DES bins to the 5 KiDS bins is a linear combination
    # map = [[   1    0    0    0]
    #        [ 1/2  1/2    0    0]
    #        [   0  2/3  1/3    0]
    #        [   0    0    1    0]
    #        [   0    0    0    1]]
    elif [ "${GoldSetTail}" == "_Hartley" ]
    then
      # KiDS prior:   0.000 0.002 0.013 0.011 -0.006
      # Hartley@KiDS: 0.008 0.015 0.014 -0.003 -0.058
      # KiDS-Hartley: -0.008 -0.013 -0.001 0.014 0.052
      python ${ROOT}/nz_shift.py \
        -i ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc \
        -s -0.008 -0.013 -0.001 0.014 0.052 \
        -o $(pwd)
    elif [ "${GoldSetTail}" == "_HartleyBest" ]
    then
      # KiDS prior:   0.000 0.002 0.013 0.011 -0.006
      # Hartley@KiDS: 0.008 0.015 0.014 -0.003 -0.058
      # KiDS+Hartley: 0.008 0.017 0.027 0.008 -0.064
      python ${ROOT}/nz_shift.py \
        -i ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc \
        -s 0.008 0.015 0.014 -0.003 -0.058 \
        -o $(pwd)
    elif [ "${GoldSetTail}" == "_HartleyWorst" ]
    then
      # KiDS prior:   0.000 0.002 0.013 0.011 -0.006
      # Hartley@KiDS: 0.012 0.014 0.004 -0.020 -0.160
      # KiDS+Hartley: 0.012 0.016 0.017 -0.009 -0.166
      python ${ROOT}/nz_shift.py \
        -i ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc \
        -s 0.012 0.014 0.004 -0.020 -0.160 \
        -o $(pwd)
    else
      rm -f *_${GoldSetLink}_blind*_Nz.asc
      ln -s ${ROOT}/${OUTPUTDIR}/*_${GoldSetLink}_blind*_Nz.asc . 
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
      sleep 200
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
    sleep 200
    #/*fend*/}}}
  done
  echo "Running Post-processing (`date`)" 
  #/*fend*/}}}
fi 

if [ "${1}" == "11a" -o "${1}" == "" ]
then
  #Run the Gold Samples with maxlike sampler for MAP /*fold*/ {{{
  P_PYTHON=${ROOT}/COSMOPIPE/INSTALL/miniconda3/bin/python3
  #install Marika's custom maxlike sampler
  installdir=${ROOT}/COSMOPIPE/INSTALL/kids1000_chains
  if [ ! -e ${installdir} ]
  then
    git clone git@bitbucket.org:marika_a/kids1000_chains.git ${installdir}
  else
    cd ${installdir}
    git pull
    cd ${ROOT}
  fi

  #Run the least square fitting
  for GoldSet in ${GOLDLIST}
  do
    maxlike_out=${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/LeastSquares/output/K1000_UNBLINDED/cosebis/
    if [ ! -e ${maxlike_out} ]
    then
      mkdir -p ${maxlike_out}
    fi
    #create runtime script that is executed in a screen
    runtime_script=${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/run_leastsquares.sh
    echo "#!/usr/bin/env bash
${P_PYTHON} ${installdir}/maxlike/maxlike_cosmosis.py \
  -i ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts//COSEBIs_chain.ini \
  -m ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/MCMC/output/K1000_UNBLINDED/cosebis/chain/output_multinest_${BLINDS}.txt \
  -o ${maxlike_out}/output_${BLINDS}.txt \
  --max_post -s \
  --best_fit_value ${maxlike_out}/best_fit_value.txt \
  --best_fit_priors ${maxlike_out}/best_fit_priors.txt \
  --maxiter 3000" > ${runtime_script}

    #Check if we can launch another run  /*fold*/ {{{
    while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_leastsquares` -ge ${MAXRUNS} ]
    do
      #If this is the first loop of the wait, then print what is running  /*fold*/ {{{
      if [ "${prompt}" != "${GoldSet}" ]
      then
        echo "Paused before starting ${GoldSet}: Maximum simultaneous runs reached (`date`)"
        prompt=${GoldSet}
      fi
      sleep 60
      #/*fend*/}}}
    done
    echo "Launching GoldSet ${GoldSet}: `ps au | grep -v 'bash -c ' | grep -v grep | grep -c run_leastsquares` -ge ${MAXRUNS} (`date`)" 
    #/*fend*/}}}
    
    #Run the main script  /*fold*/ {{{
    logfile=${maxlike_out}/cosebis_leastsquares_output.log
    screen -S CosmoWrapper_Goldset_${GoldSet}_$$.sh -d -m bash -c "nice bash ${runtime_script} > ${logfile} 2>&1"
    sleep 1
    #/*fend*/}}}
  done
  #/*fend*/}}}

  #Check if we can continue  /*fold*/ {{{
  while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_leastsquares` -ge 1 ]
  do
    #If this is the first loop of the wait, then print what is running  /*fold*/ {{{
    if [ "${prompt}" != "${GoldSet}" ]
    then
      echo "Waiting before post-processing (`date`)"
      prompt=${GoldSet}
    fi
    sleep 60
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
    echo "#### collecting ${GoldSet} ####"
    res_storage=${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}
    mcmc_out=${res_storage}/MCMC/output/K1000_UNBLINDED/cosebis/chain/output_multinest_${BLINDS}.txt
    if [ -e ${mcmc_out} ]
    then
      cp -fv ${mcmc_out} GoldSet_${GoldSet}_multinest_${BLINDS}.txt
    else
      echo "no MCMC chains found, skipping ..."
    fi

    lsq_out=${res_storage}/LeastSquares/output/K1000_UNBLINDED/cosebis/output_${BLINDS}.txt
    if [ -e ${lsq_out} ]
    then
      cp -fv ${lsq_out} GoldSet_${GoldSet}_maxlike_${BLINDS}.txt
    else
      echo "no least square fits found, skipping ..."
    fi
  done
fi 


