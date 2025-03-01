#!/usr/bin/env bash
#
# Script for running the fiducial results for the Wright 2020 (b) paper.
# Uses KV450 Catalogues to run Cosmological Chains with different GoldClasses
#

set -e

# Source the default parameters file
source cosmowrapper.param
export OMP_NUM_THREADS=${MAXTHREADS}

# Define the root directory
ROOT=`pwd`

# Check for the input directory and files
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
  elif [ ! -f ${INPUTDIR}/${SOMCOVFILE} ]
  then
    echo "ERROR: The Input Directory does not contain the shift parameter covariance file!"
    exit 1
  fi
else
  echo "ERROR: The Input Directory (containing the spectroscopic & photometric catalogues, and survey footprint mask) does not exist!"
  exit 1
fi
if [ ! -d ${COSMOFISHER} ]
then
  echo "ERROR: CosmoFisher installation path does not exist!"
fi

# Make the output directory
mkdir -p ${OUTPUTDIR}


################################### STEP 01 ###################################
# Prepare the intput data: replace the original magnitudes by the adapted ones,
# convert to LDAC table format if necessary
photcat="${PHOTCAT_ALL%.*}.cat"
SPECCAT_ALL_ADAPT=${SPECCAT_ALL//.csv/_adapt.fits}
if [ "${1}" == "1" -o "${1}" == "" ]
then
  #Construct the Spectrosopic Adapt Catalogue
  if [ ! -f ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} ]
  then
    echo "Constructing Spectroscopic Adapt Catalogue {"
    ${P_RSCRIPT} construct_adapt_catalogue.R ${INPUTDIR}/${SPECCAT_ALL} ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT}
    echo "} - Done"
  else
    echo "Spectroscopic Adapt Catalogue Already Exists! Skipping!"
  fi
  #Convert to LDAC if FITS
  if [ ! -f ${OUTPUTDIR}/${photcat} ]
  then
    has_fields=$(echo $(${DIR_LDAC}/ldacdesc${THELI} -i ${INPUTDIR}/${PHOTCAT_ALL} | grep -c FIELDS))
    if [ "${has_fields}" == "1" ]
    then
      #Link file
      ln -sv ${ROOT}/${INPUTDIR}/${PHOTCAT_ALL} ${ROOT}/${OUTPUTDIR}/${photcat}
    else
      #Create a fake LDAC table
      python fits2ldac.py ${INPUTDIR}/${PHOTCAT_ALL} -o ${OUTPUTDIR}/${photcat}
    fi
  else
    echo "Photometric Catalogue Already Exists! Skipping!"
  fi
fi
# substitute to expected file extension .cat
PHOTCAT_ALL=$photcat

# Calculate Patch Numbers
PHOTCAT_PATCHES=""
NPATCH=0
for PATCH in ${PATCHLIST}
do
  let NPATCH=$NPATCH+1
  # Combine the patch filenames
  PHOTCAT_PATCHES="${PHOTCAT_PATCHES} ${INPUTDIR}/${PHOTCAT_ALL//${PATCHALL}/${PATCH}}"
done


################################### STEP 02 ###################################
# Keep a minimum set of columns in the shear catalogue
PHOTCAT_ALL_DCOL=${PHOTCAT_ALL//.cat/_DIRcols.cat}
if [ "${1}" == "2" -o "${1}" == "" ]
then
  # Select DIRcol subset (within LDAC)
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} ]
  then
    echo -e "Constructing the DIR Column Photometry Catalogue\n"
    list=`${DIR_LDAC}/ldacdesc${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL} -t OBJECTS |
      grep "Key name" | awk -F. '{print $NF}' |
      grep -v "SeqNr\|THELI_\|_B\|ID\|MAG_GAAP_\|${WEIGHTNAME}\|autocal\|MAG_AUTO"`
    ${DIR_LDAC}/ldacdelkey${THELI} -i ${OUTPUTDIR}/${PHOTCAT_ALL} -k ${list} -o ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL}
    echo "\n- Done"
  else
    echo "DIR Column Photometry Catalogue Already Exists! Skipping!"
  fi
fi


################################### STEP 03 ###################################
# Train the SOM on the calibration data
if [ "${1}" == "3" -o "${1}" == "" ]
then
  if [ ! -f SOM_DIR.R ]
  then
    # Install the SOM_DIR code {{{
    cd ${ROOT}
    mkdir -p INSTALL
    cd ${ROOT}/INSTALL
    if [ ! -d SOM_DIR ]
    then
      git clone https://github.com/AngusWright/SOM_DIR.git
    fi
    # bash SOM_DIR/INSTALL.sh 2&>1 install_SOM_DIR.log
    # installed through the INSTALL.sh of this repository
    cd ${ROOT}
    ln -s INSTALL/SOM_DIR/R/SOM_DIR.R .
  fi
  # Construct the Fiducial SOM
  if [ ! -f ${OUTPUTDIR}/${SOMFILE} ]
  then
    echo "Constructing the Fiducial SOM"
    ${P_RSCRIPT} SOM_DIR.R \
      -r ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} -t ${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
      --toroidal --topo hexagonal --som.dim 101 101 -np -fn Inf \
      -sc ${MAXTHREADS} --only.som \
      -o ${OUTPUTDIR} -of ${SOMFILE//_SOMdata/} \
      --zr.label Z_B --zt.label ${ZLABEL} \
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
fi


################################### STEP 04 ###################################
# Define the gold sample and compute the SOM redshift distributions
if [ "${1}" == "4" -o "${1}" == "" ]
then
  # Construct the Gold Classes
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
      -cv ${WEIGHTNAME} \
      -zn ${ZLABEL} \
      >> ${OUTPUTDIR}/construct_dr4_goldclasses.log
  else
    echo "Gold Catalogues Already Exists! Skipping!"
    echo "(${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits})"
  fi
fi


################################### STEP 05 ###################################
# Merge all SOM weights and gold flags into one catalogue
PHOTCAT_ALL_GOLD=${PHOTCAT_ALL//.cat/_goldclasses.cat}
if [ "${1}" == "5" -o "${1}" == "" ]
then
  # Merge the new GoldClasses back with the original catalogue
  if [ ! -f ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} ]
  then
    echo "Constructing merged catalogue"
    GOLDFLAGLIST=""
    for goldset in ${GOLDLIST}
    do
      GOLDFLAGLIST=`echo $GOLDFLAGLIST Flag_SOM_${goldset}`
    done
    python merge_ldac_and_fits.py ${OUTPUTDIR}/${PHOTCAT_ALL} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} \
      ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} "${BLINDS}" "${GOLDFLAGLIST}"
  else
    echo "Merged catalogue already exists! Skipping!"
    echo "(${OUTPUTDIR}/${PHOTCAT_ALL_GOLD})"
  fi
fi


################################### STEP 06 ###################################
# Split the gold sample catalogues into patches
if [ "${1}" == "6" -o "${1}" == "" ]
then
  # Recreate all patchwise catalogues
  for PATCHNUM in `seq ${NPATCH}`
  do
    # Get the patch definitions
    PATCH=`echo $PATCHLIST | awk -v n=${PATCHNUM} '{print $n}'`
    TMPLO=`echo $PATCHLO   | awk -v n=${PATCHNUM} '{print $n}'`
    TMPHI=`echo $PATCHHI   | awk -v n=${PATCHNUM} '{print $n}'`
    #Filter the ALLPATCH catalogue on RA
    if [ $TMPLO -lt $TMPHI ]
    then
      # Normal RA limits: pick lo < RA <= hi
      python3 ldacfilter.py \
        -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})AND(ALPHA_J2000<=${TMPHI}));"
    else
      # RA limits cross the RA=0 boundary: pick (RA > lo | RA <= hi)
      python3 ldacfilter.py \
        -i ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} -o ${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/${PATCH}} \
        -c "((ALPHA_J2000>${TMPLO})OR(ALPHA_J2000<=${TMPHI}));"
    fi
  done
fi


################################### STEP 07 ###################################
# Install CosmoPipe
if [ "${BLINDS}" != "NONE" ]
then
  WEIGHTNAME="${WEIGHTNAME}_${BLINDS}"  # fix the weight name with the blinds
fi
if [ "${1}" == "7" -o "${1}" == "" ]
then
  #Check if conda is active
  if [ ! -z "${CONDA_DEFAULT_ENV}" ]
  then
    echo "ERROR: disable all conda environments before installing:"
    echo "       conda deactivate"
    exit 1
  fi
  # Run the CosmoPipe Installation
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

  bash INSTALL/CosmoPipe/COSMOPIPE_MASTER_INSTALL.sh \
    --noconfig \
    --packroot ${ROOT}/INSTALL/CosmoPipe/ \
    --runroot ${ROOT}/COSMOPIPE/ \
    --cosmofisher ${COSMOFISHER} \
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
    --blinding UNBLINDED \
    --weightname ${WEIGHTNAME}
fi


################################### STEP 08 ###################################
# Configure CosmoPipe for the gold samples
if [ "${1}" == "8" -o "${1}" == "" ]
then
  # Update the configure files
  for GoldSet in ${GOLDLIST}
  do
    cd ${ROOT}/COSMOPIPE
    # Setup the Gold Set configurations
    sed "s/@@GOLDSET@@/${GoldSet}/g" configure.sh > configure_${GoldSet}.sh

    # For each selection, link the shear data
    mkdir -p GoldSet_${GoldSet}/
    cd ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/
    # Remove previous links
    if [ -d PatchData ]
    then
      rm -fr PatchData
    fi

    mkdir -p PatchData/
    cd PatchData/
    cp ${ROOT}/${INPUTDIR}/${SOMCOVFILE} ./
    ln -s ${ROOT}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD//${PATCHALL}/*} .

    cd ${ROOT}/COSMOPIPE
    # sed -i.bak "s/^/  #/g" configure_${GoldSet}.sh
    sed -i.bak 's/^NTHETABINXI=/NTHETABINXI="4000"  #/g' configure_${GoldSet}.sh
    sed -i.bak "s/^DZPRIORMU=/DZPRIORMU='0.000 0.002 0.013 0.011 -0.006'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^DZPRIORSD=/DZPRIORSD='0.010 0.011 0.012 0.008 0.010'  #/g" configure_${GoldSet}.sh
    # configure m-bias
    # bin1: -0.0016485259165871778 0.008904576570799095
    # bin2: -0.003807666962710251 0.006404266275882894
    # bin3: -0.006643348533635084 0.006419427909272947
    # bin4: 0.009886030380149706 0.006303345896370764
    # bin5: 0.012568191803766029 0.00805985325907987
    sed -i.bak "s/^MBIASVALUES=/MBIASVALUES='-0.002 -0.004 -0.007  0.010  0.013'  #/g" configure_${GoldSet}.sh
    sed -i.bak "s/^MBIASERRORS=/MBIASERRORS=' 0.019  0.020  0.017  0.012  0.010'  #/g" configure_${GoldSet}.sh
    # set the gold flag name
    sed -i.bak "s/^SHEARSUBSET=/SHEARSUBSET=Flag_SOM_${GoldSet}_${BLINDS}  #/g" configure_${GoldSet}.sh
  done
  cd ${ROOT}/
fi

################################### STEP 09 ###################################
# Link gold catalogues and redshift distributions into the CosmoPipe directory
if [ "${1}" == "9" -o "${1}" == "" ]
then
  # Run the Gold Samples
  cd ${ROOT}/${INPUTDIR}/

  for GoldSet in ${GOLDLIST}
  do
    cd ${ROOT}/COSMOPIPE/
    # Run the configure script
    # Remove any previous run scripts
    if [ -f run_COSMOLOGY_PIPELINE.sh ]
    then
      rm -f run_COSMOLOGY_PIPELINE.sh
    fi

    bash configure_${GoldSet}.sh
    mv -f run_COSMOLOGY_PIPELINE.sh run_COSMOLOGY_PIPELINE_${GoldSet}.sh

    # Prepare the Redshift distributions for the run
    mkdir -p ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/
    cd ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/
    # Remove any previous Nz
    for tomo in `seq 5`
    do
      if [ -f ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}${tomo}_Nz.asc ]
      then
        rm ${PHOTCAT_ALL_DCOL//.cat/_${GoldSet}_blindNONE_TOMO}${tomo}_Nz.asc
      fi
    done

    # Copy the new SOM Redshift distributions to the Gold runs
    rm -f *_${GoldSet}_blind*_Nz.asc
    ln -s ${ROOT}/${OUTPUTDIR}/*_${GoldSet}_blind*_Nz.asc .

    # Copy the KV450 Survey Mask
    if [ ! -f ${MASKFILE} ]
    then
      ln -s ${ROOT}/${INPUTDIR}/${MASKFILE} ${MASKFILE}
    fi
    cd ${ROOT}/COSMOPIPE/
  done
fi


################################### STEP 10 ###################################
# Run CosmoPipe: measurements, data vector, covariance and multinest sampling
if [ "${1}" == "10" -o "${1}" == "" ]
then
  #Run the Gold Samples
  for GoldSet in ${GOLDLIST}
  do
    cd ${ROOT}/COSMOPIPE/
    # lower the number of live points for preliminary chains
    sed -i.bak "s/live_points = 1000/live_points = ${NLIVE}/g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/COSEBIs_chain.ini
    sed -i.bak "s/tolerance = 0.01/tolerance = ${TOLERANCE}/g" ${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/RUNTIME/scripts/COSEBIs_chain.ini

    # Check if we can launch another run
    while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge ${MAXRUNS} ]
    do
      # If this is the first loop of the wait, then print what is running
      if [ "${prompt}" != "${GoldSet}" ]
      then
        echo "Paused before starting ${GoldSet}: Maximum simultaneous runs reached (`date`)"
        prompt=${GoldSet}
      fi
      sleep 200
    done
    echo "Launching GoldSet ${GoldSet}: `ps au | grep -v 'bash -c ' | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge ${MAXRUNS} (`date`)"

    # Run the main script
    screen -S CosmoWrapper_Goldset_${GoldSet}_$$.sh -d -m bash -c "nice bash run_COSMOLOGY_PIPELINE_${GoldSet}.sh > run_COSMOLOGY_PIPELINE_${GoldSet}.log 2>&1"
    sleep 1
  done

  # Check if we can continue
  while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_COSMOLOGY_PIPELINE` -ge 1 ]
  do
    # If this is the first loop of the wait, then print what is running
    if [ "${prompt}" != "${GoldSet}" ]
    then
      echo "Waiting before post-processing (`date`)"
      prompt=${GoldSet}
    fi
    sleep 200
  done
  echo "Running Post-processing (`date`)"
fi


################################### STEP 11 ###################################
# Run CosmoPipe: minimisation for best best-fit parameter estimate
if [ "${1}" == "11" -o "${1}" == "" ]
then
  # Run the Gold Samples with maxlike sampler for MAP
  P_PYTHON=${ROOT}/COSMOPIPE/INSTALL/miniconda3/bin/python3
  # install Marika's custom maxlike sampler
  installdir=${ROOT}/COSMOPIPE/INSTALL/kids1000_chains
  if [ ! -e ${installdir} ]
  then
    git clone git@bitbucket.org:marika_a/kids1000_chains.git ${installdir}
  else
    cd ${installdir}
    git pull
    cd ${ROOT}
  fi

  # Run the least square fitting
  for GoldSet in ${GOLDLIST}
  do
    maxlike_out=${ROOT}/COSMOPIPE/GoldSet_${GoldSet}/${STORAGEDIR}/LeastSquares/output/K1000_UNBLINDED/cosebis/
    if [ ! -e ${maxlike_out} ]
    then
      mkdir -p ${maxlike_out}
    fi
    # create runtime script that is executed in a screen
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

    # Check if we can launch another run
    while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_leastsquares` -ge ${MAXRUNS} ]
    do
      # If this is the first loop of the wait, then print what is running
      if [ "${prompt}" != "${GoldSet}" ]
      then
        echo "Paused before starting ${GoldSet}: Maximum simultaneous runs reached (`date`)"
        prompt=${GoldSet}
      fi
      sleep 60
    done
    echo "Launching GoldSet ${GoldSet}: `ps au | grep -v 'bash -c ' | grep -v grep | grep -c run_leastsquares` -ge ${MAXRUNS} (`date`)"

    # Run the main script
    logfile=${maxlike_out}/cosebis_leastsquares_output.log
    screen -S CosmoWrapper_Goldset_${GoldSet}_$$.sh -d -m bash -c "nice -n 10 bash ${runtime_script} > ${logfile} 2>&1"
    sleep 1
  done

  # Check if we can continue
  while [ `ps au | grep -v "bash -c " | grep -v grep | grep -c run_leastsquares` -ge 1 ]
  do
    # If this is the first loop of the wait, then print what is running
    if [ "${prompt}" != "${GoldSet}" ]
    then
      echo "Waiting before post-processing (`date`)"
      prompt=${GoldSet}
    fi
    sleep 60
  done
  echo "Running Post-processing (`date`)"
fi


################################### STEP 12 ###################################
# Download Planck TTTEEE and create basic S8 comparison plots.
if [ "${1}" == "12" -o "${1}" == "" ]
then
  # Construct the Chain Directories
  cd ${ROOT}
  # Construct the chain directory
  mkdir -p ChainDir
  cd ChainDir
  # Get the planck contour data
  if [ ! -f base-plikHM-TTTEEE-lowl-lowE_R3.00.zip ]
  then
    wget https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip -O base-plikHM-TTTEEE-lowl-lowE_R3.00.zip
    unzip base-plikHM-TTTEEE-lowl-lowE_R3.00.zip
    cat base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE_?.txt > base_plikHM_TTTEEE_lowl_lowE.txt
    ln -s base/plikHM_TTTEEE_lowl_lowE/base_plikHM_TTTEEE_lowl_lowE.paramnames .
  fi
  # Link the completed goldsets
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
