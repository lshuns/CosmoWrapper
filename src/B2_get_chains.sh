#!/usr/bin/env bash

source vandenBusch2022.param
export OMP_NUM_THREADS=${MAXTHREADS}
ROOT=`pwd`

cp -v ${ROOT}/COSMOPIPE/GoldSet_Fid/${STORAGEDIR}/MCMC/output/K1000_UNBLINDED/cosebis/chain_noB2/output_multinest_${BLINDS}.txt \
      ChainDir/GoldSet_Fid_noB2_multinest_${BLINDS}.txt
cp -v ${ROOT}/COSMOPIPE/GoldSet_Fid/${STORAGEDIR}/LeastSquares/output/K1000_UNBLINDED/cosebis/noB2/output_${BLINDS}.txt \
      ChainDir/GoldSet_Fid_noB2_maxlike_${BLINDS}.txt

cp -v ${ROOT}/COSMOPIPE/GoldSet_Fid/${STORAGEDIR}/MCMC/output/K1000_UNBLINDED/cosebis/chain_onlyB2/output_multinest_${BLINDS}.txt \
      ChainDir/GoldSet_Fid_onlyB2_multinest_${BLINDS}.txt
cp -v ${ROOT}/COSMOPIPE/GoldSet_Fid/${STORAGEDIR}/LeastSquares/output/K1000_UNBLINDED/cosebis/onlyB2/output_${BLINDS}.txt \
      ChainDir/GoldSet_Fid_onlyB2_maxlike_${BLINDS}.txt
