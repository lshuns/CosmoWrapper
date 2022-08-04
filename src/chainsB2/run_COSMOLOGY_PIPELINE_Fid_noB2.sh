ROOT=`pwd`
GoldSet=Fid_noB2
${ROOT}/COSMOPIPE/INSTALL/miniconda3/bin/cosmosis \
    ${ROOT}/chainsB2/B2_chain_${GoldSet}.ini > ${ROOT}/chainsB2/B2_chain_${GoldSet}.log
