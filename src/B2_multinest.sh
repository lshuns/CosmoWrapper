#!/usr/bin/env bash

source vandenBusch2022.param
export OMP_NUM_THREADS=${MAXTHREADS}
ROOT=`pwd`

#Run the Gold Samples  /*fold*/ {{{
for GoldSet in Fid_noB2 Fid_onlyB2
do
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
    screen -S CosmoWrapper_Goldset_${GoldSet}_$$.sh -d -m bash -c "nice bash chainsB2/run_COSMOLOGY_PIPELINE_${GoldSet}.sh > chainsB2/run_COSMOLOGY_PIPELINE_${GoldSet}.log 2>&1"
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
