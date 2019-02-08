#!/bin/bash

#################################### Wrapper submit script for Upgrade production 
#Written by Alexis Kalogeropoulos - July 2014
#Adapted by Julie Hogan - summer 2016, jmhogan@fnal.gov

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530

startTime=`date +%s`

# Condor arguments
# Arguments = %(FILEIN)s %(OUTPUTDIR)s/%(RELPATH)s_%(PILEUP)s %(FILEOUT)s.root %(PILEUP)s
#FILEIN=$1
#OUTPUT=$2
#FILEOUT=$3
#PILEUP=$4
#SKIPEVT=$5
#MAXEVT=$6
#Arguments = %number of file

NUMBERFILE=$1
MASS=$2

FILEIN=/nfs/dust/cms/user/gunnep/GridpackGen/CMSSW_9_3_9/src/Output_inclusiveSingleAndPairS1-3-M${MASS}_${NUMBERFILE}/S13SingleAndPair-GenOnly-${NUMBERFILE}.root
#FILEIN=/nfs/dust/cms/user/gunnep/GridpackGen/CMSSW_9_3_9/src/Output_inclusiveSingleAndPairS1-3-${NUMBERFILE}/S13SingleAndPair-GenOnly-incl-${NUMBERFILE}.root
OUTPUT=user/pgunnell
FILEOUT=DelphesS13SingleAndPair-inclusive-30PU-${NUMBERFILE}.root
PILEUP=30PU
SKIPEVT=$5
MAXEVT=10

echo "Starting job on " `date`
echo "Running on " `uname -a`
echo "System release " `cat /etc/redhat-release`

if [[ $# -eq 4 ]] ; then
    echo "Setting SkipEvents to 0, no argument given"
    SKIPEVT=0
    echo "Setting MaxEvents to -1, no argument given"
    MAXEVT=-1
fi

# Set variables
#detCard=CMS_PhaseII_${PILEUP}_v03_splitter.tcl
detCard=CMS_PhaseII_${PILEUP}.tcl
energy=13
DelphesVersion=tags/3.4.2pre15
nPU=`echo $detCard | cut -d '_' -f 2 | cut -d '.' -f 1`
process=`echo $FILEIN | cut -d '_' -f 1-2`
configuration=`echo $detCard | cut -d '_' -f 1-2`
DelphesOutput=CMSP2_${nPU}_`echo $FILEIN`.root

#xrdcp -f root://cmseos.fnal.gov//store/user/snowmass/DelphesSubmissionLPCcondor/MinBias_100k.pileup .
#XRDEXIT=$?
#if [[ $XRDEXIT -ne 0 ]]; then
#    echo "exit code $XRDEXIT, failure in xrdcp of MinBias_100k.pileup"
#    exit $XRDEXIT
#fi

# echo "xrdcp input miniAOD"
# filein should already be in root://eoscms.cern.ch//store/mc/... format
# xrdcp -f ${FILEIN} delphesinput.root
# XRDEXIT=$?
# if [[ $XRDEXIT -ne 0 ]]; then
#     echo "exit code $XRDEXIT, failure in xrdcp of GEN file"
#     exit $XRDEXIT
# fi

setupTime=`date +%s`

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#run MiniAOD through Delphes

sed -i "s|MAXEVENTS|${MAXEVT}|g" cards/CMS_PhaseII/$detCard
sed -i "s|SKIPEVENTS|${SKIPEVT}|g" cards/CMS_PhaseII/$detCard

./DelphesCMSFWLite cards/CMS_PhaseII/$detCard ${FILEOUT} ${FILEIN}
DELPHESEXIT=$?
if [[ $DELPHESEXIT -ne 0 ]]; then
    echo "exit code $DELPHESEXIT, failure in DelphesCMSFWLite (maybe from xrootd)"
    exit $DELPHESEXIT
fi

DelphesTime=`date +%s`

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#metadata

echo "User: " `eval whoami`
echo "Date: " `date` 
echo 

echo "Process: " $process 
echo "Pileup Conditions: " $nPU 
echo "Configuration: " $configuration 
echo "Energy: " $energy 
echo 

echo "Input MiniAOD: " $FILEIN
echo "Skipped Events: " $SKIPEVT 
echo "Run Events: " $MAXEVT 
echo 

echo "Delphes Output: " $FILEOUT
echo "Delphes Version: " $DelphesVersion 
echo "Detector Card: " $detCard 
echo 
echo "Mass point: " $MASS


echo "Minutes spent setting up job: " `expr $setupTime / 60 - $startTime / 60` 
echo "Minutes spent running Delphes: " `expr $DelphesTime / 60 - $setupTime / 60` 
echo 

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# copy output to eos
# Running @CERN this should be fine with the root://eoscms.cern.ch eosOutDir
# Running @DESY likely need to change this copy command!
mkdir Output30PU-S13-M$MASS
mv ${FILEOUT} Output30PU-S13-M$MASS


XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
    echo "exit code $XRDEXIT, failure in xrdcp of ROOT"
    exit $XRDEXIT
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
endTime=`date +%s`
echo "Time spent copying output (s): " `expr $endTime - $DelphesTime`
echo "Total runtime (m): " `expr $endTime / 60 - $startTime / 60`



