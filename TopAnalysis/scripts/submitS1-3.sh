#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

export OUTPUT=S13SingleAndPair-GenOnly-$1.root
export OUTPUTLHE=S13-LHE-$1.root

cd /nfs/dust/cms/user/gunnep/GridpackGen/CMSSW_9_3_9/src/
mkdir Output_inclusiveSingleAndPairS1-3-M1800_$1
cd Output_inclusiveSingleAndPairS1-3-M1800_$1
cp ../LQGenerationS1-3TopMu_cfg.py .
cmsRun LQGenerationS1-3TopMu_cfg.py

exit 0;
