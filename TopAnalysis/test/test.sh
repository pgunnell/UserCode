#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cmsRun Analysis_Template_MCFastSim.py

exit 0;
