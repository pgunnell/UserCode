# UserCode
# How to download the analyzer

cmsrel CMSSW_8_0_25

cd CMSSW_8_0_25/src

cmsenv

git clone -b delphes https://github.com/pgunnell/UserCode.git 

cd UserCode/Analysis/

scramv1 b -r -j8

cd plugins

#Modify .cc and .h and run with the ../test/*.py file
