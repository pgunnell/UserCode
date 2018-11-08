UserCode
========

mkdir MyDirectory

cd MyDirectory

module use -a /afs/desy.de/group/cms/modulefiles/

module load cmssw

cmsrel CMSSW_9_4_6

cd CMSSW_9_4_6/src

cmsenv

git init

git clone -b GenLevel https://github.com/pgunnell/UserCode.git

For compiling (needed every time you make changes in a .cc/.h file):

scramv1 b -r -j8
