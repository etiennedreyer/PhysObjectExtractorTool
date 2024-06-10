conda activate appenv

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmssw-el5
cmsrel CMSSW_5_3_11_patch6
cd CMSSW_5_3_11_patch6

cd src
cmsenv
cd ../../