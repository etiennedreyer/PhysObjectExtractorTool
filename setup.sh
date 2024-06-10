cd ../../../

conda activate appenv

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmssw-el6
cmsrel CMSSW_5_3_32

cd CMSSW_5_3_32/src
cmsenv

cd PhysObjectExtractorTool/PhysObjectExtractor/
scram b