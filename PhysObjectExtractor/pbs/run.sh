#/bin/bash
source ~/.bashrc
topdir="/storage/agrp/dreyet/CMS_OpenData/FullEvent/CMSSW_5_3_32/src/PhysObjectExtractorTool"

cd ${topdir}

### setup
cd ../../../
conda activate appenv
source /cvmfs/cms.cern.ch/cmsset_default.sh

### The commands need to be in quotes, but the quotes should not touch the commands. This is very finnicky!
cmssw-el6 --command-to-run \" cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; cmsRun python/poet_cfg_genParticles.py ${infile} ${outfile} ${nevents} \"