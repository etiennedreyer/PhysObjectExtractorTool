#!/bin/bash
source /usr/wipp/conda/24.5.0u/etc/profile.d/conda.sh
topdir="/storage/agrp/dreyet/CMS_OpenData/FullEvent/CMSSW_5_3_32/src/PhysObjectExtractorTool"

cd ${topdir}

### setup
cd ../../../
conda activate apptainer
source /cvmfs/cms.cern.ch/cmsset_default.sh

### The commands need to be in quotes, but the quotes should not touch the commands. This is very finnicky!
# cmssw-el6 --command-to-run \" cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; source pbs/${tag}/run_cms_${jobID}.sh \"
### Edit (after upgrade to el9): backslashes are not needed anymore
cmssw-el6 --command-to-run " cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; source pbs/${tag}/run_cms_${jobID}.sh "