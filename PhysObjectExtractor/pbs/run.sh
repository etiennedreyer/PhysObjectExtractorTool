#!/bin/bash
export IOTHROTTLE_LIMIT=5
source /usr/wipp/conda/24.5.0u/etc/profile.d/conda.sh
topdir="/srv01/agrp/dmitrykl/projects/cmssw/"

cd ${topdir}

### setup
conda activate apptainer
source /cvmfs/cms.cern.ch/cmsset_default.sh
unset LD_PRELOAD
### The commands need to be in quotes, but the quotes should not touch the commands. This is very finnicky!
# cmssw-el6 --command-to-run \" cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; source pbs/${tag}/run_cms_${jobID}.sh \"
### Edit (after upgrade to el9): backslashes are not needed anymore
cmssw-el6 -B /storage/agrp/dmitrykl -B /run --command-to-run " cd ~/projects/cmssw; cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; source pbs/${tag}/run_cms_${jobID}.sh "
# cmssw-el6 --command-to-run " cmsrel CMSSW_5_3_32; cd CMSSW_5_3_32/src; cmsenv; cd PhysObjectExtractorTool/PhysObjectExtractor/; source pbs/${tag}/run_cms_${jobID}.sh "