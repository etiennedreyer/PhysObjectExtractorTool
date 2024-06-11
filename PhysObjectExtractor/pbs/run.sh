#/bin/bash
source ~/.bashrc
topdir="/storage/agrp/dreyet/CMS_OpenData/FullEvent/CMSSW_5_3_32/src/PhysObjectExtractorTool"

cd ${topdir}
source setup.sh

cmsRun python/poet_cfg_genParticles.py ${infile} ${outfile} ${nevents}