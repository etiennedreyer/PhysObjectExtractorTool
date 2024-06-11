import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

import os
import sys

relBase = os.environ["CMSSW_BASE"]

batchMode = False

if len(sys.argv) > 3:
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    if len(sys.argv) > 4:
        num_events = int(sys.argv[4])
    else:
        num_events = -1
else:
    batchMode = True
    num_events = -1

process = cms.Process("POET")

# ---- Configure the framework messaging system
# ---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# ---- Select the maximum number of events to process (if -1, run over all events)
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(num_events))

# ---- Load needed configuration
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

if not batchMode:
    # ---- Define the test source files to be read using the xrootd protocol (root://), or local files (file:)
    # ---- Several files can be comma-separated
    # ---- A local file, for testing, can be downloaded using, e.g., the cern open data client (https://cernopendata-client.readthedocs.io/en/latest/):
    # ---- python cernopendata-client download-files --recid 6004 --filter-range 1-1
    # ---- For running over larger number of files, comment out this section and use/uncomment the FileUtils infrastructure below
    # sourceFile='file:/eos/opendata/cms/MonteCarlo2011/Summer11LegDR/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM/START53_LV4-v1/10000/00064CCC-A218-E311-A2E9-D485646A4E1A.root'
    sourceFile = "file:" + input_file
    process.source = cms.Source(
        "PoolSource",
        fileNames=cms.untracked.vstring(
            #'file:/playground/1EC938EF-ABEC-E211-94E0-90E6BA442F24.root'
            sourceFile
        ),
    )
else:
    # ---- Alternatively, to run on larger scale, one could use index files as obtained from the Cern Open Data Portal
    # ---- and pass them into the PoolSource.  The example is for 2012 data
    files = FileUtils.loadListFromFile(
        "CMS_MonteCarlo2011_Summer11LegDR_MinBias_TuneZ2_7TeV-pythia6_GEN-SIM_START53_LV4-v1_10000_file_index.txt"
    )
    process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*files))


# ---- These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
# ---- Comment theese lines for launching at WIS
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.GlobalTag.globaltag = "START53_LV6A1::All"



# --- Jet loading configuration 

#---- Get non-PAT access to the jet flavour information
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()
from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
process.jetFlavourInfosAK5PFJets = ak5JetFlavourInfos.clone()

#---- Configure the POET jet analyzer
#---- Don't forget to run jec_cfg.py to get these .txt files!
JecString = 'START53_LV6A1_'
process.pfJetsAk5= cms.EDAnalyzer('JetAnalyzer',
                    InputCollection = cms.InputTag("ak5PFJets"),
                    isData = cms.bool(False),
                    isSim  = cms.bool(False),
                    jecL1Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L1FastJet_AK5PF.txt'), 
                    jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2Relative_AK5PF.txt'),
                    jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L3Absolute_AK5PF.txt'),
                    jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'L2L3Residual_AK5PF.txt'),
                    jecUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/'+JecString+'Uncertainty_AK5PF.txt'),
                    jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/JetResolutionInputAK5PF.txt')
                    )

process.genJetsAk7 = cms.EDAnalyzer('GenJetAnalyzer',
                    InputCollection = cms.InputTag("ak7GenJets")
                    )

process.genJetsAk5 = cms.EDAnalyzer('GenJetAnalyzer',
                    InputCollection = cms.InputTag("ak5GenJets")
                    )
# ---- Configure the PhysObjectExtractor modules!

# ---- More information about InputCollections at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
process.events = cms.EDAnalyzer('EventAnalyzer')

process.gens = cms.EDAnalyzer(
    "GenParticleAnalyzer",
    # ---- Collect particles with specific "status:pdgid"
    # ---- Check PDG ID in the PDG.
    # ---- if 0:0, collect them all
    input_particle=cms.vstring("1:0"),
    InputJetCollection = cms.InputTag("ak5GenJets")
)

process.pfcs = cms.EDAnalyzer(
    "ParticleFlowAnalyzer",
    input_particle=cms.vstring("0:0"),
    InputJetCollection = cms.InputTag("ak5PFJets")
)

process.vtxs = cms.EDAnalyzer(
    "VertexAnalyzer"
)


# --- Trigger

process.trigger = cms.EDAnalyzer('TriggerAnalyzer',
                              processName = cms.string("HLT"),
                              #---- These are example triggers for 2011 DoubleMu dataset
                              #---- Wildcards * and ? are accepted (with usual meanings)
                               #---- If left empty, all triggers will run              
#                              triggerPatterns = cms.vstring("HLT_L2DoubleMu23_NoVertex_v*","HLT_Mu13_Mu8_v*", "HLT_DoubleMu45_v*", "HLT_Mu8_Jet40_v*", "HLT_TripleMu5_v*"), 
                            #   triggerPatterns = cms.vstring("HLT_L2DoubleMu23_NoVertex_v*","HLT_Mu13_Mu8_v*"),
                              triggerPatterns = cms.vstring("HLT_Jet300_*"),
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
                              )


# --- PileUp
process.pu = cms.EDAnalyzer('PileupEventAnalyzer')

# ---- Configure the output ROOT filename
process.TFileService = cms.Service("TFileService", fileName=cms.string(output_file))

# ---- Finally run everything!
# ---- Separation by * implies that processing order is important.
# ---- separation by + implies that any order will work
# ---- One can put in or take out the needed processes
process.p = cms.Path(process.events+process.gens+process.pfcs+process.pfJetsAk5+process.genJetsAk5+process.genJetsAk7+process.vtxs+process.trigger + process.pu)
# process.p = cms.Path(process.events+process.gens+process.pfcs+process.pfJetsAk5+process.genJetsAk5+process.genJetsAk7)
# process.p = cms.Path(process.genparticles)