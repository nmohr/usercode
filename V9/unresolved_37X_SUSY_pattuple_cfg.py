#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 37X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV9
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.30 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
    annotation = cms.untracked.string('SUSY pattuple definition')
)

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
     '/store/data/Commissioning10/MinimumBias/RECO/SD_Mu-Jun9thSkim_v1/0041/FE82CCF6-197F-DF11-AEA5-001A92971B8A.root'
     #'/store/relval/CMSSW_3_7_0_pre5/RelValProdTTbar/GEN-SIM-RECO/MC_37Y_V4-v1/0023/BA92C6D3-8863-DF11-B3AF-002618943939.root'
     #'/store/relval/CMSSW_3_6_0_pre6/RelValProdTTbar/GEN-SIM-RECO/MC_36Y_V4-v1/0011/82DAA1BE-B344-DF11-A116-00304867C0C4.root'
     #'/store/relval/CMSSW_3_7_0_pre5/MinimumBias/RECO/GR_R_37X_V4_RelVal_col_10-v1/0024/CAD8C3ED-8E63-DF11-B60D-00261894393B.root'
     #'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0007/FE90A396-233C-DF11-8106-002618943898.root'
     #'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr1Skim_GOODCOLL-v1/0140/E27B88D1-8040-DF11-B3FC-00261894391B.root'
    ]
process.maxEvents.input = 100
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'GR_R_37X_V4::All'

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands
#Apply SUSYPAT, parameters are: mcInfo, HLT menu, Jet energy corrections, mcVersion ('35x' for 35x samples, empty string for 36X samples),JetCollections
addDefaultSUSYPAT(process,False,'HLT','Spring10','',['IC5Calo','AK5PF']) 
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
############################## END SUSYPAT specifics ####################################


#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )

#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.susyPatDefaultSequence )
#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
