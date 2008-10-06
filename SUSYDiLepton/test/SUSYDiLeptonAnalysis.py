import FWCore.ParameterSet.Config as cms

process = cms.Process('Analysis')
process.options = cms.untracked.PSet(
     wantSummary = cms.untracked.bool(True)
)
process.MessageLogger = cms.Service('MessageLogger',
  test_baseCuts_LM9t175_sftsdkpyt = cms.untracked.PSet(
     INFO = cms.untracked.PSet(
          limit = cms.untracked.int32(0)
     ),
     FwkReport = cms.untracked.PSet(
          reportEvery = cms.untracked.int32(1000),
          limit = cms.untracked.int32(1000)
     ),
     default = cms.untracked.PSet(
          limit = cms.untracked.int32(100)
     ),
     Root_NoDictionary = cms.untracked.PSet(
          limit = cms.untracked.int32(0)
     ),
     FwkJob = cms.untracked.PSet(
          limit = cms.untracked.int32(0)
     ),
     FwkSummary = cms.untracked.PSet(
          reportEvery = cms.untracked.int32(1),
          limit = cms.untracked.int32(10000000)
     ),
     threshold = cms.untracked.string('INFO')
  ),
  destinations = cms.untracked.vstring('test_baseCuts_LM9t175_sftsdkpyt')
)
process.source = cms.Source('PoolSource', 
     fileNames = cms.untracked.vstring('file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_7.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_1.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_2.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_6.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_3.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_4.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_8.root','file:/opt/user/nmohr/skimCSA07/data/LM9t175_sftsdkpyt/skimPATlepHLT_5.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.TFileService = cms.Service('TFileService', fileName = cms.string('test.baseCuts.LM9t175_sftsdkpyt.root'))

process.DiLeptonAnalysis = cms.EDAnalyzer('SUSYDiLeptonAnalysis',

debug = cms.untracked.bool(False),
mcInfo = cms.untracked.bool(True),

mcSource = cms.InputTag("genParticles"),
backmapSource = cms.InputTag("genParticles"),
muonSource = cms.InputTag("selectedLayer1Muons"),
electronSource = cms.InputTag("selectedLayer1Electrons"),
metSource = cms.InputTag("selectedLayer1METs"),
jetSource = cms.InputTag("selectedLayer1Jets"),

CSA_weighted = cms.untracked.bool(False),

rej_Cuts = cms.untracked.bool(False), 
rej_JetCut = cms.untracked.bool(True),
rej_METCut = cms.untracked.bool(True), 
rej_TwoLeptonCut = cms.untracked.bool(True),
rej_bTagCut = cms.untracked.bool(False),
user_nJets = cms.untracked.int32(4), #if 0 sum of 4Jets is checked
user_nbJets = cms.untracked.int32(2), #exactly the number of b_jets
user_pt1JetMin = cms.untracked.double(120.),
user_pt2JetMin = cms.untracked.double(80.),
user_pt3JetMin = cms.untracked.double(80.),
user_pt4JetMin = cms.untracked.double(80.),
user_METMin = cms.untracked.double(50.),
user_bTagDiscriminator = cms.untracked.double(0.4),

acc_MuonPt = cms.untracked.double(10.), 
acc_MuonEta = cms.untracked.double(2.), 
iso_MuonIso = cms.untracked.double(1.2),
    
acc_ElectronPt = cms.untracked.double(10.), 
acc_ElectronEta = cms.untracked.double(2.), 
iso_ElectronIso = cms.untracked.double(1.2), 
    
acc_JetPt = cms.untracked.double(30.), 
acc_JetEta = cms.untracked.double(2.5) 
)

process.p = cms.Path(process.DiLeptonAnalysis)
