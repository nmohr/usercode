import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')
process.load('SUSYBSMAnalysis.SusyCAF.ntuple_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
                           fileName = cms.string( 'testTree.root' ),
                           closeFileFast = cms.untracked.bool(True)  ) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/disk1/user/nmohr/E81DF107-84AB-DE11-A089-001E4F3F28DC.root'),
    secondaryFileNames = cms.untracked.vstring())
