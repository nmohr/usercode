import FWCore.ParameterSet.Config as cms

susycafelectron = cms.EDProducer("SusyCAF_GsfElectron",
                            InputTag = cms.InputTag('gsfElectrons'),
                            Prefix = cms.string('electron'),
                            Suffix = cms.string('gsf')
                            )
