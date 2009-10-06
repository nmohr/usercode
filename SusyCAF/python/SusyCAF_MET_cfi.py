import FWCore.ParameterSet.Config as cms

susycafmet = cms.EDProducer("SusyCAF_CaloMET",
                            InputTag = cms.InputTag('met'),
                            Prefix = cms.string('met'),
                            Suffix = cms.string('calo')
                            )
