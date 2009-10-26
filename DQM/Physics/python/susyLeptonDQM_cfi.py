import FWCore.ParameterSet.Config as cms

susyLeptonDQM = cms.EDAnalyzer("RecoSusyLeptonDQM",

    moduleName     = cms.untracked.string('Physics/Susy/Lepton'),

    muonCollection = cms.InputTag('muons'),
    muon_pt_cut    = cms.double(  7.0 ),
    muon_eta_cut   = cms.double(  2.5 ),
    muon_nHits_cut   = cms.double(  11 ),
    muon_nChi2_cut   = cms.double(  10 ),
    muon_d0_cut   = cms.double(  0.2 ),

    electronCollection = cms.InputTag('gsfElectrons'),
    elec_pt_cut    = cms.double(  7.0 ),
    elec_eta_cut   = cms.double(  2.5 ),
    elec_mva_cut   = cms.double(  0.1 ),
    elec_d0_cut   = cms.double(  0.2 ),
    
    jetCollection = cms.InputTag('antikt5CaloJets'),
    jet_pt_cut    = cms.double(  30.0 ),
    jet_eta_cut    = cms.double(  3.0 ),
    jet_min_emf_cut    = cms.double(  0.05 ),
    jet_max_emf_cut    = cms.double(  0.95 ),
    jet_sum_pt_cut    = cms.double(  100.0 ),
    
    metCollection = cms.InputTag('met'),
    met_cut    = cms.double(  50.0 ),
    
    vertexCollection = cms.InputTag('offlinePrimaryVertices')

)

susyLeptonAnalyzer = cms.Sequence(susyLeptonDQM)
