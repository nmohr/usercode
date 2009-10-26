import FWCore.ParameterSet.Config as cms

susyHadronDQM = cms.EDAnalyzer("RecoSusyHadronDQM",

    moduleName     = cms.untracked.string('Physics/Susy/Hadron'),

    muonCollection = cms.InputTag('muons'),
    muon_pt_cut    = cms.double(  7.0 ),
    muon_eta_cut   = cms.double(  2.5 ),
    muon_nHits_cut = cms.double(  11 ),
    muon_nChi2_cut = cms.double(  10 ),
    muon_d0_cut    = cms.double(  0.2 ),
    muon_iso_cut   = cms.double(  5. ),

    electronCollection = cms.InputTag('gsfElectrons'),
    elec_pt_cut    = cms.double(  7.0 ),
    elec_eta_cut   = cms.double(  2.5 ),
    elec_mva_cut   = cms.double(  0.1 ),
    elec_d0_cut    = cms.double(  0.2 ),
    elec_iso_cut   = cms.double(  5. ),
    
    jetCollection   = cms.InputTag('antikt5CaloJets'),
    jet_pt_cut      = cms.double(  30.0 ),
    jet_eta_cut     = cms.double(   3.0 ),
    jet_min_emf_cut = cms.double(   0.05 ),
    jet_max_emf_cut = cms.double(   0.95 ),
    jet1_pt_cut     = cms.double( 110.0 ),
    jet2_pt_cut     = cms.double(  80.0 ),
    jet3_pt_cut     = cms.double(  30.0 ),
    N_jets_cut      = cms.int32(    3 ),
    
    metCollection = cms.InputTag('met'),
    met_cut       = cms.double( 100.0 ),
    ht_cut        = cms.double( 300.0 ),
    mht_cut       = cms.double( 100.0 ),
    deltaPhi_cut  = cms.double(   0.3 ),
    
    vertexCollection = cms.InputTag('offlinePrimaryVertices')

)

susyHadronAnalyzer = cms.Sequence(susyHadronDQM)
