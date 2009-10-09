#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

DEFINE_SEAL_MODULE();

#include "SUSYBSMAnalysis/SusyCAF/interface/SusyTree.h"
#include "SUSYBSMAnalysis/SusyCAF/interface/SusyCAF_Event.h"
#include "SUSYBSMAnalysis/SusyCAF/interface/SusyCAF_MET.h"
#include "SUSYBSMAnalysis/SusyCAF/interface/SusyCAF_Electron.h"
#include "SUSYBSMAnalysis/SusyCAF/interface/SusyCAF_Muon.h"

typedef SusyCAF_MET<reco::CaloMET> SusyCAF_CaloMET;
typedef SusyCAF_Electron<reco::GsfElectron> SusyCAF_GsfElectron;
typedef SusyCAF_Muon<reco::Muon> SusyCAF_RecoMuon;

DEFINE_ANOTHER_FWK_MODULE(SusyTree);
DEFINE_ANOTHER_FWK_MODULE(SusyCAF_Event);
DEFINE_ANOTHER_FWK_MODULE(SusyCAF_CaloMET);
DEFINE_ANOTHER_FWK_MODULE(SusyCAF_GsfElectron);
DEFINE_ANOTHER_FWK_MODULE(SusyCAF_RecoMuon);
