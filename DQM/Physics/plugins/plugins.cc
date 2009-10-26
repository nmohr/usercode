#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DQM/Physics/src/BPhysicsOniaDQM.h"
#include "DQM/Physics/src/EwkDQM.h"
#include "DQM/Physics/src/QcdPhotonsDQM.h"
#include "DQM/Physics/src/QcdHighPtDQM.h"
#include "DQM/Physics/src/TopDiLeptonDQM.h"
#include "DQM/Physics/src/SusyLeptonDQM.h"
#include "DQM/Physics/src/SusyHadronDQM.h"

typedef SusyLeptonDQM< reco::Muon, reco::GsfElectron, reco::CaloJet, reco::CaloMET > RecoSusyLeptonDQM;
typedef SusyHadronDQM< reco::Muon, reco::GsfElectron, reco::CaloJet, reco::CaloMET > RecoSusyHadronDQM;

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(BPhysicsOniaDQM);
DEFINE_ANOTHER_FWK_MODULE(EwkDQM);
DEFINE_ANOTHER_FWK_MODULE(QcdPhotonsDQM);
DEFINE_ANOTHER_FWK_MODULE(QcdHighPtDQM);
DEFINE_ANOTHER_FWK_MODULE(TopDiLeptonDQM);
DEFINE_ANOTHER_FWK_MODULE(RecoSusyLeptonDQM);
DEFINE_ANOTHER_FWK_MODULE(RecoSusyHadronDQM);
