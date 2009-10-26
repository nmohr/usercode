#ifndef SusyHadronDQM_H
#define SusyHadronDQM_H

#include <string>
#include <vector>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class TH1F;
class TH2F;

class PtGreater {
   public:
      template<typename T> bool operator ()(const T& i, const T& j) {
         return (i.pt() > j.pt());
      }
};

template<typename Mu, typename Ele, typename Jet, typename Met>
class SusyHadronDQM: public edm::EDAnalyzer {

   public:

      explicit SusyHadronDQM(const edm::ParameterSet&);
      ~SusyHadronDQM();

   protected:

      void beginRun(const edm::Run&);
      void endRun(const edm::Run&);

   private:

      void initialize();
      virtual void beginJob(const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual bool goodSusyElectron(const Ele*);
      virtual bool goodSusyMuon(const Mu*);
      virtual void endJob();

      edm::ParameterSet parameters_;
      DQMStore * dbe_;

      std::string moduleName_;

      edm::InputTag muons_;
      double muon_pt_cut_;
      double muon_eta_cut_;
      double muon_nHits_cut_;
      double muon_nChi2_cut_;
      double muon_d0_cut_;
      double muon_iso_cut_;

      edm::InputTag electrons_;
      double elec_pt_cut_;
      double elec_eta_cut_;
      double elec_mva_cut_;
      double elec_d0_cut_;
      double elec_iso_cut_;

      edm::InputTag jets_;
      double jet_pt_cut_;
      double jet1_pt_cut_;
      double jet2_pt_cut_;
      double jet3_pt_cut_;
      double jet_eta_cut_;
      double jet_min_emf_cut_;
      double jet_max_emf_cut_;
      int N_jets_cut_;

      edm::InputTag met_;
      double met_cut_;
      double mht_cut_;
      double ht_cut_;
      double deltaPhi_cut_;

      edm::InputTag vertex_;
      math::XYZPoint bs;

      MonitorElement * h_N_jets_;
      MonitorElement * h_ht_;
      MonitorElement * h_pt_jet1_;
      MonitorElement * h_pt_jet2_;
      MonitorElement * h_deltaPhi_mht_jets_;
      MonitorElement * h_mht_;
      MonitorElement * h_pt_muons_;
      MonitorElement * h_pt_elecs_;

};

template<typename Mu, typename Ele, typename Jet, typename Met>
SusyHadronDQM<Mu, Ele, Jet, Met>::SusyHadronDQM(const edm::ParameterSet& pset) {

   parameters_ = pset;
   initialize();

   moduleName_ = pset.getUntrackedParameter<std::string> ("moduleName");

   muons_ = pset.getParameter<edm::InputTag> ("muonCollection");
   muon_pt_cut_ = pset.getParameter<double> ("muon_pt_cut");
   muon_eta_cut_ = pset.getParameter<double> ("muon_eta_cut");
   muon_nHits_cut_ = pset.getParameter<double> ("muon_nHits_cut");
   muon_nChi2_cut_ = pset.getParameter<double> ("muon_nChi2_cut");
   muon_d0_cut_ = pset.getParameter<double> ("muon_d0_cut");
   muon_iso_cut_ = pset.getParameter<double> ("muon_iso_cut");

   electrons_ = pset.getParameter<edm::InputTag> ("electronCollection");
   elec_pt_cut_ = pset.getParameter<double> ("elec_pt_cut");
   elec_eta_cut_ = pset.getParameter<double> ("elec_eta_cut");
   elec_mva_cut_ = pset.getParameter<double> ("elec_mva_cut");
   elec_d0_cut_ = pset.getParameter<double> ("elec_d0_cut");
   elec_iso_cut_ = pset.getParameter<double> ("elec_iso_cut");

   jets_ = pset.getParameter<edm::InputTag> ("jetCollection");
   jet_pt_cut_ = pset.getParameter<double> ("jet_pt_cut");
   jet1_pt_cut_ = pset.getParameter<double> ("jet1_pt_cut");
   jet2_pt_cut_ = pset.getParameter<double> ("jet2_pt_cut");
   jet3_pt_cut_ = pset.getParameter<double> ("jet3_pt_cut");
   jet_eta_cut_ = pset.getParameter<double> ("jet_eta_cut");
   jet_min_emf_cut_ = pset.getParameter<double> ("jet_min_emf_cut");
   jet_max_emf_cut_ = pset.getParameter<double> ("jet_max_emf_cut");
   N_jets_cut_ = pset.getParameter<int> ("N_jets_cut");

   met_ = pset.getParameter<edm::InputTag> ("metCollection");
   ht_cut_ = pset.getParameter<double> ("ht_cut");
   met_cut_ = pset.getParameter<double> ("met_cut");
   mht_cut_ = pset.getParameter<double> ("mht_cut");
   deltaPhi_cut_ = pset.getParameter<double> ("deltaPhi_cut");

   vertex_ = pset.getParameter<edm::InputTag> ("vertexCollection");

}

template<typename Mu, typename Ele, typename Jet, typename Met>
SusyHadronDQM<Mu, Ele, Jet, Met>::~SusyHadronDQM() {

}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::initialize() {

}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::beginJob(const edm::EventSetup& evt) {

   dbe_ = edm::Service<DQMStore>().operator->();

   dbe_->setCurrentFolder(moduleName_);

   h_N_jets_ = dbe_->book1D("N_jets", "N_jets", 10, 0., 10.);
   h_ht_ = dbe_->book1D("ht", "ht", 100, 0., 2000.);
   h_pt_jet1_ = dbe_->book1D("pt_jet1", "pt_jet1", 100, 0., 1000.);
   h_pt_jet2_ = dbe_->book1D("pt_jet2", "pt_jet2", 100, 0., 1000.);
   h_deltaPhi_mht_jets_ = dbe_->book1D("deltaPhi_mht_jets", "deltaPhi_mht_jets", 100, 0., 2.);
   h_mht_ = dbe_->book1D("mht", "mht", 100, 0., 2000.);
   h_pt_muons_ = dbe_->book1D("pt_muons", "pt_muons", 100, 0., 200.);
   h_pt_elecs_ = dbe_->book1D("pt_elecs", "pt_elecs", 100, 0., 200.);

}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::beginRun(const edm::Run& run) {

}

template<typename Mu, typename Ele, typename Jet, typename Met>
bool SusyHadronDQM<Mu, Ele, Jet, Met>::goodSusyElectron(const Ele* ele) {
   //   if (ele->pt() < elec_pt_cut_)
   //      return false;
   if (fabs(ele->eta()) > elec_eta_cut_)
      return false;
   if (ele->mva() < elec_mva_cut_)
      return false;
   if (fabs(ele->gsfTrack()->dxy(bs)) > elec_d0_cut_)
      return false;
   return true;
}

template<typename Mu, typename Ele, typename Jet, typename Met>
bool SusyHadronDQM<Mu, Ele, Jet, Met>::goodSusyMuon(const Mu* mu) {
   //   if (mu->pt() < muon_pt_cut_)
   //      return false;
   if (fabs(mu->eta()) > muon_eta_cut_)
      return false;
   if (!mu->isGlobalMuon())
      return false;
   if (mu->innerTrack()->numberOfValidHits() < muon_nHits_cut_)
      return false;
   if (mu->globalTrack()->normalizedChi2() > muon_nChi2_cut_)
      return false;
   if (fabs(mu->innerTrack()->dxy(bs)) > muon_d0_cut_)
      return false;
   return true;
}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::analyze(const edm::Event& evt, const edm::EventSetup& iSetup) {

   edm::Handle<std::vector<Mu> > muons;
   evt.getByLabel(muons_, muons);

   edm::Handle<std::vector<Ele> > elecs;
   evt.getByLabel(electrons_, elecs);

   edm::Handle<std::vector<Jet> > jets;
   evt.getByLabel(jets_, jets);
   //// sorted jets
   //   edm::Handle<std::vector<Jet> > cJets;
   //   evt.getByLabel(jets_, cJets);
   //   std::vector<Jet> jets = *cJets;
   //   std::sort(jets.begin(), jets.end(), PtGreater());

   edm::Handle<std::vector<Met> > mets;
   evt.getByLabel(met_, mets);

   edm::Handle<reco::VertexCollection> vertices;
   evt.getByLabel(vertex_, vertices);

   float HT = 0.;
   math::PtEtaPhiMLorentzVector vMHT(0., 0., 0., 0.);
   int nJets = 0;
   float jet1_pt = 0;
   float jet2_pt = 0;
   for (typename std::vector<Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i) {
      if (jet_i->pt() < jet_pt_cut_)
         continue;
      if (jet_i->eta() > fabs(jet_eta_cut_))
         continue;
      if (jet_i->emEnergyFraction() < (jet_min_emf_cut_))
         continue;
      if (jet_i->emEnergyFraction() > (jet_max_emf_cut_))
         continue;
      if (nJets == 0)
         jet1_pt = jet_i->pt();
      if (nJets == 1)
         jet2_pt = jet_i->pt();
      ++nJets;
      HT += jet_i->pt();
      vMHT -= jet_i->p4();
   }
   float MHT = vMHT.pt();

   int i = 0;
   float minDeltaPhi = 9999.;
   for (typename std::vector<Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i) {
      ++i;
      if (i <= 3) {
         double deltaPhi = reco::deltaPhi(jet_i->phi(),vMHT.phi());
         if (fabs(deltaPhi) < minDeltaPhi)
            minDeltaPhi = fabs(deltaPhi);
      }
   }

   float MET = 0.;
   for (typename std::vector<Met>::const_iterator met_i = mets->begin(); met_i != mets->end(); ++met_i) {
      MET = met_i->pt();
      break;
   }

   for (reco::VertexCollection::const_iterator vertex = vertices->begin(); vertex != vertices->end(); ++vertex) {
      bs = vertex->position();
      break;
   }

   float leadingMuPt = 0;
   for (typename std::vector<Mu>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i) {
      if (!goodSusyMuon(&(*mu_i)))
         continue;

      reco::MuonIsolation Iso_muon = mu_i->isolationR03();
      float muIso = Iso_muon.emEt + Iso_muon.hadEt + Iso_muon.sumPt;

      if (muIso < muon_iso_cut_) {
         if (mu_i->pt() > leadingMuPt)
            leadingMuPt = mu_i->pt();
      }
   }

   float leadingElecPt = 0;
   for (typename std::vector<Ele>::const_iterator ele_i = elecs->begin(); ele_i != elecs->end(); ++ele_i) {
      if (!goodSusyElectron(&(*ele_i)))
         continue;

      float elecIso = ele_i->dr03TkSumPt() + ele_i->dr03EcalRecHitSumEt() + ele_i->dr03HcalTowerSumEt();

      if (elecIso < elec_iso_cut_) {
         if (ele_i->pt() > leadingElecPt)
            leadingElecPt = ele_i->pt();
      }
   }

   //// Fill N-1 hsitograms
   for (int i = 0; i < 8; ++i) {
      if (nJets >= N_jets_cut_ || i == 0) {
         if (HT >= ht_cut_ || i == 1) {
            if (jet1_pt >= jet1_pt_cut_ || i == 2) {
               if (jet2_pt >= jet2_pt_cut_ || i == 3) {
                  if (minDeltaPhi >= deltaPhi_cut_ || i == 4) {
                     if (MHT >= mht_cut_ || i == 5) {
                        if (leadingMuPt <= muon_pt_cut_ || i == 6) {
                           if (leadingElecPt <= elec_pt_cut_ || i == 7) {
                              if (i == 0)
                                 h_N_jets_->Fill(nJets);
                              if (i == 1)
                                 h_ht_->Fill(HT);
                              if (i == 2)
                                 h_pt_jet1_->Fill(jet1_pt);
                              if (i == 3)
                                 h_pt_jet2_->Fill(jet2_pt);
                              if (i == 4)
                                 h_deltaPhi_mht_jets_->Fill(minDeltaPhi);
                              if (i == 5)
                                 h_mht_->Fill(MHT);
                              if (i == 6)
                                 h_pt_muons_->Fill(leadingMuPt);
                              if (i == 7)
                                 h_pt_elecs_->Fill(leadingElecPt);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::endRun(const edm::Run& run) {

}

template<typename Mu, typename Ele, typename Jet, typename Met>
void SusyHadronDQM<Mu, Ele, Jet, Met>::endJob() {

}

#endif
