/*
 *  $Date: 2009/09/16 21:10:31 $
 *  $Revision: 1.1 $
 *  \author N. Mohr - niklas.mohr@cern.ch
 */


#ifndef SusyLeptonDQM_H
#define SusyLeptonDQM_H

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
#include "DataFormats/Math/interface/deltaR.h"

class TH1F;
class TH2F;

template< typename Mu, typename Ele, typename Jet, typename Met >
class SusyLeptonDQM : public edm::EDAnalyzer {

  public:

    explicit SusyLeptonDQM(const edm::ParameterSet&);
    ~SusyLeptonDQM();

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

    edm::InputTag electrons_;
    double elec_pt_cut_;
    double elec_eta_cut_;
    double elec_mva_cut_;
    double elec_d0_cut_;
    
    edm::InputTag jets_;
    double jet_pt_cut_;
    double jet_eta_cut_;
    double jet_min_emf_cut_;
    double jet_max_emf_cut_;
    double jet_sum_pt_cut_;
    
    edm::InputTag met_;
    double met_cut_;

    edm::InputTag vertex_;
    math::XYZPoint bs;

    MonitorElement * N_muons_;
    MonitorElement * pt_muons_;
    MonitorElement * eta_muons_;
    MonitorElement * phi_muons_;
    MonitorElement * Iso_muons_;

    MonitorElement * N_elecs_;
    MonitorElement * pt_elecs_;
    MonitorElement * eta_elecs_;
    MonitorElement * phi_elecs_;
    MonitorElement * Iso_elecs_;
    
    MonitorElement * Sum_pt_jets_;
    MonitorElement * Met_;
    
    MonitorElement * dR_emu_;
    
    MonitorElement * mass_OS_mumu_;
    MonitorElement * mass_OS_ee_;
    MonitorElement * mass_OS_emu_;
    MonitorElement * mass_SS_mumu_;
    MonitorElement * mass_SS_ee_;
    MonitorElement * mass_SS_emu_;

    MonitorElement * Muon_monitor_;
    MonitorElement * Electron_monitor_;
    MonitorElement * OSee_monitor_;
    MonitorElement * OSemu_monitor_;
    MonitorElement * OSmumu_monitor_;
    MonitorElement * SSee_monitor_;
    MonitorElement * SSemu_monitor_;
    MonitorElement * SSmumu_monitor_;
    MonitorElement * TriMuon_monitor_;
};

template< typename Mu, typename Ele, typename Jet, typename Met >
SusyLeptonDQM<Mu,Ele,Jet,Met>::SusyLeptonDQM(const edm::ParameterSet& pset) {

    parameters_ = pset;
    initialize();

    moduleName_     = pset.getUntrackedParameter<std::string>("moduleName");

    muons_          = pset.getParameter<edm::InputTag>("muonCollection");
    muon_pt_cut_    = pset.getParameter<double>("muon_pt_cut");
    muon_eta_cut_   = pset.getParameter<double>("muon_eta_cut");
    muon_nHits_cut_ = pset.getParameter<double>("muon_nHits_cut");
    muon_nChi2_cut_ = pset.getParameter<double>("muon_nChi2_cut");
    muon_d0_cut_    = pset.getParameter<double>("muon_d0_cut");

    electrons_      = pset.getParameter<edm::InputTag>("electronCollection");
    elec_pt_cut_    = pset.getParameter<double>("elec_pt_cut");
    elec_eta_cut_   = pset.getParameter<double>("elec_eta_cut");
    elec_mva_cut_   = pset.getParameter<double>("elec_mva_cut");
    elec_d0_cut_    = pset.getParameter<double>("elec_d0_cut");
  
    jets_           = pset.getParameter<edm::InputTag>("jetCollection");
    jet_pt_cut_     = pset.getParameter<double>("jet_pt_cut");
    jet_eta_cut_    = pset.getParameter<double>("jet_eta_cut");
    jet_min_emf_cut_= pset.getParameter<double>("jet_min_emf_cut");
    jet_max_emf_cut_= pset.getParameter<double>("jet_max_emf_cut");
    jet_sum_pt_cut_ = pset.getParameter<double>("jet_sum_pt_cut");
  
    met_            = pset.getParameter<edm::InputTag>("metCollection");
    met_cut_        = pset.getParameter<double>("met_cut");
    
    vertex_         = pset.getParameter<edm::InputTag>("vertexCollection");

}


template< typename Mu, typename Ele, typename Jet, typename Met >
SusyLeptonDQM<Mu,Ele,Jet,Met>::~SusyLeptonDQM() {

}


template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::initialize() {

}


template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::beginJob(const edm::EventSetup& evt) {

    dbe_ = edm::Service<DQMStore>().operator->();

    dbe_->setCurrentFolder(moduleName_);

    N_muons_    = dbe_->book1D("N_muons",    "N_muons",           10,  0.,  10.);
    pt_muons_   = dbe_->book1D("pt_muons",   "pt_muons",          300,  0., 300.);
    eta_muons_  = dbe_->book1D("eta_muons",  "eta_muons",         100, -2.5,   2.5);
    phi_muons_  = dbe_->book1D("phi_muons",  "phi_muons",         100, -4.,   4.);
    Iso_muons_  = dbe_->book1D("Iso_muons",  "Iso_muons",         100, 0., 25.);

    N_elecs_    = dbe_->book1D("N_elecs",    "N_elecs",           10,  0.,  10.);
    pt_elecs_   = dbe_->book1D("pt_elecs",   "pt_elecs",          300,  0., 300.);
    eta_elecs_  = dbe_->book1D("eta_elecs",  "eta_elecs",         100, -2.5, 2.5);
    phi_elecs_  = dbe_->book1D("phi_elecs",  "phi_elecs",         100, -4.,   4.);
    Iso_elecs_  = dbe_->book1D("Iso_elecs",  "Iso_elecs",         100, 0., 25.);
  
    Sum_pt_jets_  = dbe_->book1D("Sum_pt_jets",   "Sum_pt_jets",          200,  20., 2020.);
    Met_          = dbe_->book1D("Met",   "Met",          200,  0., 1000.);
  
    dR_emu_       = dbe_->book1D("deltaR_emu",   "deltaR_emu",          100,  0., 10.);
    
    mass_OS_mumu_ = dbe_->book1D("mass_OS_mumu",     "mass_OS_mumu",     100, 0., 300.);
    mass_OS_ee_ = dbe_->book1D("mass_OS_ee",     "mass_OS_ee",     100, 0., 300.);
    mass_OS_emu_ = dbe_->book1D("mass_OS_emu",     "mass_OS_emu",     100, 0., 300.);
    mass_SS_mumu_ = dbe_->book1D("mass_SS_mumu",     "mass_SS_mumu",     100, 0., 300.);
    mass_SS_ee_ = dbe_->book1D("mass_SS_ee",     "mass_SS_ee",     100, 0., 300.);
    mass_SS_emu_ = dbe_->book1D("mass_SS_emu",     "mass_SS_emu",     100, 0., 300.);

    Muon_monitor_ = dbe_->book2D("Single_Muon_Selection", "Single_Muon_Selection", 100, 0., 1000., 100, 0., 1000.);
    Electron_monitor_ = dbe_->book2D("Single_Electron_Selection", "Single_Electron_Selection", 100, 0., 1000., 100, 0., 1000.);
    OSee_monitor_ = dbe_->book2D("OS_Electron_Selection", "OS_Electron_Selection", 100, 0., 1000., 100, 0., 1000.);
    OSemu_monitor_ = dbe_->book2D("OS_ElectronMuon_Selection", "OS_ElectronMuon_Selection", 100, 0., 1000., 100, 0., 1000.);
    OSmumu_monitor_ = dbe_->book2D("OS_Muon_Selection", "OS_Muon_Selection", 100, 0., 1000., 100, 0., 1000.);
    SSee_monitor_ = dbe_->book2D("SS_Electron_Selection", "SS_Electron_Selection", 100, 0., 1000., 100, 0., 1000.);
    SSemu_monitor_ = dbe_->book2D("SS_ElectronMuon_Selection", "SS_ElectronMuon_Selection", 100, 0., 1000., 100, 0., 1000.);
    SSmumu_monitor_ = dbe_->book2D("SS_Muon_Selection", "SS_Muon_Selection", 100, 0., 1000., 100, 0., 1000.);
    TriMuon_monitor_ = dbe_->book2D("Tri_Muon_Selection", "Tri_Muon_Selection", 100, 0., 1000., 100, 0., 1000.);
}
    
template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::beginRun(const edm::Run& run) {

}

template< typename Mu, typename Ele, typename Jet, typename Met >
bool SusyLeptonDQM<Mu,Ele,Jet,Met>::goodSusyElectron(const Ele* ele){
    if( ele->pt()         < elec_pt_cut_  )  return false;
    if( fabs(ele->eta())  > elec_eta_cut_ )  return false;
    if( ele->mva()  < elec_mva_cut_ )  return false;
    if( fabs(ele->gsfTrack()->dxy(bs)) > elec_d0_cut_ ) return false;
    return true;
}

template< typename Mu, typename Ele, typename Jet, typename Met >
bool SusyLeptonDQM<Mu,Ele,Jet,Met>::goodSusyMuon(const Mu* mu){
    if( mu->pt()        < muon_pt_cut_  )  return false;
    if( fabs(mu->eta()) > muon_eta_cut_ )  return false;
    if( !mu->isGlobalMuon() )  return false;
    if( mu->innerTrack()->numberOfValidHits()   < muon_nHits_cut_ )  return false;
    if( mu->globalTrack()->normalizedChi2()     > muon_nChi2_cut_ )  return false;
    if( fabs(mu->innerTrack()->dxy(bs)) > muon_d0_cut_ ) return false;
    return true;
}

template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::analyze(const edm::Event& evt, const edm::EventSetup& iSetup) {

    edm::Handle<std::vector<Mu> > muons;
    evt.getByLabel(muons_, muons);

    edm::Handle<std::vector<Ele> > elecs;
    evt.getByLabel(electrons_, elecs);
  
    edm::Handle<std::vector<Jet> > jets;
    evt.getByLabel(jets_, jets);
  
    edm::Handle<std::vector<Met> > mets;
    evt.getByLabel(met_, mets);

    edm::Handle<reco::VertexCollection> vertices;
    evt.getByLabel(vertex_, vertices);
    
 
    float sumPt = 0.; 
    for(typename std::vector<Jet>::const_iterator jet_i = jets->begin(); jet_i!= jets->end(); ++jet_i) {
        if (jet_i->pt() < jet_pt_cut_ ) continue;
        if (jet_i->eta() > fabs(jet_eta_cut_) ) continue;
        if (jet_i->eta() > fabs(jet_eta_cut_) ) continue;
        if (jet_i->emEnergyFraction() < (jet_min_emf_cut_) ) continue;
        if (jet_i->emEnergyFraction() > (jet_max_emf_cut_) ) continue;
        sumPt += jet_i->pt(); 
    }

    Sum_pt_jets_->Fill(sumPt);
  
    float MET = 0.; 
    for(typename std::vector<Met>::const_iterator met_i = mets->begin(); met_i!= mets->end(); ++met_i) {
        MET = met_i->pt();
        break; 
    }
  
    Met_->Fill(MET);

    for (reco::VertexCollection::const_iterator vertex = vertices->begin(); vertex != vertices->end(); ++vertex) {
        bs = vertex->position();
        break;
    }

    int nMuons = 0;
    int nSSmumu = 0;
    int nOSmumu = 0;
    int nSSemu = 0;
    int nOSemu = 0;
    float inv = 0.;
    float dR = 0.;

    for(typename std::vector<Mu>::const_iterator mu_i = muons->begin(); mu_i!= muons->end(); ++mu_i) {
        if (!goodSusyMuon(&(*mu_i))) continue; 
        ++nMuons;
      
        pt_muons_->Fill(  mu_i->pt() );
        eta_muons_->Fill( mu_i->eta() );
        phi_muons_->Fill( mu_i->phi() );

        reco::MuonIsolation muIso = mu_i->isolationR03();
        Iso_muons_->Fill( muIso.emEt+muIso.hadEt+muIso.sumPt );

        //Muon muon pairs 
        for(typename std::vector<Mu>::const_iterator mu_j = muons->begin(); mu_j!= muons->end(); ++mu_j) {
            if( mu_i>=mu_j ) continue;
            if (!goodSusyMuon(&(*mu_j))) continue; 
	        
            inv = (mu_i->p4()+mu_j->p4()).M();
            if( mu_i->charge()*mu_j->charge() > 0){
                ++nSSmumu;
                mass_SS_mumu_->Fill(inv);
            }
            if( mu_i->charge()*mu_j->charge() < 0){
                ++nOSmumu;
                mass_OS_mumu_->Fill(inv);
            }
        }
       
        //Electron muon pairs 
        for(typename std::vector<Ele>::const_iterator ele_j = elecs->begin(); ele_j!= elecs->end(); ++ele_j) {
            if (!goodSusyElectron(&(*ele_j))) continue; 
            inv = (mu_i->p4()+ele_j->p4()).M();
            dR = deltaR(*mu_i,*ele_j);
            dR_emu_->Fill(dR);
            if( mu_i->charge()*ele_j->charge() > 0){
                ++nSSemu;
                mass_SS_emu_->Fill(inv);
            }
            if( mu_i->charge()*ele_j->charge() < 0){
                ++nOSemu;
                mass_OS_emu_->Fill(inv);
            }
        }
    }
      
    N_muons_->Fill( nMuons );
  
    int nElectrons = 0;
    int nSSee = 0;
    int nOSee = 0;
    for(typename std::vector<Ele>::const_iterator ele_i = elecs->begin(); ele_i!= elecs->end(); ++ele_i) {
        if (!goodSusyElectron(&(*ele_i))) continue; 
        nElectrons++;

        pt_elecs_->Fill( ele_i->pt() );
        eta_elecs_->Fill(ele_i->eta());
        phi_elecs_->Fill(ele_i->phi());
      
        Iso_elecs_->Fill( ele_i->dr03TkSumPt()+ele_i->dr03EcalRecHitSumEt()+ele_i->dr03HcalTowerSumEt());
        
        //Electron electron pairs 
        for(typename std::vector<Ele>::const_iterator ele_j = elecs->begin(); ele_j!= elecs->end(); ++ele_j) {
            if( ele_i>=ele_j ) continue;
            if (!goodSusyElectron(&(*ele_j))) continue; 
            
            inv = (ele_i->p4()+ele_j->p4()).M();
            if( ele_i->charge()*ele_j->charge() > 0){
                ++nSSee;
                mass_SS_ee_->Fill(inv);
            }
            if( ele_i->charge()*ele_j->charge() < 0){
                ++nOSee;
                mass_OS_ee_->Fill(inv);
            }
        }
    }

    N_elecs_->Fill( nElectrons );

    if( MET > met_cut_ && sumPt > jet_sum_pt_cut_){
        if( nMuons>=1 ){ Muon_monitor_->Fill(sumPt,MET); }
        if( nElectrons>=1 ){ Electron_monitor_->Fill(sumPt,MET); }
        if( nOSee >=1 ){ OSee_monitor_->Fill(sumPt,MET); }
        if( nOSemu >=1 ){ OSemu_monitor_->Fill(sumPt,MET); }
        if( nOSmumu >=1 ){ OSmumu_monitor_->Fill(sumPt,MET); }
        if( nSSee >=1 ){ SSee_monitor_->Fill(sumPt,MET); }
        if( nSSemu >=1 ){ SSemu_monitor_->Fill(sumPt,MET); }
        if( nSSmumu >=1 ){ SSmumu_monitor_->Fill(sumPt,MET); }
    }
    if(nMuons>=3){
        TriMuon_monitor_->Fill(sumPt,MET);
    }

}
    
template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::endRun(const edm::Run& run) {

}

template< typename Mu, typename Ele, typename Jet, typename Met >
void SusyLeptonDQM<Mu,Ele,Jet,Met>::endJob() {

}

#endif
