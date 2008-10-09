#ifndef SUSYDiLeptonAnalysis_h
#define SUSYDiLeptonAnalysis_h

/** \class SUSYDiLeptonAnalysis
 *
 *  
 *  This class is an EDAnalyzer for PAT
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 17.03.2008$
 *  $Revision: 0.1 $
 *  for CMSSW_1_6_10
 *  \author Niklas Mohr  --  niklas.mohr@cern.ch
 *
 */

// system include files
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <map>
#include <string>
#include <utility> 
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class SUSYDiLeptonAnalysis : public edm::EDAnalyzer {
    public:
    
    /// Constructor
    SUSYDiLeptonAnalysis(const edm::ParameterSet &iConfig);

    /// Destructor
    virtual ~SUSYDiLeptonAnalysis();

    /// Perform the real analysis
    void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

   private:
    //Histograms
        TH1F**       hLeptonMult;
        TH1F**       hElectronMult;
        TH1F**       hMuonMult;
        TH1F**       hJetMult;
        TH1F**       hbJetMult;

        TH1F**       hElectronIso;
        TH1F**       hMuonIso;

        TH1F**       hMissingET;
        TH1F**       hMissingETmc;
        TH1F**       hEtSum;
        TH1F**       hHT;
        TH2F**       h2dMETptSumJets;
        TH2F**       h2dMETHT;
        TH2F**       h2dMETEtSum;

        TH1F**       hJetPt;
        TH1F**       hJetEta;
        TH1F**       hJetPhi;
        TH1F**       hPtJet1;
        TH1F**       hPtJet2;
        TH1F**       hPtJet3;
        TH1F**       hPtJet4;

        TH1F**       hInvMSFOS;
        TH1F**       hInvMOFOS;
        TH1F**       hInvMass;
        TH1F**       hInvMMuon;
        TH1F**       hInvMMuonSS;
        TH1F**       hInvMElectron;
        TH1F**       hInvMElectronSS;
        
        TH1F**       hInvMbbllSFOS;
        TH1F**       hInvMbbllOFOS;

        TH2F**       h2dMETInvMassSFOS;
        TH2F**       h2dHTInvMassSFOS;
        TH2F**       h2dPtJetsInvMassSFOS;
        TH2F**       h2dPtJet4InvMassSFOS;

        TH2F**       h2dMETInvMassOFOS;
        TH2F**       h2dHTInvMassOFOS;
        TH2F**       h2dPtJetsInvMassOFOS;
        TH2F**       h2dPtJet4InvMassOFOS;

        TH1F**       hMuonPt;
        TH1F**       hMuonEta;
        TH1F**       hMuonPhi; 
	TH1F**       hGenMuonPt;
        TH1F**       hGenMuonEta;
        
        TH1F**       hElectronPt;
        TH1F**       hElectronEta;
        TH1F**       hElectronPhi;
        TH1F**       hGenElectronPt;
        TH1F**       hGenElectronEta;

        TH2F**       h2dMuonEtaPt;
        TH2F**       h2dIsoMuonEtaPt;
        TH2F**       h2dGenMuonEtaPt;
        TH2F**       h2dElectronEtaPt;
        TH2F**       h2dIsoElectronEtaPt;
        TH2F**       h2dGenElectronEtaPt;
        
	TH2F**       h2dTnPMuonEtaPt;
	TH2F**       h2dTnPIsoMuonEtaPt;
	
        TH1F**       hTrigger;

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;
    
    // Switch for mc analysis
    bool mcInfo;

    //std::string rootFileName;
    edm::InputTag mcSrc;
    edm::InputTag backmapSrc;
    edm::InputTag muonSrc;
    edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag jetObj;

    //CSA07 weighted dataset
    bool weighted;

    //Cuts for the analysis
    bool rej_Cuts;
    bool rej_JetCut;
    bool rej_METCut;
    bool rej_TwoLeptonCut;
    bool rej_bTagCut;
    int cut_nJets;
    int cut_nbJets;
    double cut_ptfirstJet;
    double cut_ptsecondJet;
    double cut_ptthirdJet;
    double cut_ptfourthJet;
    double cut_MET;
    double cut_bTagDiscriminator;

    double cut_MuonPt;
    double cut_MuonEta;
    double cut_MuonIso;
    
    double cut_ElectronPt;
    double cut_ElectronEta;
    double cut_ElectronIso;
    
    double cut_JetPt;
    double cut_JetEta;
    
    double cut_TnPlow;
    double cut_TnPhigh;
    
    std::string bJetAlgo;

  //Global counters
  int numTotEvents;
  int numTotEventsAfterCuts;

  int numTotEvents_tt;
  int numTotEventsAfterCuts_tt;
  int numTotEvents_W;
  int numTotEventsAfterCuts_W;
  int numTotEvents_Z;
  int numTotEventsAfterCuts_Z;

  int numTotElectrons;
  int numTotCleanElectrons;
  int numTotIsolatedElectrons;

  int numTotMuons;
  int numTotCleanMuons;
  int numTotIsolatedMuons;

  int numTotJets;
  int numTotCleanJets;

  inline void InitHisto(TFileDirectory *fs, const int process);

  virtual bool CheckCuts(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<pat::Jet> >&, const edm::Handle< std::vector<pat::MET> >&);

  virtual bool JetCut(const edm::Handle< std::vector<pat::Jet> >&);
  virtual bool METCut(const edm::Handle< std::vector<pat::MET> >&);
  virtual bool TwoLeptonCut(const edm::Handle< std::vector<pat::Muon> >& , const edm::Handle< std::vector<pat::Electron> >&);

  virtual void Analysis(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<pat::Jet> >&, const edm::Handle< std::vector<pat::MET> >&, double weight, const int process);

  virtual void MuonTnP(const edm::Handle< std::vector<pat::Muon> >&, double weight, const int process);
  
  virtual void TriggerMonitor(const edm::Handle< edm::TriggerResults >&, double weight, const int process);
  
  virtual void MCAnalysis(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<reco::GenParticle> >& , double weight, const int process); 

  virtual double square(double);
//  virtual double GetInvMass(const pat::Particle&, const pat::Particle&);

  virtual bool IsCleanMuon(const pat::Muon&);
  virtual bool IsCleanElectron(const pat::Electron&);
  virtual bool IsCleanJet(const pat::Jet& );
  virtual bool IsIsolatedMuon(const pat::Muon&);
  virtual bool IsIsolatedElectron(const pat::Electron&);
  virtual void PrintStatistics(void);

};
#endif
