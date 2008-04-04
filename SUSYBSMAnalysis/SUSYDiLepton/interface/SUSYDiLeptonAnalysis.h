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

using namespace edm;
using namespace reco;
using namespace std;

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
        TH1D*       hLeptonMult;
        TH1D*       hElectronMult;
        TH1D*       hMuonMult;
        TH1D*       hJetMult;

        TH1D*       hCleanElectronEoverP;
        TH1D*       hCleanMElectronEoverP;
        TH1D*       hCleanElectronEoverPinv;
        TH1D*       hCleanMElectronEoverPinv;
        TH1D*       hCleanElectronHoverE;
        TH1D*       hCleanMElectronHoverE;
        TH1D*       hCleanMuondpp;
        TH1D*       hCleanMMuondpp;
        TH1D*       hCleanMuonchi2;
        TH1D*       hCleanMMuonchi2;
        TH1D*       hCleanMuonNhits;
        TH1D*       hCleanMMuonNhits;
        TH1D*       hElectronIso;
        TH1D*       hMElectronIso;
        TH1D*       hElectrondeltaR;
        TH1D*       hMuonIso;
        TH1D*       hMMuonIso;
        TH1D*       hMuondeltaR;

        TH1D*       hMissingET;
        TH1D*       hMissingETmc;
        TH1D*       hEtSum;
        TH1D*       hHT;
        TH2D*       h2dMETptSumJets;
        TH2D*       h2dMETHT;
        TH2D*       h2dMETEtSum;

        TH1D*       hJetPt;
        TH1D*       hJetEta;
        TH1D*       hJetPhi;
        TH1D*       hPtJet1;
        TH1D*       hPtJet2;
        TH1D*       hPtJet3;
        TH1D*       hPtJet4;

        TH1D*       hInvMSFOS;
        TH1D*       hInvMOFOS;
        TH1D*       hInvMass;
        TH1D*       hInvMMuon;
        TH1D*       hInvMElectron;

        TH1D*       hMuonPt;
        TH1D*       hMuonEta;
        TH1D*       hMuonPhi;
        
        TH1D*       hElectronPt;
        TH1D*       hElectronEta;
        TH1D*       hElectronPhi;

        TH2D*       h2dLeptEtaPt;

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;

    std::string rootFileName;
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
    int cut_nJets;
    double cut_ptfirstJet;
    double cut_ptsecondJet;
    double cut_ptthirdJet;
    double cut_ptfourthJet;
    double cut_MET;

    double cut_MuonIso;
    double cut_ElectronIso;
    

  virtual bool CheckCuts(const edm::Handle< std::vector<pat::Muon> >, const edm::Handle< std::vector<pat::Electron> >, const edm::Handle< std::vector<pat::Jet> >, const edm::Handle< std::vector<pat::MET> >);

  virtual bool JetCut(const edm::Handle< std::vector<pat::Jet> >);
  virtual bool METCut(const edm::Handle< std::vector<pat::MET> >);
  virtual bool TwoLeptonCut(const edm::Handle< std::vector<pat::Muon> >, const edm::Handle< std::vector<pat::Electron> >);

  virtual void Analysis(const edm::Handle< std::vector<pat::Muon> >, const edm::Handle< std::vector<pat::Electron> >, const edm::Handle< std::vector<pat::Jet> >, const edm::Handle< std::vector<pat::MET> >, double weight);

  virtual double square(double);
  virtual double GetInvMass(const pat::Particle&, const pat::Particle&);

  virtual bool IsClean(const pat::Particle&);
  virtual bool IsIsolatedMuon(const pat::Muon&);
  virtual bool IsIsolatedElectron(const pat::Electron&);

};
#endif
