#ifndef SUSYDiLeptonAnalysis_h
#define SUSYDiLeptonAnalysis_h

/** \class SUSYDiLeptonAnalysis
 *
 *  
 *  This class is an EDAnalyzer for PAT
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 2008/12/9 08:07:03 $
 *  $Revision: 1.4 $
 *  for CMSSW_2_2_3
 *  \author Niklas Mohr  --  niklas.mohr@cern.ch
 *
 */

// system include files
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <vector>
#include <map>
#include <utility> 
#include <memory>
#include <iostream>
#include <fstream>
#include <string>

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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

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

#include "NiklasMohr/SUSYDiLepton/interface/Combinations.h"


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
    TH1F**       hElectronTrackIso;
    TH1F**       hElectronCaloIso;
    TH1F**       hMuonIso;
    TH1F**       hMuonTrackIso;
    TH1F**       hMuonCaloIso;

    TH1F**       hMissingET;
    TH1F**       hMissingETmc;
    TH1F**       hEtSum;
    TH1F**       halphaT;
    TH1F**       hHT;
    TH1F**       hMHT;
    TH2F**       h2dMETEtSumJets;
    TH2F**       h2dMETHT;
    TH2F**       h2dMETMHT;

    TH1F**       hJetEt;
    TH1F**       hJetEta;
    TH1F**       hJetPhi;
    TH1F**       hEtJet1;
    TH1F**       hEtJet2;
    TH1F**       hEtJet3;
    TH1F**       hEtJet4;

    TH1F**       hInvMSFOS;
    TH1F**       hInvMOFOS;
    TH1F**       hInvMass;
    TH1F**       hInvMMuon;
    TH1F**       hInvMMuonSS;
    TH1F**       hInvMElectron;
    TH1F**       hInvMElectronSS;
    TH1F**       hInvMassMC;
        
    TH1F**       hInvMbbllSFOS;
    TH1F**       hInvMbbllOFOS;

    TH2F**       h2dMETInvMassSFOS;
    TH2F**       h2dHTInvMassSFOS;
    TH2F**       h2dEtJetsInvMassSFOS;
    TH2F**       h2dEtJet4InvMassSFOS;
    TH2F**       h2dIsoInvMassSFOS;

    TH2F**       h2dMETInvMassOFOS;
    TH2F**       h2dHTInvMassOFOS;
    TH2F**       h2dEtJetsInvMassOFOS;
    TH2F**       h2dEtJet4InvMassOFOS;
    TH2F**       h2dIsoInvMassOFOS;

    TH1F**       hMuonPt;
    TH1F**       hMuon1Pt;
    TH1F**       hMuon2Pt;
    TH1F**       hMuonEta;
    TH1F**       hMuonPhi; 
    TH1F**       hGenMuonPt;
    TH1F**       hGenMuonEta;
        
    TH1F**       hMuonChi2;
    TH1F**       hMuond0;
    TH1F**       hMuond0Sig;
    TH1F**       hMuonnHits;
    TH2F**       hMuonIsod0;
        
    TH1F**       hElectronPt;
    TH1F**       hElectron1Pt;
    TH1F**       hElectron2Pt;
    TH1F**       hElectronEta;
    TH1F**       hElectronPhi;
    TH1F**       hElectrond0;
    TH1F**       hGenElectronPt;
    TH1F**       hGenElectronEta;
        
    TH1F**       hElectronEoverP;
    TH1F**       hElectronfBrem;
    TH1F**       hElectronHoverE;
    TH1F**       hElectrondeltaPhiIn;
    TH1F**       hElectrondeltaEtaIn;
    TH2F**       hElectronIsod0;

    TH2F**       h2dMuonEtaPt;
    TH2F**       h2dMatchedMuonEtaPt;
    TH2F**       h2dGenMuonEtaPt;
    TH2F**       h2dElectronEtaPt;
    TH2F**       h2dMatchedElectronEtaPt;
    TH2F**       h2dGenElectronEtaPt;
        
    TH2F**       h2dTnPProbeSigEtaPt;
    TH2F**       h2dTnPPassSigEtaPt;
    TH2F**       h2dTnPProbeSSSigEtaPt;
    TH2F**       h2dTnPPassSSSigEtaPt;
    TH2F**       h2dTnPProbeSB1EtaPt;
    TH2F**       h2dTnPPassSB1EtaPt;
    TH2F**       h2dTnPProbeSSSB1EtaPt;
    TH2F**       h2dTnPPassSSSB1EtaPt;
    TH2F**       h2dTnPProbeSB2EtaPt;
    TH2F**       h2dTnPPassSB2EtaPt;
    TH2F**       h2dTnPProbeSSSB2EtaPt;
    TH2F**       h2dTnPPassSSSB2EtaPt;
    TH1F**       hTnPProbeInvMass;
    TH1F**       hTnPPassInvMass;
    TH1F**       hTnPSSInvMass;
	
    TH1F**       hTrigger;
    TH1F**       hWeight;
    
    //Trees
    TTree*       treeOFOS;
    TTree*       treeSFOS;
    TTree*       treeMuon;
    TTree*       treeElec;
    float invMOFOS;
    float invMSFOS;
    float invMMuon;
    float invMElec;
    float invweight;

    // The file which will store the histos
    TFile *theFile;

    // Switch for debug output
    bool debug;
    
    // Switch for mc analysis
    bool mcInfo;

    //std::string rootFileName;
    edm::InputTag mcSrc;
    edm::InputTag beamSpotSrc;
    edm::InputTag muonSrc;
    edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag jetObj;
    edm::InputTag trackSrc;

    //CSA07 weighted dataset
    bool weighted;
    
    //Set weight via cfg
    double externalWeight;
    bool Signal_Analysis;

    //Cuts for the analysis
    bool rej_Cuts;
    bool rej_JetCut;
    bool rej_METCut;
    std::string rej_LeptonCut;
    bool rej_bTagCut;
    int cut_nPos_Elec;
    int cut_nNeg_Elec;
    int cut_nPos_Muon;
    int cut_nNeg_Muon;
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
    double cut_MuonChi2;
    double cut_Muond0;
    double cut_MuonnHits;
    double cut_MuonIso;
    
    std::string ElectronID;
    double cut_ElectronPt;
    double cut_ElectronEta;
    double cut_Electrond0;
    double cut_ElectronIso;
    
    double cut_JetPt;
    double cut_JetEta;
   
    std::string methodTnP; 
    double cut_TnPChi;
    double cut_TnPnHits;
    double cut_TnPDr;
    double cut_TnPTagPt;
    double cut_TnPTagEta;
    double cut_TnPlowSB1;
    double cut_TnPhighSB1;
    double cut_TnPlowSB2;
    double cut_TnPhighSB2;
    double cut_TnPlow;
    double cut_TnPhigh;

    
    double electron_Scale;
    double jet_Scale;
    
    std::string bJetAlgo;

    math::XYZPoint bs;

    static const int nEtaBins=5;
    static const float boundEta[nEtaBins+1]; 
    static const int nPtBins=4;
    static const float boundPt[nPtBins+1]; 
    float Muon_Eff[nEtaBins][nPtBins];
    float Electron_Eff[nEtaBins][nPtBins];

    std::vector<unsigned int> mlaMinDpt; //!< list of good jet indices forming pseudo-jet a (MinDEt)
    std::vector<unsigned int> mlbMinDpt; //!< list of good jet indices forming pseudo-jet b (MinDEt)
    
    int general;
    int clean;
    int all;
    int effcor;
    int unmatched;
    int promt;
    int decay;
  
    //Global counters
    int numTotEvents;
    int numTotEventsAfterCuts;

    int numTotElectrons;
    int numTotCleanElectrons;
    int numTotIsolatedElectrons;

    int numTotMuons;
    int numTotCleanMuons;
    int numTotIsolatedMuons;

    int numTotJets;
    int numTotCleanJets;

    inline void InitHisto(TFileDirectory *fs, const int process);
    inline void ReadEfficiency();

    virtual bool CheckCuts(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<pat::Jet> >&, const edm::Handle< std::vector<pat::MET> >&);

    virtual bool JetCut(const edm::Handle< std::vector<pat::Jet> >&);
    virtual bool METCut(const edm::Handle< std::vector<pat::MET> >&, const edm::Handle< std::vector<pat::Jet> >&);
    virtual bool LeptonCut(const edm::Handle< std::vector<pat::Muon> >& , const edm::Handle< std::vector<pat::Electron> >&);

    virtual void Analysis(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<pat::Jet> >&, const edm::Handle< std::vector<pat::MET> >&, double weight, const int process);

    virtual void MuonTnP(const edm::Handle< std::vector<pat::Muon> >&,const edm::Handle< std::vector<reco::Track> >& , double weight, const int process);
    virtual void ElectronTnP(const edm::Handle< std::vector<pat::Electron> >&,const edm::Handle< std::vector<reco::Track> >& , double weight, const int process);
  
    virtual void TriggerMonitor(const edm::Handle< edm::TriggerResults >&, double weight, const int process);
  
    virtual void ElectronMonitor(const pat::Electron*, const int, double, const int);
    virtual void MuonMonitor(const pat::Muon*, const int, double, const int);

    virtual void MuonInvMonitor(const double, const double, const double, const double, const double, const double, double, const int);
    virtual void OFOSInvMonitor(const double, const double, const double, const double, const double, const double, double, const int);
    virtual void ElectronInvMonitor(const double, const double, const double, const double, const double, const double, double, const int);
  
    virtual bool MCAnalysis(const edm::Handle< std::vector<pat::Muon> >&, const edm::Handle< std::vector<pat::Electron> >&, const edm::Handle< std::vector<reco::GenParticle> >& , double weight, const int process); 

    virtual bool FindMinCombo(std::vector<const reco::Candidate* >, std::vector<unsigned int>&, std::vector<unsigned int>&);
    virtual double CalcalphaT(std::vector<const reco::Candidate* >, std::vector<unsigned int>, std::vector<unsigned int> );

    virtual double square(double);
    //virtual double GetInvMass(double, double, double, double, double, double, double, double);

    virtual double getMuonWeight(const pat::Muon*);
    virtual double getElectronWeight(const pat::Electron*);
    virtual bool IsCleanMuon(const pat::Muon&);
    virtual bool IsCleanElectron(const pat::Electron&);
    virtual bool IsCleanJet(const pat::Jet& );
    virtual bool IsIsolatedMuon(const pat::Muon&);
    virtual bool IsIsolatedElectron(const pat::Electron&);
    virtual void PrintStatistics(void);

    //template < class T > int SUSYDiLeptonAnalysis::GetLeptKind(const T *);

};

const float SUSYDiLeptonAnalysis::boundEta[] = {-2.,-1.2,-0.4,0.4,1.2,2};
const float SUSYDiLeptonAnalysis::boundPt[] = {0.,40.,45.,50.,1000.};

#endif
