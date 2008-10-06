/** \class SUSYDiLeptonAnalysis
 *
 *  
 *  This class is an EDAnalyzer for PAT 
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 17.05.2008$
 *  $Revision: 0.3 $ for CMSSW 1_6_11/2_0_6
 *
 *  \author: Niklas Mohr -- niklas.mohr@cern.ch
 *  
 */

#include "TestSystem/SUSYDiLepton/interface/SUSYDiLeptonAnalysis.h"

using namespace edm;
using namespace reco;
using namespace std;

//Constructor
SUSYDiLeptonAnalysis::SUSYDiLeptonAnalysis(const edm::ParameterSet &iConfig)
{
  //now do what ever initialization is needed
  //Debug flag
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  
  //Monte carlo information
  mcInfo            = iConfig.getUntrackedParameter<bool>   ("mcInfo");

  //Input collections
  mcSrc             = iConfig.getParameter<edm::InputTag> ("mcSource");
  backmapSrc        = iConfig.getParameter<edm::InputTag> ("backmapSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
 
  //Cuts
  weighted         = iConfig.getUntrackedParameter<bool>   ("CSA_weighted");
  rej_Cuts         = iConfig.getUntrackedParameter<bool>   ("rej_Cuts");
  rej_JetCut       = iConfig.getUntrackedParameter<bool>   ("rej_JetCut");
  rej_METCut       = iConfig.getUntrackedParameter<bool>   ("rej_METCut");
  rej_TwoLeptonCut = iConfig.getUntrackedParameter<bool>   ("rej_TwoLeptonCut");
  rej_bTagCut      = iConfig.getUntrackedParameter<bool>   ("rej_bTagCut");
  cut_nJets        = iConfig.getUntrackedParameter<int> ("user_nJets");
  cut_nbJets       = iConfig.getUntrackedParameter<int> ("user_nbJets");
  cut_ptfirstJet   = iConfig.getUntrackedParameter<double> ("user_pt1JetMin");
  cut_ptsecondJet  = iConfig.getUntrackedParameter<double> ("user_pt2JetMin");
  cut_ptthirdJet   = iConfig.getUntrackedParameter<double> ("user_pt3JetMin");
  cut_ptfourthJet  = iConfig.getUntrackedParameter<double> ("user_pt4JetMin");
  cut_MET	   = iConfig.getUntrackedParameter<double> ("user_METMin");
  cut_bTagDiscriminator = iConfig.getUntrackedParameter<double> ("user_bTagDiscriminator");

  //bJetAlgo  = iConfig.getUntrackedParameter<string> ("user_") 
  bJetAlgo = "jetProbabilityJetTags";
  
  cut_MuonPt	   = iConfig.getUntrackedParameter<double> ("acc_MuonPt");
  cut_MuonEta	   = iConfig.getUntrackedParameter<double> ("acc_MuonEta");
  cut_MuonIso	   = iConfig.getUntrackedParameter<double> ("iso_MuonIso");
  
  cut_ElectronPt  = iConfig.getUntrackedParameter<double> ("acc_ElectronPt");
  cut_ElectronEta  = iConfig.getUntrackedParameter<double> ("acc_ElectronEta");
  cut_ElectronIso  = iConfig.getUntrackedParameter<double> ("iso_ElectronIso");
  
  cut_JetPt  = iConfig.getUntrackedParameter<double> ("acc_JetPt");
  cut_JetEta  = iConfig.getUntrackedParameter<double> ("acc_JetEta");
  
  //cut_TnPlow  = iConfig.getUntrackedParameter<double> ("tnp_low_inv");
  //cut_TnPhigh  = iConfig.getUntrackedParameter<double> ("tnp_high_inv");
  cut_TnPlow = 50;
  cut_TnPhigh = 120;
 
  //initialize global counters
  numTotEvents = 0;
  numTotEventsAfterCuts = 0;

  numTotEvents_tt = 0;
  numTotEventsAfterCuts_tt = 0;
  
  numTotEvents_W = 0;
  numTotEventsAfterCuts_W = 0;
  
  numTotEvents_Z = 0;
  numTotEventsAfterCuts_Z = 0;

  numTotElectrons = 0;
  numTotCleanElectrons = 0;
  numTotIsolatedElectrons = 0;
  
  numTotMuons = 0;
  numTotCleanMuons = 0;
  numTotIsolatedMuons = 0;
  
  numTotJets = 0;
  numTotCleanJets = 0;

  const int nHistos=4;
 
  // Create the root file
  edm::Service<TFileService> theFile;

  // book histograms for multiplicities of leptons and jets
  hLeptonMult = new TH1F * [nHistos];
  hElectronMult = new TH1F * [nHistos];
  hMuonMult = new TH1F * [nHistos];
  hJetMult = new TH1F * [nHistos];

  //histograms for lepton isolation cuts
  hElectronIso = new TH1F * [nHistos];
  hMuonIso = new TH1F * [nHistos];
        
  //histograms for the invariant mass of the leptons
  hInvMSFOS = new TH1F * [nHistos];
  hInvMOFOS = new TH1F * [nHistos];
  hInvMass = new TH1F * [nHistos];
  hInvMElectron = new TH1F * [nHistos];
  hInvMMuon = new TH1F * [nHistos];
  
  hInvMbbllSFOS = new TH1F * [nHistos];
  hInvMbbllOFOS = new TH1F * [nHistos];

  h2dMETInvMassSFOS = new TH2F * [nHistos];
  h2dHTInvMassSFOS = new TH2F * [nHistos];
  h2dPtJetsInvMassSFOS = new TH2F * [nHistos];
  h2dPtJet4InvMassSFOS = new TH2F * [nHistos];

  h2dMETInvMassOFOS = new TH2F * [nHistos];
  h2dHTInvMassOFOS = new TH2F * [nHistos];
  h2dPtJetsInvMassOFOS = new TH2F * [nHistos];
  h2dPtJet4InvMassOFOS = new TH2F * [nHistos];

  //muon histograms
  hMuonPt = new TH1F * [nHistos];
  hMuonEta = new TH1F * [nHistos];
  hMuonPhi = new TH1F * [nHistos];
  hGenMuonPt = new TH1F * [nHistos];
  hGenMuonEta = new TH1F * [nHistos];

  //electron histograms
  hElectronPt = new TH1F * [nHistos];
  hElectronEta = new TH1F * [nHistos];
  hElectronPhi = new TH1F * [nHistos];
  hGenElectronPt = new TH1F * [nHistos];
  hGenElectronEta = new TH1F * [nHistos];

  h2dMuonEtaPt = new TH2F * [nHistos];
  h2dIsoMuonEtaPt = new TH2F * [nHistos];
  h2dGenMuonEtaPt = new TH2F * [nHistos];
  h2dElectronEtaPt = new TH2F * [nHistos];
  h2dIsoElectronEtaPt = new TH2F * [nHistos];
  h2dGenElectronEtaPt = new TH2F * [nHistos];
  
  //TnP
  h2dTnPMuonEtaPt = new TH2F * [nHistos];
  h2dTnPIsoMuonEtaPt = new TH2F * [nHistos];
 
  //histograms for Missing ET
  hMissingET = new TH1F * [nHistos];
  hMissingETmc =  new TH1F * [nHistos];
  hEtSum = new TH1F * [nHistos];
  hHT = new TH1F * [nHistos];
  h2dMETptSumJets = new TH2F * [nHistos];
  h2dMETHT = new TH2F * [nHistos];
  h2dMETEtSum = new TH2F * [nHistos];

  //histograms for jets
  hJetPt = new TH1F * [nHistos];
  hJetEta = new TH1F * [nHistos];
  hJetPhi = new TH1F * [nHistos];
  hPtJet1 = new TH1F * [nHistos];
  hPtJet2 = new TH1F * [nHistos];
  hPtJet3 = new TH1F * [nHistos];
  hPtJet4 = new TH1F * [nHistos];
  
  hTrigger = new TH1F * [nHistos];

  TFileDirectory General = theFile->mkdir( "General" );
  TFileDirectory Wjets = theFile->mkdir( "Wjets" );
  TFileDirectory Zjets = theFile->mkdir( "Zjets" );
  TFileDirectory ttjets = theFile->mkdir( "ttjets" );

  InitHisto(&General,0);
  InitHisto(&Wjets,1);
  InitHisto(&Zjets,2);
  InitHisto(&ttjets,3);

}

void inline SUSYDiLeptonAnalysis::InitHisto(TFileDirectory *theFile, const int process)
{
  //TFileDirectory test = theFile->mkdir( "test" );
  // book histograms for multiplicities of leptons and jets
  hLeptonMult[process] = theFile->make<TH1F>( "LeptonMultiplicity", "Multiplicity of electrons + muons", 15, 0.0, 15.0);
  hElectronMult[process] = theFile->make<TH1F>( "ElectronMultiplicity", "Multiplicity of electrons", 10, 0.0, 10.0);
  hMuonMult[process] = theFile->make<TH1F>( "MuonMultiplicity", "Multiplicity of muons", 10, 0.0, 10.0);
  hJetMult[process] = theFile->make<TH1F>( "JetMultiplicity", "Multiplicity of jets", 30, 0.0, 30.0);

  //histograms for lepton isolation cuts
  hElectronIso[process] = theFile->make<TH1F>( "ElectronIso", "Isolation of electrons", 1000, 0.0, 10.0); 
  hMuonIso[process] = theFile->make<TH1F>( "MuonIso", "Isolation of muons", 1000, 0.0, 10.0);
        
  //histograms for the invariant mass of the leptons
  hInvMSFOS[process] = theFile->make<TH1F>( "Invariant mass of SFOS lepton pairs", "Invariant mass of SFOS lepton pairs", 300, 0, 300);
  hInvMOFOS[process] = theFile->make<TH1F>( "Invariant mass of OFOS lepton pairs", "Invariant mass of OFOS lepton pairs", 300, 0, 300);
  hInvMass[process] = theFile->make<TH1F>( "Invariant mass of lepton pairs", "Invariant mass of lepton pairs", 300, 0, 300);
  hInvMElectron[process] = theFile->make<TH1F>( "Invariant mass of electron pairs", "Invariant mass of electron pairs", 300, 0, 300);
  hInvMMuon[process] = theFile->make<TH1F>( "Invariant mass of muon pairs", "Invariant mass of muon pairs", 300, 0, 300);
  
  //histograms for bbll inv mass
  hInvMbbllSFOS[process] = theFile->make<TH1F>( "Invariant mass of SFOS bbll", "Invariant mass of SFOS bbll", 1500, 0, 1500);
  hInvMbbllOFOS[process] = theFile->make<TH1F>( "Invariant mass of OFOS bbll", "Invariant mass of OFOS bbll", 1500, 0, 1500);

  h2dMETInvMassSFOS[process] = theFile->make<TH2F>( "MET - Invariant mass SFOS", "MET - Invariant mass SFOS", 300, 0., 300., 2000, 0., 1000.0);
  h2dHTInvMassSFOS[process] = theFile->make<TH2F>( "HT - Invariant mass SFOS", "HT - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
  h2dPtJetsInvMassSFOS[process] = theFile->make<TH2F>( "SumJets - Invariant mass SFOS", "SumJets - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
  h2dPtJet4InvMassSFOS[process] = theFile->make<TH2F>( "PtJet4 - Invariant mass SFOS", "PtJet4 - Invariant mass SFOS", 300, 0., 300., 2500, 0., 1500.0);

  h2dMETInvMassOFOS[process] = theFile->make<TH2F>( "MET - Invariant mass OFOS", "MET - Invariant mass OFOS", 300, 0., 300., 2000, 0., 1000.0);
  h2dHTInvMassOFOS[process] = theFile->make<TH2F>( "HT - Invariant mass OFOS", "HT - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
  h2dPtJetsInvMassOFOS[process] = theFile->make<TH2F>( "SumJets - Invariant mass OFOS", "SumJets - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
  h2dPtJet4InvMassOFOS[process] = theFile->make<TH2F>( "PtJet4 - Invariant mass OFOS", "PtJet4 - Invariant mass OFOS", 300, 0., 300., 1500, 0., 1500.0);

  //muon histograms
  hMuonPt[process] = theFile->make<TH1F>( "muon pt", "muon pt", 1000, 0.0, 1000.0);
  hMuonEta[process] = theFile->make<TH1F>( "muon eta", "muon eta", 250, -2.5, 2.5);
  hMuonPhi[process] = theFile->make<TH1F>( "muon phi", "muon phi", 350, -3.5, 3.5);
  hGenMuonPt[process] = theFile->make<TH1F>( "Generator muon pt", "Generator muon pt", 1000, 0.0, 1000.0);
  hGenMuonEta[process] = theFile->make<TH1F>( "Generator muon eta", "Generator muon eta", 250, -2.5, 2.5);

  //electron histograms
  hElectronPt[process] = theFile->make<TH1F>( "electron pt", "electron pt", 1000, 0.0, 1000.0);
  hElectronEta[process] = theFile->make<TH1F>( "electron eta", "electron eta", 250, -2.5, 2.5);
  hElectronPhi[process] = theFile->make<TH1F>( "electron phi", "electron phi", 350, -3.5, 3.5);
  hGenElectronPt[process] = theFile->make<TH1F>( "Generator electron pt", "Generator electron pt", 1000, 0.0, 1000.0);
  hGenElectronEta[process] = theFile->make<TH1F>( "Generator electron eta", "Generator electron eta", 250, -2.5, 2.5);

  h2dMuonEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of muons", "Eta - Pt of muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dIsoMuonEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of isolated muons", "Eta - Pt of isolated muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dGenMuonEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of MC muons", "Eta - Pt of generator muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dElectronEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of electrons", "Eta - Pt of electrons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dIsoElectronEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of isolated electrons", "Eta - Pt of isolated electrons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dGenElectronEtaPt[process] = theFile->make<TH2F>( "Eta-Pt of MC electrons", "Eta - Pt of generator electrons", 1000, 0.0, 1000.0, 250, -2.0, 2.5); 
  
  //TnP
  h2dTnPMuonEtaPt[process] = theFile->make<TH2F>( "Eta-Pt TnP muons", "Eta - Pt of TnP muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  h2dTnPIsoMuonEtaPt[process] = theFile->make<TH2F>( "Eta-Pt isolated TnP  muons", "Eta - Pt of isolated TnP  muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
 
  
  //histograms for Missing ET
  hMissingET[process] = theFile->make<TH1F>( "Missing transverse energy", "Missing transverse energy", 2000, 0.0, 1000.0);
  hMissingETmc[process] =  theFile->make<TH1F>( "Missing transverse energy MC", "Missing transverse energy MC", 2000, 0.0, 1000.0);
  hEtSum[process] = theFile->make<TH1F>( "ETsum", "Transverse energy sum ET", 2000, 0.0, 2000.0);
  hHT[process] = theFile->make<TH1F>( "HT", "Transverse energy sum HT", 2000, 0.0, 2000.0);
  h2dMETptSumJets[process] = theFile->make<TH2F>( "MET - sumJets", "MET - sumJets", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
  h2dMETHT[process] = theFile->make<TH2F>( "MET - HT", "MET - HT", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
  h2dMETEtSum[process] = theFile->make<TH2F>( "MET - ET", "MET - ET", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);

  //histograms for jets
  hJetPt[process] = theFile->make<TH1F>( "jet Et", "jet Et", 1500, 0.0, 1500.0);
  hJetEta[process] = theFile->make<TH1F>( "jet eta", "jet eta", 300, -3., 3.);
  hJetPhi[process] = theFile->make<TH1F>( "jet phi", "jet phi", 350, -3.5, 3.5);
  hPtJet1[process] = theFile->make<TH1F>( "Et first jet", "Et spectrum of the 1st jet", 1500, 0.0, 1500.0);
  hPtJet2[process] = theFile->make<TH1F>( "Et second jet", "Et spectrum of the 2nd jet", 1500, 0.0, 1500.0);
  hPtJet3[process] = theFile->make<TH1F>( "Et third jet", "Et spectrum of the 3rd jet", 1500, 0.0, 1500.0);
  hPtJet4[process] = theFile->make<TH1F>( "Et fourth jet", "Et spectrum of the 4th jet", 1500, 0.0, 1500.0);
  
  hTrigger[process] = theFile->make<TH1F>( "Trigger paths", "Trigger paths", 160, 0, 160);

}

//Destructor
SUSYDiLeptonAnalysis::~SUSYDiLeptonAnalysis()
{
   PrintStatistics();
   if (debug) cout << "************* Finished analysis" << endl;
} 

//SUSYDiLeptonAnalysis::IsCleanMuon
//-----------------------------------------------
//Check if Electron is clean and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanMuon(const pat::Muon& muon)
{
 //if (muon.combinedMuon()->normalizedChi2()<5&&fabs(muon.track()->d0())<0.25&&muon.track()->numberOfValidHits() >= 7&&muon.pt()>cut_MuonPt&&abs(muon.eta())<cut_MuonEta){return true;} else return false;
 if (muon.track()->normalizedChi2()<5&&fabs(muon.track()->d0())<0.25&&muon.track()->numberOfValidHits() >= 7&&muon.pt()>cut_MuonPt&&abs(muon.eta())<cut_MuonEta){return true;} else return false;
}

//SUSYDiLeptonAnalysis::IsCleanElectron
//-----------------------------------------------
//Check if Electron is clean and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanElectron(const pat::Electron& electron)
{
 if (electron.leptonID("tight")==1.0&&electron.pt()>cut_ElectronPt&&abs(electron.eta())<cut_ElectronEta){return true;} else return false;
}

//SUSYDiLeptonAnalysis::IsCleanJet
//-----------------------------------------------
//Check if Jet is clean and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanJet(const pat::Jet& jet)
{
 if (pat::Flags::test(jet, pat::Flags::Overlap::Electrons)&&jet.pt()>cut_JetPt&&abs(jet.eta())<cut_JetEta){return true;} else return false;
}

//SUSYDiLeptonAnalysis::IsIsolatedMuon
//-----------------------------------------------
//Check if Muon is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedMuon(const pat::Muon& muon)
{
 if (debug) cout << "Muon tracker isolation: "<< muon.trackIso() << endl;
 if (debug) cout << "Muon pt: "<< muon.pt() << endl;
 if (muon.trackIso()<cut_MuonIso){return true;} else return false;
}

//SUSYDiLeptonAnalysis::IsIsolatedElectron
//-----------------------------------------------
//Check if Electron is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedElectron(const pat::Electron& electron)
{
 if (debug) cout << "Electron tracker isolation: "<< electron.trackIso() << endl;
 if (debug) cout << "Electron pt: "<< electron.pt() << endl;
 if (electron.trackIso()<cut_ElectronIso){return true;} else return false;
}

//SUSYDiLeptonAnalysis::JetCut
//-----------------------------------------------
//check if four jets above cuts are existent
bool SUSYDiLeptonAnalysis::JetCut(const edm::Handle< std::vector<pat::Jet> >& jets){
  if (debug) cout << "Checking Jet cut: " << endl;
  int n_Jet=0;
  int n_bJet=0;
  double pt4Jets=0;
  bool firstJet=false;
  bool secondJet=false;
  bool thirdJet=false;
  bool fourthJet=false;
  bool bJets=false;
  for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
	//Test only clean jets in acceptance
	if(IsCleanJet(*jet_i)){
		++n_Jet;
		if(n_Jet<=4){pt4Jets+=jet_i->pt();}
		if(n_Jet==1 && jet_i->pt()>cut_ptfirstJet){firstJet=true;}
		if(n_Jet==2 && jet_i->pt()>cut_ptsecondJet){secondJet=true;}
		if(n_Jet==3 && jet_i->pt()>cut_ptthirdJet){thirdJet=true;}
		if(n_Jet==4 && jet_i->pt()>cut_ptfourthJet){fourthJet=true;}
		if(jet_i->bDiscriminator(bJetAlgo)>cut_bTagDiscriminator){
			++n_bJet;
		}
	}
   }
   if (debug&&firstJet) cout <<" 1st Jet above cut " << endl;
   if (debug&&secondJet) cout <<" 2nd Jet above cut " << endl;
   if (debug&&thirdJet) cout <<" 3rd Jet above cut " << endl;
   if (debug&&fourthJet) cout <<" 4th Jet above cut " << endl;
   //Check if 2 bjets are in event
   if(rej_bTagCut&&n_bJet==cut_nbJets){bJets = true;}
   if(!rej_bTagCut){bJets = true;}
   //Check case of 1 Jet
   if (cut_nJets==0){
	if((pt4Jets>cut_ptfirstJet+cut_ptsecondJet+cut_ptthirdJet+cut_ptfourthJet)&&bJets) {return true;} else return false;
   }
   if (cut_nJets==1){
	if(firstJet&&bJets) {return true;} else return false;
   }
   //Check case of 2 Jets
   if (cut_nJets==2){
	if(firstJet&&secondJet&&bJets) {return true;} else return false;
   }
   //Check case of 3 Jets
   if (cut_nJets==3){
	if(firstJet&&secondJet&&thirdJet&&bJets) {return true;} else return false;
   }
   //Check case of 4 Jets
   if (cut_nJets==4){
	if(firstJet&&secondJet&&thirdJet&&fourthJet&&bJets) {return true;} else return false;
   }
   if (cut_nJets!=0||cut_nJets!=1||cut_nJets!=2||cut_nJets!=3||cut_nJets!=4) {return false;}

}

//SUSYDiLeptonAnalysis::METCut
//-----------------------------------------------
//check if MET above cut is existent
bool SUSYDiLeptonAnalysis::METCut(const edm::Handle< std::vector<pat::MET> >& met){
  if (debug) cout << "Checking MET cut: " << endl;
  bool MET=false;
  for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); ++met_i){      
	if (met_i->pt()>cut_MET) {MET=true;} 
  }
  if (MET) {return true;} else return false;
}

//SUSYDiLeptonAnalysis::TwoLeptonCut
//-----------------------------------------------
//check if Two leptons are the event
bool SUSYDiLeptonAnalysis::TwoLeptonCut(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons){
  int n_Muons = 0, n_Electrons = 0;
  if (debug) cout << "Checking Two Lepton cut: " << endl;
  for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
	if(IsCleanMuon(*mu_i)&&IsIsolatedMuon(*mu_i)){++n_Muons;}
  }
  for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
	if(IsCleanElectron(*ele_i)&&IsIsolatedElectron(*ele_i)){++n_Electrons;}
  }
  if(n_Muons+n_Electrons>=2) {return true;} else return false;
}

//SUSYDiLeptonAnalysis::CheckCuts
//-----------------------------------------------
//Check all cuts
bool SUSYDiLeptonAnalysis::CheckCuts(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle< std::vector<pat::MET> >& met){
   if (rej_Cuts){
	bool passedJetCut=false;
	bool passedMETCut=false;
	bool passedTwoLeptonCut=false;
   	if (rej_JetCut){
		if (JetCut(jets)){
			if (debug) cout << "Jet cut passed" << endl;
			passedJetCut = true;
   		} 
        } else passedJetCut = true;
   	if (rej_METCut){
		if (METCut(met)){
			if (debug) cout << "MET cut passed" << endl;
        		passedMETCut = true;
		} 
        } else passedMETCut = true;
   	if (rej_TwoLeptonCut){
   		if (TwoLeptonCut(muons,electrons)){
			if (debug) cout << "Two Lepton cut passed" << endl;
        		passedTwoLeptonCut = true;
   		}
	} else passedTwoLeptonCut = true; 
   	if (passedJetCut&&passedMETCut&&passedTwoLeptonCut){return true;} else return false;
   }
   else return true;  
}

//SUSYDiLeptonAnalysis::analysis
//-----------------------------------------------
//Filling of all histograms and calculation of kinematics
void SUSYDiLeptonAnalysis::Analysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle< std::vector<pat::MET> >& met, double weight, const int process){
  
   double MET=0;
   for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); ++met_i){
       	if (debug) cout <<"MET pt = "<< met_i->pt() << endl;
        MET = met_i->pt();
   }

   hMissingET[process]->Fill(MET,weight);
  
  int n_Jet=0;
  int n_bJet=0;
  std::vector< pat::Jet > bJets;
  double pt4Jets=0;
  for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
       	if (debug) cout <<"Jet pt = "<< jet_i->pt() << endl;
	++numTotJets;
	if(IsCleanJet(*jet_i)){
		++numTotCleanJets;
		++n_Jet;
		//Plots of leading Jets
		if(n_Jet==1){hPtJet1[process]->Fill(jet_i->pt(),weight);
			     pt4Jets+=jet_i->pt();}
		if(n_Jet==2){hPtJet2[process]->Fill(jet_i->pt(),weight);
			     pt4Jets+=jet_i->pt();}
		if(n_Jet==3){hPtJet3[process]->Fill(jet_i->pt(),weight);
			     pt4Jets+=jet_i->pt();}
		if(n_Jet==4){hPtJet4[process]->Fill(jet_i->pt(),weight);
			     pt4Jets+=jet_i->pt();}
		//Jet base plots
		if(jet_i->bDiscriminator(bJetAlgo)>cut_bTagDiscriminator){
			++n_bJet;
			bJets.push_back( *jet_i );
		}
		hJetPt[process]->Fill(jet_i->pt(),weight);
		hJetEta[process]->Fill(jet_i->eta(),weight);
		hJetPhi[process]->Fill(jet_i->phi(),weight);
	}
   }
   hJetMult[process]->Fill(n_Jet,weight);
   if (debug) cout <<" Number of Jets = " << n_Jet << endl;
   
   //Muon histograms
   int n_Muons = 0;
   for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
       	if (debug) cout <<"mu eta = "<< mu_i->eta() << endl;
	++numTotMuons;
	if(IsCleanMuon(*mu_i)){
	++numTotCleanMuons;
        h2dMuonEtaPt[process]->Fill(mu_i->pt(),mu_i->eta(),weight);
	}
	if(IsCleanMuon(*mu_i)&&IsIsolatedMuon(*mu_i)){
	++numTotIsolatedMuons;
	++n_Muons;
	//Muon base plots
	hMuonPt[process]->Fill(mu_i->pt(),weight);
	hMuonEta[process]->Fill(mu_i->eta(),weight);
	hMuonPhi[process]->Fill(mu_i->phi(),weight);
        h2dIsoMuonEtaPt[process]->Fill(mu_i->pt(),mu_i->eta(),weight);	
	//Muon isolation
        hMuonIso[process]->Fill(mu_i->trackIso(),weight);
	//Invariant mass plots
	//Muon pairs
   	for (std::vector<pat::Muon>::const_iterator mu_j = muons->begin(); mu_j != muons->end(); ++mu_j){
		if(IsCleanMuon(*mu_j)&&IsIsolatedMuon(*mu_j)){
       		if( mu_i->charge()==+1 && mu_j->charge()==-1){
			if (debug) cout <<"Invariant Mass mu+ mu-: " << (mu_i->p4()+mu_j->p4()).M() <<endl;
			hInvMass[process]->Fill((mu_i->p4()+mu_j->p4()).M(),weight);
			if(n_bJet>=2){
			for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                		for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
					if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+mu_j->p4()).M(),weight);}
				}
			}
			}
			hInvMMuon[process]->Fill((mu_i->p4()+mu_j->p4()).M(),weight);
			hInvMSFOS[process]->Fill((mu_i->p4()+mu_j->p4()).M(),weight);
			h2dMETInvMassSFOS[process]->Fill((mu_i->p4()+mu_j->p4()).M(),MET,weight);
			h2dPtJetsInvMassSFOS[process]->Fill((mu_i->p4()+mu_j->p4()).M(),pt4Jets,weight);
		}
		}
	}
   	//Wrong pairings for different flavour subtraction
	for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
		if(IsCleanElectron(*ele_j)&&IsIsolatedElectron(*ele_j)){
       		if(( mu_i->charge()==+1 && ele_j->charge()==-1)|(mu_i->charge()==-1 && ele_j->charge()==+1)){
			if (debug) cout <<"Invariant Mass mu+ e-, mu- e+: " << (mu_i->p4()+ele_j->p4()).M() <<endl;
			hInvMass[process]->Fill((mu_i->p4()+ele_j->p4()).M(),-weight);
			if(n_bJet>=2){
			for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                		for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
					if(bjet_i<bjet_j){hInvMbbllOFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+ele_j->p4()).M(),weight);}
				}
			}
			}
			hInvMOFOS[process]->Fill((mu_i->p4()+ele_j->p4()).M(),weight);
			h2dMETInvMassOFOS[process]->Fill((mu_i->p4()+ele_j->p4()).M(),pt4Jets,weight);
			h2dPtJetsInvMassOFOS[process]->Fill((mu_i->p4()+ele_j->p4()).M(),pt4Jets,weight);
		}
		}
	}
	}
   }
  //Muon multiplicity 
  hMuonMult[process]->Fill(n_Muons,weight);

  //Loop over electrons 
  int n_Electrons = 0;
  for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
       	if (debug) cout <<"ele eta = "<< ele_i->eta() << endl;
	++numTotElectrons;
	if(IsCleanElectron(*ele_i)){++numTotCleanElectrons;
        h2dElectronEtaPt[process]->Fill(ele_i->pt(),ele_i->eta(),weight);
	}
	if(IsCleanElectron(*ele_i)&&IsIsolatedElectron(*ele_i)){
	++numTotIsolatedElectrons;
	//Electron base plot
	hElectronPt[process]->Fill(ele_i->pt(),weight);
	hElectronEta[process]->Fill(ele_i->eta(),weight);
	hElectronPhi[process]->Fill(ele_i->phi(),weight);
        h2dIsoElectronEtaPt[process]->Fill(ele_i->pt(),ele_i->eta(),weight);
	//Electron isolation
        hElectronIso[process]->Fill(ele_i->trackIso(),weight);
	//Invariant mass plots
	//Electron pairs
   	for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
		if(IsCleanElectron(*ele_j)&&IsIsolatedElectron(*ele_j)){
       		if( ele_i->charge()==-1 && ele_j->charge()==+1){
			if (debug) cout <<"Invariant Mass e+ e-: " << (ele_i->p4()+ele_j->p4()).M() <<endl;
			hInvMass[process]->Fill((ele_i->p4()+ele_j->p4()).M(),weight);
			if(n_bJet>=2){
			for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                		for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
				
					if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+ele_i->p4()+ele_j->p4()).M(),weight);}
				}
			}
			}
			hInvMElectron[process]->Fill((ele_i->p4()+ele_j->p4()).M(),weight);
			hInvMSFOS[process]->Fill((ele_i->p4()+ele_j->p4()).M(),weight);
			h2dMETInvMassSFOS[process]->Fill((ele_i->p4()+ele_j->p4()).M(),MET,weight);
			h2dPtJetsInvMassSFOS[process]->Fill((ele_i->p4()+ele_j->p4()).M(),pt4Jets,weight);
		}
		}
	}
	}
   }
  
  //Electron multiplicity 
  hElectronMult[process]->Fill(n_Electrons,weight);
  hLeptonMult[process]->Fill(n_Muons+n_Electrons,weight);
  
  
}

void SUSYDiLeptonAnalysis::MuonTnP(const edm::Handle< std::vector<pat::Muon> >& muons, double weight, const int process){
   for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
	if(IsCleanMuon(*mu_i)&&IsIsolatedMuon(*mu_i)){
	//h2dTnPMuonEtaPt[process]->Fill(mu_i->pt(),mu_i->eta(),weight);
	//h2dTnPIsoMuonEtaPt[process]->Fill(mu_i->pt(),mu_i->eta(),weight);
	for (std::vector<pat::Muon>::const_iterator mu_j = muons->begin(); mu_j != muons->end(); ++mu_j){
                if(IsCleanMuon(*mu_j) && (mu_i->p4()+mu_j->p4()).M()>cut_TnPlow && (mu_i->p4()+mu_j->p4()).M()<cut_TnPhigh ){
			h2dTnPMuonEtaPt[process]->Fill(mu_j->pt(),mu_j->eta(),weight);
	                if(IsIsolatedMuon(*mu_j)) {h2dTnPIsoMuonEtaPt[process]->Fill(mu_j->pt(),mu_j->eta(),weight);}	
                }
        }
	}
   }
}

void SUSYDiLeptonAnalysis::TriggerMonitor(const edm::Handle< edm::TriggerResults>& trigger, double weight, const int process){
   
   //const int nFilters = trigger->size();
   //int nAccept = 0;
   TriggerNames names(*trigger);
   std::vector< std::string > hlNames;
   for (TriggerNames::Strings::const_iterator j = names.triggerNames().begin(); j !=names.triggerNames().end(); ++j ) { 
     	hlNames.push_back(*j);
   }
	
   hTrigger[process]->Fill(hlNames[0].c_str(),weight);
   //LogPrint("Trigger") << "Event" << "\n";
   for (unsigned int i=0; i<trigger->size()-1; ++i){
	hTrigger[process]->Fill(hlNames[i].c_str(),0);
   	if((*trigger)[i].accept()){
		//++nAccept;
   		hTrigger[process]->Fill(hlNames[i].c_str(),weight);
		//LogPrint("Event")  << i << " : " << hlNames[i] ;
	}
   }
   //LogPrint("Event") << "Number of triggers: " << nFilters << "\n"
   //		     << "Number of accepted triggers: " << nAccept;
}

//SUSYDiLeptonAnalysis::analyze
//-----------------------------------------------
//Real analysis
void SUSYDiLeptonAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
   //retrieve the weight from CSA soups
   double weight = 1.;
   int processID = 0;
   bool muonTnP=true;
   if (weighted){
   	edm::Handle< double > weightHandle;
   	iEvent.getByLabel("csaweightproducer","weight", weightHandle);
   	weight = * weightHandle;  
   	
	edm::Handle< int > processHandle;
   	iEvent.getByLabel("csaweightproducer","AlpgenProcessID", processHandle);
	processID = int((*processHandle)/1000);
  	if (debug) {
		cout << "----------- Starting to analyze -------------------" << endl;
       		cout <<" Alpgen Process ID = " << (*processHandle) << endl;
        	cout <<" Weight = " << weight << endl;
   	}
	
   }   
  
   // retrieve the PAT-objects
   //Muons
   edm::Handle< std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);
  
   //Electrons
   edm::Handle< std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);
   
   //Jets
   edm::Handle< std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);
   
   //MET
   edm::Handle< std::vector<pat::MET> > met;
   iEvent.getByLabel(metSrc, met);
   
   edm::Handle< edm::TriggerResults > trigger;
   iEvent.getByLabel("TriggerResults", trigger);
   TriggerMonitor(trigger,weight,processID);   

   if (mcInfo){
   	//MC gen Particle
	edm::Handle< std::vector<reco::GenParticle> > genParticles;
   	iEvent.getByLabel(mcSrc, genParticles);

  	MCAnalysis(muons,electrons,genParticles,weight,processID); 
   }

   if((processID==1||processID==2||processID==3)||(weight<50.)){
   ++numTotEvents;
   if(processID==1){++numTotEvents_W;}
   if(processID==2){++numTotEvents_Z;}
   if(processID==3){++numTotEvents_tt;}
   if(CheckCuts(muons, electrons, jets, met)){
   	++numTotEventsAfterCuts;
   	if(processID==1){++numTotEventsAfterCuts_W;}
   	if(processID==2){++numTotEventsAfterCuts_Z;}
   	if(processID==3){++numTotEventsAfterCuts_tt;}
   	Analysis(muons, electrons, jets, met, weight, processID);
   	if(muonTnP){ MuonTnP(muons, weight, processID);}
   } 
   }
}

void SUSYDiLeptonAnalysis::MCAnalysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<reco::GenParticle> >& genParticles, double weight, const int process){
  int pid = 0;
  float metx = 0;
  float mety = 0;
  for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
        pid = abs(p_i->pdgId());
        if (p_i->status()==1){
                //Muons (13) with status 1
                if (p_i->pt()>cut_MuonPt&&abs(p_i->eta())<cut_MuonEta&&pid==13){
                        hGenMuonPt[process]->Fill(p_i->pt(),weight);
                        hGenMuonEta[process]->Fill(p_i->eta(),weight);
                        h2dGenMuonEtaPt[process]->Fill(p_i->pt(),p_i->eta(),weight);
               }
                //Electron (11) with status 1
                if (p_i->pt()>cut_ElectronPt&&abs(p_i->eta())<cut_ElectronEta&&abs(p_i->pdgId())==11&&p_i->status()==1){
                        hGenElectronPt[process]->Fill(p_i->pt(),weight);
                        hGenElectronEta[process]->Fill(p_i->eta(),weight);
                        h2dGenElectronEtaPt[process]->Fill(p_i->pt(),p_i->eta(),weight);
                }
                if (    pid == 12 || pid == 13 || pid == 14 || pid == 16 || 
                        pid == 1000022 || pid == 2000012 || pid == 2000014 ||
                        pid == 2000016 || pid == 1000039 || pid == 5000039 ||
                        pid == 4000012 || pid == 9900012 || pid == 9900014 ||
                        pid == 9900016 || pid == 39 ){
                        metx += p_i->px();
                        mety += p_i->py();
                }
        }
  }
  hMissingETmc[process]->Fill(sqrt(metx*metx+mety*mety),weight);
}

//SUSYDiLeptonAnalysis::square
//-----------------------------------------------
double SUSYDiLeptonAnalysis::square(double input)
{ 
  return (input*input);
}

/*//SUSYDiLeptonAnalysis::GetInvMass
//-----------------------------------------------
//Computes invariant mass of 2 particles lept1 and lept2
double SUSYDiLeptonAnalysis::GetInvMass(const pat::Particle& lept1, const pat::Particle& lept2)
{
double radicand = square((lept1).energy()+(lept2).energy()) - ( square((lept1).px()+(lept2).px()) + square((lept1).py()+(lept2).py()) + square((lept1).pz()+(lept2).pz()));

  if ( !(radicand < 0) ) return (sqrt(radicand));
  else return -(sqrt(-radicand));
}*/


//SUSYDiLeptonAnalysis::PrintStatics
//-----------------------------------------------
void SUSYDiLeptonAnalysis::PrintStatistics(void)
{ 
int numTotNonIsolatedElectrons = numTotCleanElectrons - numTotIsolatedElectrons;
int numTotNonIsolatedMuons = numTotCleanMuons - numTotIsolatedMuons;
 
edm::LogPrint("Summary")  << "Total number of events processed = " << numTotEvents << "\n"
 			  << "Total number of events accepted = " << numTotEventsAfterCuts << "\n" << "\n"
			  << "Accepted good objects: " << "\n"
 			  << "Total number of electrons = " << numTotElectrons 
				 << " per event = " << (float)numTotElectrons / (float)numTotEventsAfterCuts << "\n"
		          << "   Clean                  = " << numTotCleanElectrons
				 << " per event = " << (float)numTotCleanElectrons / (float)numTotEventsAfterCuts << "\n"
 			  << "   Isolated               = " << numTotIsolatedElectrons 
				 << " per event = " << (float)numTotIsolatedElectrons / (float)numTotEventsAfterCuts << "\n"
			  << "   Non isolated           = " << numTotNonIsolatedElectrons << "\n"
			  << "Total number of muons     = " << numTotMuons
				<< " per event = " << (float)numTotMuons / (float)numTotEventsAfterCuts << "\n"
			  << "   Clean                  = " << numTotCleanMuons 
				<< " per event = " << (float)numTotCleanMuons / (float)numTotEventsAfterCuts << "\n"
			  << "   Isolated               = " << numTotIsolatedMuons
				<< " per event = " << (float)numTotIsolatedMuons / (float)numTotEventsAfterCuts << "\n"
			  << "   Non isolated           = " << numTotNonIsolatedMuons << "\n"
			  << "Total number of jets      = " << numTotJets
			      	<< " per event = " << (float)numTotJets / (float)numTotEventsAfterCuts << "\n"
			  << "   Clean                  = " << numTotCleanJets
			      	<< " per event = " << (float)numTotCleanJets / (float)numTotEventsAfterCuts << "\n";

edm::LogPrint("Consists of")<< "Total number of W+jets events = " << numTotEvents_W << "\n"
 			    << "Total number of W+jets events accepted = " << numTotEventsAfterCuts_W << "\n" << "\n"
			    << "Total number of Z+jets events = " << numTotEvents_Z << "\n"
 			    << "Total number of Z+jets events accepted = " << numTotEventsAfterCuts_Z << "\n" << "\n"
			    << "Total number of tt+jets events = " << numTotEvents_tt << "\n"
 			    << "Total number of tt+jets events accepted = " << numTotEventsAfterCuts_tt << "\n" << "\n";

}

DEFINE_FWK_MODULE(SUSYDiLeptonAnalysis);
