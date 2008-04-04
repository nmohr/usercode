/** \class SUSYDiLeptonAnalysis
 *
 *  
 *  This class is an EDAnalyzer for PAT 
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 17.03.2008$
 *  $Revision: 0.1 $ for CMSSW 1_6_10
 *
 *  \author: Niklas Mohr -- niklas.mohr@cern.ch
 *  
 */

#include "SUSYBSMAnalysis/SUSYDiLepton/interface/SUSYDiLeptonAnalysis.h"

using namespace edm;
using namespace reco;
using namespace std;

//Constructor
SUSYDiLeptonAnalysis::SUSYDiLeptonAnalysis(const edm::ParameterSet &iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
 
  weighted         = iConfig.getUntrackedParameter<bool>   ("CSA_weighted");
  rej_Cuts         = iConfig.getUntrackedParameter<bool>   ("rej_Cuts");
  rej_JetCut       = iConfig.getUntrackedParameter<bool>   ("rej_JetCut");
  rej_METCut       = iConfig.getUntrackedParameter<bool>   ("rej_METCut");
  rej_TwoLeptonCut = iConfig.getUntrackedParameter<bool>   ("rej_TwoLeptonCut");
  cut_nJets        = iConfig.getUntrackedParameter<int> ("user_nJets");
  cut_ptfirstJet   = iConfig.getUntrackedParameter<double> ("user_pt1JetMin");
  cut_ptsecondJet  = iConfig.getUntrackedParameter<double> ("user_pt2JetMin");
  cut_ptthirdJet   = iConfig.getUntrackedParameter<double> ("user_pt3JetMin");
  cut_ptfourthJet  = iConfig.getUntrackedParameter<double> ("user_pt4JetMin");
  cut_MET	   = iConfig.getUntrackedParameter<double> ("user_METMin");

  cut_MuonIso	   = iConfig.getUntrackedParameter<double> ("iso_MuonIso");
  cut_ElectronIso	   = iConfig.getUntrackedParameter<double> ("iso_ElectronIso");
  
  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");

   // book histograms for multiplicities of leptons and jets
        hLeptonMult = new TH1D( "LeptonMult", "Multiplicity of leptons", 15, 0.0, 15.0);
        hElectronMult = new TH1D( "ElectronMult", "Multiplicity of electrons", 10, 0.0, 10.0);
        hMuonMult = new TH1D( "MuonMult", "Multiplicity of muons", 10, 0.0, 10.0);
        hJetMult = new TH1D( "JetMult", "Multiplicity of jets", 30, 0.0, 30.0);

        //histograms for lepton isolation cuts
        hCleanElectronEoverP = new TH1D( "ElectronEoverP", "Ecal energy over tracker momentum of electrons", 2500, 0.0, 2.5);
        hCleanMElectronEoverP = new TH1D( "MatchedElectronEoverP", "Ecal energy over tracker momentum of matched electrons", 2500, 0.0, 2.5);
        hCleanElectronEoverPinv = new TH1D( "ElectronEoverPinv", "(Ecal energy over tracker momentum)^-1 of matched electrons", 300, -0.15, 0.15);
        hCleanMElectronEoverPinv = new TH1D( "MatchedElectronEoverPinv", "(Ecal energy over tracker momentum)^-1 of matched electrons", 300, -0.15, 0.15);
        hCleanElectronHoverE = new TH1D( "ElectronHoverE", "Hcal energy over Ecal energy of electrons", 200, 0.0, 0.2);
        hCleanMElectronHoverE = new TH1D( "MatchedElectronHoverE", "Hcal energy over Ecal energy of matched electrons", 200, 0.0, 0.2);
        hCleanMuondpp = new TH1D( "Muondpp", "dp over p of muons", 100, 0.0, 1.0);
        hCleanMMuondpp = new TH1D( "MatchedMuondpp", "dp over p of matched muons", 100, 0.0, 1.0);
        hCleanMuonchi2 = new TH1D( "Muonchi2", "Chi^2 of muons", 1500, 0.0, 15.0);
        hCleanMMuonchi2 = new TH1D( "MatchedMuonchi2", "Chi^2 of matched muons", 1500, 0.0, 15.0);
        hCleanMuonNhits = new TH1D( "MuonNhits", "# Hits of muons", 55, 5.0, 60.0);
        hCleanMMuonNhits = new TH1D( "MatchedMuonNhits", "# Hits of matched muons", 55, 5.0, 60.0);
        
        hElectronIso = new TH1D( "ElectronIso", "Isolation of electrons", 1000, 0.0, 100.0); 
        hMElectronIso = new TH1D( "MatchedElectronIso", "Isolation of matched electrons", 1000, 0.0, 100.0);
        hElectrondeltaR = new TH1D( "ElectrondeltaR", "deltaR of electron pairs", 500, 0.0, 5.0);
        hMuonIso = new TH1D( "MuonIso", "Isolation of muons", 1000, 0.0, 100.0);
        hMMuonIso = new TH1D( "MatchedMuonIso", "Isolation of matched muons", 1000, 0.0, 100.0);
        hMuondeltaR = new TH1D( "MuondeltaR", "deltaR of muon pairs", 500, 0.0, 5.0);
        
        //histograms for the invariant mass of the leptons
        hInvMSFOS = new TH1D( "Invariant mass of SFOS lepton pairs", "Invariant mass of SFOS lepton pairs", 300, 0, 300);
        hInvMOFOS = new TH1D( "Invariant mass of OFOS lepton pairs", "Invariant mass of OFOS lepton pairs", 300, 0, 300);
        hInvMass = new TH1D( "Invariant mass", "Invariant mass of lepton pairs", 300, 0, 300);
        hInvMElectron = new TH1D( "Invariant mass of electron pairs", "Invariant mass of electron pairs", 300, 0, 300);
        hInvMMuon = new TH1D( "Invariant mass of muon pairs", "Invariant mass of muon pairs", 300, 0, 300);

        //muon histograms
        hMuonPt = new TH1D( "muon pt", "muon pt", 1000, 0.0, 1000.0);
        hMuonEta = new TH1D( "muon eta", "muon eta", 100, -3, 3);
        hMuonPhi = new TH1D( "muon phi", "muon phi", 100, -3.5, 3.5);

        //electron histograms
        hElectronPt = new TH1D( "electron pt", "electron pt", 1000, 0.0, 1000.0);
        hElectronEta = new TH1D( "electron eta", "electron eta", 100, -3, 3);
        hElectronPhi = new TH1D( "electron phi", "electron phi", 100, -3.5, 3.5);

        h2dLeptEtaPt = new TH2D( "LeptEtaPt", "Eta - Pt", 1000, 0.0, 1000.0, 100, -3.0, 3.0); 
 
  
        //histograms for Missing ET
        hMissingET = new TH1D( "Missing transverse energy", "Missing transverse energy", 2000, 0.0, 1000.0);
        hMissingETmc =  new TH1D( "Missing transverse energy MC", "Missing transverse energy MC", 2000, 0.0, 1000.0);
        hEtSum = new TH1D( "ETsum", "Transverse energy sum", 2000, 0.0, 2000.0);
        hHT = new TH1D( "HT", "Transverse energy sum", 2000, 0.0, 2000.0);
        h2dMETptSumJets = new TH2D( "MET-sumJets", "MET - Transverse energy sum of jets", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
        h2dMETHT = new TH2D( "MET-HT", "MET - Transverse energy sum", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
        h2dMETEtSum = new TH2D( "MET-ETSum", "MET - effective mass", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);

        //histograms for jets

        hJetPt = new TH1D( "jet pt", "jet pt", 1500, 0.0, 1500.0);
        hJetEta = new TH1D( "jet eta", "jet eta", 100, -3., 3.);
        hJetPhi = new TH1D( "jet phi", "jet phi", 100, -3.5, 3.5);
        hPtJet1 = new TH1D( "pt first jet", "Pt spectrum of 1st jet", 1500, 0.0, 1500.0);
        hPtJet2 = new TH1D( "pt second jet", "Pt spectrum of 2nd jet", 1500, 0.0, 1500.0);
        hPtJet3 = new TH1D( "pt third jet", "Pt spectrum of 3rd jet", 1500, 0.0, 1500.0);
        hPtJet4 = new TH1D( "pt fourth jet", "Pt spectrum of 4th jet", 1500, 0.0, 1500.0);

}

//Destructor
SUSYDiLeptonAnalysis::~SUSYDiLeptonAnalysis()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[SeedQualityAnalysis] Destructor called" << endl;
   theFile->cd();
   theFile->Write();
   
   //Release the memory
        delete       hLeptonMult;
        delete       hElectronMult;
        delete       hMuonMult;
        delete       hJetMult;

        delete       hCleanElectronEoverP;
        delete       hCleanMElectronEoverP;
        delete       hCleanElectronEoverPinv;
        delete       hCleanMElectronEoverPinv;
        delete       hCleanElectronHoverE;
        delete       hCleanMElectronHoverE;
        delete       hCleanMuondpp;
        delete       hCleanMMuondpp;
        delete       hCleanMuonchi2;
        delete       hCleanMMuonchi2;
        delete       hCleanMuonNhits;
        delete       hCleanMMuonNhits;
        delete       hElectronIso;
        delete       hMElectronIso;
        delete       hElectrondeltaR;
        delete       hMuonIso;
        delete       hMMuonIso;
        delete       hMuondeltaR;

        delete       hMissingET;
        delete       hMissingETmc;
        delete       hEtSum;
        delete       hHT;
        delete       h2dMETptSumJets;
        delete      h2dMETHT;
        delete      h2dMETEtSum;

        delete       hJetPt;
        delete       hJetEta;
        delete       hJetPhi;
        delete       hPtJet1;
        delete       hPtJet2;
        delete       hPtJet3;
        delete       hPtJet4;

        delete       hInvMSFOS;
        delete       hInvMOFOS;
        delete       hInvMass;
        delete       hInvMMuon;
        delete       hInvMElectron;

        delete       hMuonPt;
        delete       hMuonEta;
        delete       hMuonPhi;
        
        delete       hElectronPt;
        delete       hElectronEta;
        delete       hElectronPhi;

        delete       h2dLeptEtaPt;
 
   //Close the Root file
   theFile->Close();
   if (debug) cout << "************* Finished writing histograms to file" << endl;

} 

//SUSYDiLeptonAnalysis::IsClean
//-----------------------------------------------
//Check if Lepton is clean
bool SUSYDiLeptonAnalysis::IsClean(const pat::Particle& lept)
{
 return true;
}

//SUSYDiLeptonAnalysis::IsIsolatedMuon
//-----------------------------------------------
//Check if Muon is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedMuon(const pat::Muon& muon)
{
 if (debug) cout << "Muon tracker isolation: "<< muon.trackIso() << endl;
 if (debug) cout << "Muon pt: "<< muon.pt() << endl;
 if (muon.trackIso()/muon.pt()<cut_MuonIso){return true;} else return false;
}

//SUSYDiLeptonAnalysis::IsIsolatedElectron
//-----------------------------------------------
//Check if Electron is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedElectron(const pat::Electron& electron)
{
 if (debug) cout << "Electron tracker isolation: "<< electron.trackIso() << endl;
 if (debug) cout << "Electron pt: "<< electron.pt() << endl;
 return true;
}

//SUSYDiLeptonAnalysis::JetCut
//-----------------------------------------------
//check if four jets above cuts are existent
bool SUSYDiLeptonAnalysis::JetCut(const edm::Handle< std::vector<pat::Jet> > jets){
  if (debug) cout << "Checking Jet cut: " << endl;
  int n_Jet=0;
  bool firstJet=false;
  bool secondJet=false;
  bool thirdJet=false;
  bool fourthJet=false;
  for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); jet_i++){       
	++n_Jet;
	//Plots of leading Jets
	if(n_Jet==1&&(*jet_i).pt()>cut_ptfirstJet){firstJet=true;}
	if(n_Jet==2&&(*jet_i).pt()>cut_ptsecondJet){secondJet=true;}
	if(n_Jet==3&&(*jet_i).pt()>cut_ptthirdJet){thirdJet=true;}
	if(n_Jet==4&&(*jet_i).pt()>cut_ptfourthJet){fourthJet=true;}
   }
   if (debug&&firstJet) cout <<" 1st Jet above cut " << endl;
   if (debug&&secondJet) cout <<" 2nd Jet above cut " << endl;
   if (debug&&thirdJet) cout <<" 3rd Jet above cut " << endl;
   if (debug&&fourthJet) cout <<" 4th Jet above cut " << endl;
   //Check case of 1 Jet
   if (cut_nJets==1){
	if(firstJet) {return true;} else return false;
   }
   //Check case of 2 Jets
   if (cut_nJets==2){
	if(firstJet&&secondJet) {return true;} else return false;
   }
   //Check case of 3 Jets
   if (cut_nJets==3){
	if(firstJet&&secondJet&&thirdJet) {return true;} else return false;
   }
   //Check case of 4 Jets
   if (cut_nJets==4){
	if(firstJet&&secondJet&&thirdJet&&fourthJet) {return true;} else return false;
   }
   if (cut_nJets!=1||cut_nJets!=2||cut_nJets!=3||cut_nJets!=4) {return false;}

}

//SUSYDiLeptonAnalysis::METCut
//-----------------------------------------------
//check if MET above cut is existent
bool SUSYDiLeptonAnalysis::METCut(const edm::Handle< std::vector<pat::MET> > met){
  if (debug) cout << "Checking MET cut: " << endl;
  bool MET=false;
  for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); met_i++){      
	if ((*met_i).pt()>cut_MET) {MET=true;} 
  }
  if (MET) {return true;} else return false;
}

//SUSYDiLeptonAnalysis::TwoLeptonCut
//-----------------------------------------------
//check if Two leptons are the event
bool SUSYDiLeptonAnalysis::TwoLeptonCut(const edm::Handle< std::vector<pat::Muon> > muons, const edm::Handle< std::vector<pat::Electron> > electrons){
  int n_Muons = 0, n_Electrons = 0;
  if (debug) cout << "Checking Two Lepton cut: " << endl;
  for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
	if(IsIsolatedMuon((*mu_i))){++n_Muons;}
  }
  for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
	if(IsIsolatedElectron((*ele_i))){++n_Electrons;}
  }
  if(n_Muons+n_Electrons>=2) {return true;} else return false;
}

//SUSYDiLeptonAnalysis::CheckCuts
//-----------------------------------------------
//Check all cuts
bool SUSYDiLeptonAnalysis::CheckCuts(const edm::Handle< std::vector<pat::Muon> > muons, const edm::Handle< std::vector<pat::Electron> > electrons, const edm::Handle< std::vector<pat::Jet> > jets, const edm::Handle< std::vector<pat::MET> > met){
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
   	if (METCut(met)){
		if (debug) cout << "MET cut passed" << endl;
        	passedMETCut = true; 
        } else passedJetCut = true;
   	if (TwoLeptonCut(muons,electrons)){
		if (debug) cout << "Two Lepton cut passed" << endl;
        	passedTwoLeptonCut = true;
   	} 
   	if (passedJetCut&&passedMETCut&&passedTwoLeptonCut){return true;} else return false;
   }
   else return true;  
}

//SUSYDiLeptonAnalysis::analysis
//-----------------------------------------------
//Filling of all histograms and calculation of kinematics
void SUSYDiLeptonAnalysis::Analysis(const edm::Handle< std::vector<pat::Muon> > muons, const edm::Handle< std::vector<pat::Electron> > electrons, const edm::Handle< std::vector<pat::Jet> > jets, const edm::Handle< std::vector<pat::MET> > met, double weight){
   //Muon histograms
   for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); mu_i++){
       	if (debug) cout <<"mu eta = "<<(*mu_i).eta()<<" mother = "<<(*mu_i).mother()<<endl;
	//Muon base plots
	hMuonPt->Fill((*mu_i).pt(),weight);
	hMuonEta->Fill((*mu_i).eta(),weight);
	hMuonPhi->Fill((*mu_i).phi(),weight);
	//Invariant mass plots
	//Muon pairs
   	for (std::vector<pat::Muon>::const_iterator mu_j = muons->begin(); mu_j != muons->end(); mu_j++){
       		if((*mu_i).charge()==+1 && (*mu_j).charge()==-1){
			if (debug) cout <<"Invariant Mass mu+ mu-: " << GetInvMass((*mu_i),(*mu_j)) <<endl;
			hInvMass->Fill(GetInvMass((*mu_i),(*mu_j)),weight);
			hInvMMuon->Fill(GetInvMass((*mu_i),(*mu_j)),weight);
			hInvMSFOS->Fill(GetInvMass((*mu_i),(*mu_j)),weight);
		}
	}
   	//Wrong pairungs for different flavour subtraction
	for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ele_j++){
       		if(((*mu_i).charge()==+1 && (*ele_j).charge()==-1)|((*mu_i).charge()==-1 && (*ele_j).charge()==+1)){
			if (debug) cout <<"Invariant Mass mu+ e-, mu- e+: " << GetInvMass((*mu_i),(*ele_j)) <<endl;
			hInvMass->Fill(GetInvMass((*mu_i),(*ele_j)),-weight);
			hInvMOFOS->Fill(GetInvMass((*mu_i),(*ele_j)),weight);
		}
	}
   }
   
  for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ele_i++){
       	if (debug) cout <<"ele eta = "<<(*ele_i).eta()<<" mother = "<<(*ele_i).mother()<<endl;
	//Electron base plot
	hElectronPt->Fill((*ele_i).pt(),weight);
	hElectronEta->Fill((*ele_i).eta(),weight);
	hElectronPhi->Fill((*ele_i).phi(),weight);
	//Invariant mass plots
	//Electron pairs
   	for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ele_j++){
       		if((*ele_i).charge()==-1 && (*ele_j).charge()==+1){
			if (debug) cout <<"Invariant Mass e+ e-: " << GetInvMass((*ele_i),(*ele_j)) <<endl;
			hInvMass->Fill(GetInvMass((*ele_i),(*ele_j)),weight);
			hInvMElectron->Fill(GetInvMass((*ele_i),(*ele_j)),weight);
			hInvMSFOS->Fill(GetInvMass((*ele_i),(*ele_j)),weight);
		}
	}
   }
  
  int n_Jet=0;
  for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); jet_i++){       
       	if (debug) cout <<"Jet pt = "<<(*jet_i).pt()<<" mother = "<<(*jet_i).mother()<<endl;
	++n_Jet;
	//Plots of leading Jets
	if(n_Jet==1){hPtJet1->Fill((*jet_i).pt(),weight);}
	if(n_Jet==2){hPtJet2->Fill((*jet_i).pt(),weight);}
	if(n_Jet==3){hPtJet3->Fill((*jet_i).pt(),weight);}
	if(n_Jet==4){hPtJet4->Fill((*jet_i).pt(),weight);}
	//Jet base plots
	hJetPt->Fill((*jet_i).pt(),weight);
	hJetEta->Fill((*jet_i).eta(),weight);
	hJetPhi->Fill((*jet_i).phi(),weight);
   }
   if (debug) cout <<" Number of Jets = " << n_Jet << endl;
  
  for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); met_i++){
       	if (debug) cout <<"MET pt = "<<(*met_i).pt()<<" mother = "<<(*met_i).mother()<<endl;
	hMissingET->Fill((*met_i).pt(),weight);
   }
}

//SUSYDiLeptonAnalysis::analyze
//-----------------------------------------------
//Real analysis
void SUSYDiLeptonAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
   //retrieve the weight from CSA soups
   double weight = 1.;
   if (weighted){
   	edm::Handle< double > weightHandle;
   	iEvent.getByLabel("csaweightproducer","weight", weightHandle);
   	weight = * weightHandle;  
   }   

   // retrieve the reco-objects
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
   if (debug) cout << "----------- Starting to analyze -------------------" << endl;
   if(CheckCuts(muons, electrons, jets, met)){
   	Analysis(muons, electrons, jets, met, weight);
   }
 
}

//SUSYDiLeptonAnalysis::square
//-----------------------------------------------
double SUSYDiLeptonAnalysis::square(double input)
{ 
  return (input*input);
}

//SUSYDiLeptonAnalysis::GetInvMass
//-----------------------------------------------
//Computes invariant mass of 2 particles lept1 and lept2
double SUSYDiLeptonAnalysis::GetInvMass(const pat::Particle& lept1, const pat::Particle& lept2)
{
double radicand = square((lept1).energy()+(lept2).energy()) - ( square((lept1).px()+(lept2).px()) + square((lept1).py()+(lept2).py()) + square((lept1).pz()+(lept2).pz()));

  if ( !(radicand < 0) ) return(sqrt(radicand));
  else return(sqrt(-radicand));
}


DEFINE_FWK_MODULE(SUSYDiLeptonAnalysis);
