/** \class SUSYDiLeptonAnalysis
 *
 *  
 *  This class is an EDAnalyzer for PAT 
 *  Layer 0 and Layer 1 output
 *
 *  $Date: 2008/18/18 08:07:04 $
 *  $Revision: 1.4 $ for CMSSW 2_2_3
 *
 *  \author: Niklas Mohr -- niklas.mohr@cern.ch
 *  
 */

#include "NiklasMohr/SUSYDiLepton/interface/SUSYDiLeptonAnalysis.h"

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
    mcInfo            = iConfig.getUntrackedParameter<bool>   ("mcInfo",false);

    //Input collections
    mcSrc             = iConfig.getParameter<edm::InputTag> ("mcSource");
    beamSpotSrc        = iConfig.getParameter<edm::InputTag> ("beamSpotSource");
    muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
    electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
    metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
    jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
    trackSrc          = iConfig.getParameter<edm::InputTag> ("trackSource");
 
    //Cuts
    weighted         = iConfig.getUntrackedParameter<bool>   ("CSA_weighted");
    externalWeight   = iConfig.getUntrackedParameter<double>   ("external_Weight",1.);
    Signal_Analysis  = iConfig.getUntrackedParameter<bool>   ("Signal_Analysis",false);
    rej_Cuts         = iConfig.getUntrackedParameter<bool>   ("rej_Cuts");
    rej_JetCut       = iConfig.getUntrackedParameter<bool>   ("rej_JetCut");
    rej_METCut       = iConfig.getUntrackedParameter<bool>   ("rej_METCut");
    rej_LeptonCut    = iConfig.getUntrackedParameter<string>   ("rej_LeptonCut");
    rej_bTagCut      = iConfig.getUntrackedParameter<bool>   ("rej_bTagCut");
    cut_nPos_Elec    = iConfig.getUntrackedParameter<int> ("user_nPos_Electrons");
    cut_nNeg_Elec    = iConfig.getUntrackedParameter<int> ("user_nNeg_Electrons");
    cut_nPos_Muon    = iConfig.getUntrackedParameter<int> ("user_nPos_Muons");
    cut_nNeg_Muon    = iConfig.getUntrackedParameter<int> ("user_nNeg_Muons");
    cut_nJets        = iConfig.getUntrackedParameter<int> ("user_nJets");
    cut_nbJets       = iConfig.getUntrackedParameter<int> ("user_nbJets");
    cut_ptfirstJet   = iConfig.getUntrackedParameter<double> ("user_pt1JetMin");
    cut_ptsecondJet  = iConfig.getUntrackedParameter<double> ("user_pt2JetMin");
    cut_ptthirdJet   = iConfig.getUntrackedParameter<double> ("user_pt3JetMin");
    cut_ptfourthJet  = iConfig.getUntrackedParameter<double> ("user_pt4JetMin");
    cut_MET	   = iConfig.getUntrackedParameter<double> ("user_METMin");
    cut_bTagDiscriminator = iConfig.getUntrackedParameter<double> ("user_bTagDiscriminator");
  
    bJetAlgo  = iConfig.getUntrackedParameter<string> ("user_bJetAlgo");
  
    cut_MuonPt	   = iConfig.getUntrackedParameter<double> ("acc_MuonPt");
    cut_MuonEta	   = iConfig.getUntrackedParameter<double> ("acc_MuonEta");
    cut_MuonChi2   = iConfig.getUntrackedParameter<double> ("acc_MuonChi2");
    cut_Muond0 	   = iConfig.getUntrackedParameter<double> ("acc_Muond0");
    cut_MuonnHits  = iConfig.getUntrackedParameter<double> ("acc_MuonnHits");
    cut_MuonHCALIso = iConfig.getUntrackedParameter<double> ("iso_MuonHCALIso",9999999999.);
    cut_MuonECALIso = iConfig.getUntrackedParameter<double> ("iso_MuonECALIso",9999999999.);
    cut_MuonIso	   = iConfig.getUntrackedParameter<double> ("iso_MuonIso");
  
    ElectronID = iConfig.getUntrackedParameter<string> ("user_ElectronID");
    cut_ElectronPt  = iConfig.getUntrackedParameter<double> ("acc_ElectronPt");
    cut_ElectronEta  = iConfig.getUntrackedParameter<double> ("acc_ElectronEta");
    cut_Electrond0   = iConfig.getUntrackedParameter<double> ("acc_Electrond0");
    cut_ElectronIso  = iConfig.getUntrackedParameter<double> ("iso_ElectronIso");
  
    cut_JetPt  = iConfig.getUntrackedParameter<double> ("acc_JetPt");
    cut_JetEta  = iConfig.getUntrackedParameter<double> ("acc_JetEta");
    cut_JetEMF  = iConfig.getUntrackedParameter<double> ("acc_JetEMF",0.);
    cut_JetHEF  = iConfig.getUntrackedParameter<double> ("acc_JetHEF",0.);
  
    methodTnP     = iConfig.getUntrackedParameter<string>   ("tnp_method");
    cut_TnPChi  = iConfig.getUntrackedParameter<double>   ("tnp_probe_chi2",100);
    cut_TnPnHits  = iConfig.getUntrackedParameter<double>   ("tnp_probe_chi2",2);
    cut_TnPDr  = iConfig.getUntrackedParameter<double> ("tnp_tag_dr");
    cut_TnPTagPt  = iConfig.getUntrackedParameter<double> ("tnp_tag_pt");
    cut_TnPTagEta  = iConfig.getUntrackedParameter<double> ("tnp_tag_eta");
    cut_TnPlowSB1  = iConfig.getUntrackedParameter<double> ("tnp_low_SB1");
    cut_TnPhighSB1  = iConfig.getUntrackedParameter<double> ("tnp_high_SB1");
    cut_TnPlow  = iConfig.getUntrackedParameter<double> ("tnp_low_inv");
    cut_TnPhigh  = iConfig.getUntrackedParameter<double> ("tnp_high_inv");
    cut_TnPlowSB2  = iConfig.getUntrackedParameter<double> ("tnp_low_SB2");
    cut_TnPhighSB2  = iConfig.getUntrackedParameter<double> ("tnp_high_SB2");
    
    electron_Scale  = iConfig.getUntrackedParameter<double>   ("electron_Scale",0.);
    jet_Scale  = iConfig.getUntrackedParameter<double>   ("jet_Scale",0.);

    //initialize global counters
    numTotEvents = 0;
    numTotEventsAfterCuts = 0;

    numTotElectrons = 0;
    numTotCleanElectrons = 0;
    numTotIsolatedElectrons = 0;
  
    numTotMuons = 0;
    numTotCleanMuons = 0;
    numTotIsolatedMuons = 0;
  
    numTotJets = 0;
    numTotCleanJets = 0;

    const int nHistos=7;

    // Create the root file
    edm::Service<TFileService> theFile;

    // book histograms for multiplicities of leptons and jets
    hLeptonMult = new TH1F * [nHistos];
    hElectronMult = new TH1F * [nHistos];
    hMuonMult = new TH1F * [nHistos];
    hJetMult = new TH1F * [nHistos];
    hbJetMult = new TH1F * [nHistos];

    //histograms for lepton isolation cuts
    hElectronIso = new TH1F * [nHistos];
    hElectronTrackIso = new TH1F * [nHistos];
    hElectronCaloIso = new TH1F * [nHistos];
    hMuonIso = new TH1F * [nHistos];
    hMuonTrackIso = new TH1F * [nHistos];
    hMuonCaloIso = new TH1F * [nHistos];
    
    //histograms for the invariant mass of the leptons
    hInvMSFOS = new TH1F * [nHistos];
    hInvMOFOS = new TH1F * [nHistos];
    hInvMass = new TH1F * [nHistos];
    hInvMElectron = new TH1F * [nHistos];
    hInvMElectronSS = new TH1F * [nHistos];
    hInvMMuon = new TH1F * [nHistos];
    hInvMMuonSS = new TH1F * [nHistos];
    hInvMassMC = new TH1F * [nHistos];
  
    //bbll histos
    hInvMbbllSFOS = new TH1F * [nHistos];
    hInvMbbllOFOS = new TH1F * [nHistos];

    //2D histos
    h2dMETInvMassSFOS = new TH2F * [nHistos];
    h2dHTInvMassSFOS = new TH2F * [nHistos];
    h2dEtJetsInvMassSFOS = new TH2F * [nHistos];
    h2dEtJet4InvMassSFOS = new TH2F * [nHistos];
    h2dIsoInvMassSFOS = new TH2F * [nHistos];

    h2dMETInvMassOFOS = new TH2F * [nHistos];
    h2dHTInvMassOFOS = new TH2F * [nHistos];
    h2dEtJetsInvMassOFOS = new TH2F * [nHistos];
    h2dEtJet4InvMassOFOS = new TH2F * [nHistos];
    h2dIsoInvMassOFOS = new TH2F * [nHistos];

    //muon histograms
    hMuonPt = new TH1F * [nHistos];
    hMuon1Pt = new TH1F * [nHistos];
    hMuon2Pt = new TH1F * [nHistos];
    hMuonEta = new TH1F * [nHistos];
    hMuonPhi = new TH1F * [nHistos];
    hGenMuonPt = new TH1F * [nHistos];
    hGenMuonEta = new TH1F * [nHistos];
    hMuonChi2 = new TH1F * [nHistos];
    hMuond0 = new TH1F * [nHistos];
    hMuond0Sig = new TH1F * [nHistos];
    hMuonnHits = new TH1F * [nHistos];
    hMuonIsod0 = new TH2F * [nHistos];

    //electron histograms
    hElectronPt = new TH1F * [nHistos];
    hElectron1Pt = new TH1F * [nHistos];
    hElectron2Pt = new TH1F * [nHistos];
    hElectronEta = new TH1F * [nHistos];
    hElectronPhi = new TH1F * [nHistos];
    hElectrond0 = new TH1F * [nHistos];
    hGenElectronPt = new TH1F * [nHistos];
    hGenElectronEta = new TH1F * [nHistos];
    hElectronEoverP = new TH1F * [nHistos];
    hElectronfBrem = new TH1F * [nHistos];
    hElectronHoverE = new TH1F * [nHistos];
    hElectrondeltaPhiIn = new TH1F * [nHistos];
    hElectrondeltaEtaIn = new TH1F * [nHistos];
    hElectronIsod0 = new TH2F * [nHistos];

    h2dMuonEtaPt = new TH2F * [nHistos];
    h2dMatchedMuonEtaPt = new TH2F * [nHistos];
    h2dGenMuonEtaPt = new TH2F * [nHistos];
    h2dElectronEtaPt = new TH2F * [nHistos];
    h2dMatchedElectronEtaPt = new TH2F * [nHistos];
    h2dGenElectronEtaPt = new TH2F * [nHistos];
  
    //TnP
    h2dTnPProbeSigEtaPt = new TH2F * [nHistos];
    h2dTnPPassSigEtaPt = new TH2F * [nHistos];
    h2dTnPProbeSSSigEtaPt = new TH2F * [nHistos];
    h2dTnPPassSSSigEtaPt = new TH2F * [nHistos];
    h2dTnPProbeSB1EtaPt = new TH2F * [nHistos];
    h2dTnPPassSB1EtaPt = new TH2F * [nHistos];
    h2dTnPProbeSSSB1EtaPt = new TH2F * [nHistos];
    h2dTnPPassSSSB1EtaPt = new TH2F * [nHistos];
    h2dTnPProbeSB2EtaPt = new TH2F * [nHistos];
    h2dTnPPassSB2EtaPt = new TH2F * [nHistos];
    h2dTnPProbeSSSB2EtaPt = new TH2F * [nHistos];
    h2dTnPPassSSSB2EtaPt = new TH2F * [nHistos];
    hTnPProbeInvMass = new TH1F * [nHistos];
    hTnPPassInvMass = new TH1F * [nHistos];
    hTnPSSInvMass = new TH1F * [nHistos];
 
    //histograms for Missing ET
    hMissingET = new TH1F * [nHistos];
    hMissingETmc =  new TH1F * [nHistos];
    hEtSum = new TH1F * [nHistos];
    halphaT = new TH1F * [nHistos];
    hHT = new TH1F * [nHistos];
    hMHT = new TH1F * [nHistos];
    h2dMETEtSumJets = new TH2F * [nHistos];
    h2dMETHT = new TH2F * [nHistos];
    h2dMETMHT = new TH2F * [nHistos];

    //histograms for jets
    hJetEt = new TH1F * [nHistos];
    hJetEta = new TH1F * [nHistos];
    hJetPhi = new TH1F * [nHistos];
    hEtJet1 = new TH1F * [nHistos];
    hEtJet2 = new TH1F * [nHistos];
    hEtJet3 = new TH1F * [nHistos];
    hEtJet4 = new TH1F * [nHistos];
  
    hTrigger = new TH1F * [nHistos];
    hWeight = new TH1F * [nHistos];

    TFileDirectory General = theFile->mkdir( "General" );
    TFileDirectory Clean = theFile->mkdir( "Clean" );
    TFileDirectory All = theFile->mkdir( "All" );
    TFileDirectory Effcor = theFile->mkdir( "Efficiency corrected" );
    TFileDirectory Unmatched = theFile->mkdir( "Unmatched" );
    TFileDirectory Promt = theFile->mkdir( "Promt" );
    TFileDirectory Decay = theFile->mkdir( "Decay" );

    //Trees for unbinned maximum likelihood fit
    TFileDirectory Tree = theFile->mkdir( "Trees" );
    treeOFOS = Tree.make<TTree>("OFOS tree", "OFOS tree"); 
    treeOFOS->Branch("inv",&invMOFOS,"invMOFOS/F");
    treeOFOS->Branch("weight",&invweight,"invweight/F");
    treeSFOS = Tree.make<TTree>("SFOS tree", "SFOS tree"); 
    treeSFOS->Branch("inv",&invMSFOS,"invMSFOS/F");
    treeSFOS->Branch("weight",&invweight,"invweight/F");
    treeElec = Tree.make<TTree>("Electron tree", "Electron tree"); 
    treeElec->Branch("inv",&invMElec,"invMElec/F");
    treeElec->Branch("weight",&invweight,"invweight/F");
    treeMuon = Tree.make<TTree>("Muon tree", "Muon tree"); 
    treeMuon->Branch("inv",&invMMuon,"invMMuon/F");
    treeMuon->Branch("weight",&invweight,"invweight/F");

    general = 0;
    clean = 1;
    all = 2;
    effcor = 3;
    unmatched = 4;
    promt = 5;
    decay = 6;

    InitHisto(&General,general);
    InitHisto(&Clean,clean);
    InitHisto(&All,all);
    InitHisto(&Effcor,effcor);
    if (mcInfo){
        InitHisto(&Unmatched,unmatched);
        InitHisto(&Promt,promt);
        InitHisto(&Decay,decay);
    }
   
    //Read the efficiencies from the files 
    ReadEfficiency();
}

//Initialize all histos including their boundaries
void inline SUSYDiLeptonAnalysis::InitHisto(TFileDirectory *theFile, const int process)
{
    //Multiplicity plots
    TFileDirectory Multiplicity = theFile->mkdir("Multiplicity"); 
    hLeptonMult[process] = Multiplicity.make<TH1F>( "LeptonMultiplicity", "Multiplicity of electrons + muons", 15, 0.0, 15.0);
    hElectronMult[process] = Multiplicity.make<TH1F>( "ElectronMultiplicity", "Multiplicity of electrons", 10, 0.0, 10.0);
    hMuonMult[process] = Multiplicity.make<TH1F>( "MuonMultiplicity", "Multiplicity of muons", 10, 0.0, 10.0);
    hJetMult[process] = Multiplicity.make<TH1F>( "JetMultiplicity", "Multiplicity of jets", 30, 0.0, 30.0);
    hbJetMult[process] = Multiplicity.make<TH1F>( "bJetMultiplicity", "Multiplicity of b jets", 15, 0.0, 15.0);

        
    TFileDirectory InvMass = theFile->mkdir("Invariant Mass"); 
    //histograms for the invariant mass of the leptons
    hInvMSFOS[process] = InvMass.make<TH1F>( "Invariant mass of SFOS lepton pairs", "Invariant mass of SFOS lepton pairs", 300, 0, 300);
    hInvMOFOS[process] = InvMass.make<TH1F>( "Invariant mass of OFOS lepton pairs", "Invariant mass of OFOS lepton pairs", 300, 0, 300);
    hInvMass[process] = InvMass.make<TH1F>( "Invariant mass of lepton pairs", "Invariant mass of lepton pairs", 300, 0, 300);
    hInvMElectron[process] = InvMass.make<TH1F>( "Invariant mass of electron pairs", "Invariant mass of electron pairs", 300, 0, 300);
    hInvMElectronSS[process] = InvMass.make<TH1F>( "Invariant mass of same sign electron pairs", "Invariant mass of same sign electron pairs", 300, 0, 300);
    hInvMMuon[process] = InvMass.make<TH1F>( "Invariant mass of muon pairs", "Invariant mass of muon pairs", 300, 0, 300);
    hInvMMuonSS[process] = InvMass.make<TH1F>( "Invariant mass of same sign muon pairs", "Invariant mass of same sign muon pairs", 300, 0, 300);
    hInvMassMC[process] = InvMass.make<TH1F>( "Invariant mass of signal decays", "Invariant mass of signal decays", 300, 0, 300);
  
    //histograms for bbll inv mass
    hInvMbbllSFOS[process] = InvMass.make<TH1F>( "Invariant mass of SFOS bbll", "Invariant mass of SFOS bbll", 1500, 0, 1500);
    hInvMbbllOFOS[process] = InvMass.make<TH1F>( "Invariant mass of OFOS bbll", "Invariant mass of OFOS bbll", 1500, 0, 1500);

    h2dMETInvMassSFOS[process] = InvMass.make<TH2F>( "MET - Invariant mass SFOS", "MET - Invariant mass SFOS", 300, 0., 300., 1000, 0., 1000.0);
    h2dHTInvMassSFOS[process] = InvMass.make<TH2F>( "HT - Invariant mass SFOS", "HT - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJetsInvMassSFOS[process] = InvMass.make<TH2F>( "SumJets - Invariant mass SFOS", "SumJets - Invariant mass SFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJet4InvMassSFOS[process] = InvMass.make<TH2F>( "PtJet4 - Invariant mass SFOS", "PtJet4 - Invariant mass SFOS", 300, 0., 300., 2000, 0., 1000.0);
    h2dIsoInvMassSFOS[process] = InvMass.make<TH2F>( "Iso - Invariant mass SFOS", "Iso - Invariant mass SFOS", 300, 0., 300., 600, 0., 6.0);

    h2dMETInvMassOFOS[process] = InvMass.make<TH2F>( "MET - Invariant mass OFOS", "MET - Invariant mass OFOS", 300, 0., 300., 1000, 0., 1000.0);
    h2dHTInvMassOFOS[process] = InvMass.make<TH2F>( "HT - Invariant mass OFOS", "HT - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJetsInvMassOFOS[process] = InvMass.make<TH2F>( "SumJets - Invariant mass OFOS", "SumJets - Invariant mass OFOS", 300, 0., 300., 2000, 0., 2000.0);
    h2dEtJet4InvMassOFOS[process] = InvMass.make<TH2F>( "PtJet4 - Invariant mass OFOS", "PtJet4 - Invariant mass OFOS", 300, 0., 300., 2000, 0., 1000.0);
    h2dIsoInvMassOFOS[process] = InvMass.make<TH2F>( "Iso - Invariant mass OFOS", "Iso - Invariant mass OFOS", 300, 0., 300., 600, 0., 6.0);

    TFileDirectory Muons = theFile->mkdir("Muons"); 
    //muon histograms
    hMuonPt[process] = Muons.make<TH1F>( "muon pt", "muon pt", 1000, 0.0, 1000.0);
    hMuon1Pt[process] = Muons.make<TH1F>( "muon 1 pt", "pt of first muon", 1000, 0.0, 1000.0);
    hMuon2Pt[process] = Muons.make<TH1F>( "muon 2 pt", "pt of second muon", 1000, 0.0, 1000.0);
    hMuonEta[process] = Muons.make<TH1F>( "muon eta", "muon eta", 250, -2.5, 2.5);
    hMuonPhi[process] = Muons.make<TH1F>( "muon phi", "muon phi", 350, -3.5, 3.5);
    //Muon quality variables
    hMuonChi2[process] = Muons.make<TH1F>( "muon track chi2 / dof", "muon track chi2 / dof", 200, 0.0, 20.0);
    hMuond0[process] = Muons.make<TH1F>( "muon track d0", "muon track d0", 400, -0.2, 0.2);
    hMuond0Sig[process] = Muons.make<TH1F>( "muon track d0 significance", "muon track d0 significance", 1000, 0.0, 100.0);
    hMuonnHits[process] = Muons.make<TH1F>( "muon track hits", "muon track number of valid hits", 50, 0.0, 50.0);
    hMuonIsod0[process] = Muons.make<TH2F>( "muon iso d0", "muon iso d0", 300, 0.0 , 3.0, 200, 0.0, 0.2);
    //histograms for lepton isolation cuts
    hMuonIso[process] = Muons.make<TH1F>( "MuonIso", "Isolation of muons", 300, 0.0, 3.0);
    hMuonTrackIso[process] = Muons.make<TH1F>( "MuonTrackIso", "Isolation of muons in tracker", 1000, 0.0, 10.0);
    hMuonCaloIso[process] = Muons.make<TH1F>( "MuonCaloIso", "Isolation of muons in calorimeter", 1000, 0.0, 10.0);
    hGenMuonPt[process] = Muons.make<TH1F>( "Generator muon pt", "Generator muon pt", 1000, 0.0, 1000.0);
    hGenMuonEta[process] = Muons.make<TH1F>( "Generator muon eta", "Generator muon eta", 250, -2.5, 2.5);

    TFileDirectory Electrons = theFile->mkdir("Electrons"); 
    //electron histograms
    hElectronPt[process] = Electrons.make<TH1F>( "electron pt", "electron pt", 1000, 0.0, 1000.0);
    hElectron1Pt[process] = Electrons.make<TH1F>( "electron 1 pt", "pt of first electron", 1000, 0.0, 1000.0);
    hElectron2Pt[process] = Electrons.make<TH1F>( "electron 2 pt", "pt of second electron", 1000, 0.0, 1000.0);
    hElectronEta[process] = Electrons.make<TH1F>( "electron eta", "electron eta", 250, -2.5, 2.5);
    hElectronPhi[process] = Electrons.make<TH1F>( "electron phi", "electron phi", 350, -3.5, 3.5);
    hElectrond0[process] = Electrons.make<TH1F>( "electron track d0", "electron track d0", 400, -0.2, 0.2);
    //histograms for lepton isolation cuts
    hElectronIso[process] = Electrons.make<TH1F>( "ElectronIso", "Isolation of electrons", 300, 0.0, 3.0);
    hElectronTrackIso[process] = Electrons.make<TH1F>( "ElectronTrackIso", "Isolation of electrons in tracker", 1000, 0.0, 10.0); 
    hElectronCaloIso[process] = Electrons.make<TH1F>( "ElectronCaloIso", "Isolation of electrons in calorimeter", 1000, 0.0, 10.0); 
    hGenElectronPt[process] = Electrons.make<TH1F>( "Generator electron pt", "Generator electron pt", 1000, 0.0, 1000.0);
    hGenElectronEta[process] = Electrons.make<TH1F>( "Generator electron eta", "Generator electron eta", 250, -2.5, 2.5);
    //electron variables
    hElectronEoverP[process] = Electrons.make<TH1F>( "electron E over P", "electron E over P", 250, 0.0, 2.5);
    hElectronfBrem[process] = Electrons.make<TH1F>( "electron fBrem", "electron fBrem", 110, 0.0, 1.1);
    hElectronHoverE[process] = Electrons.make<TH1F>( "electron H over E", "electron H over E", 200, 0.0, 0.2);
    hElectrondeltaPhiIn[process] = Electrons.make<TH1F>( "electron deltaPhiIn", "electron deltaPhiIn", 200, -0.1, 0.1);
    hElectrondeltaEtaIn[process] = Electrons.make<TH1F>( "electron deltaEtaIn", "electron deltaEtaIn", 400, -0.2, 0.2);
    hElectronIsod0[process] = Electrons.make<TH2F>( "electron iso d0", "electron iso d0", 300, 0.0 , 3.0, 200, 0.0, 0.2);

    TFileDirectory TnP = theFile->mkdir("TnP"); 
    h2dMuonEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of muons", "Eta - Pt of muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dMatchedMuonEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of matched muons", "Eta - Pt of matched muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dGenMuonEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of MC muons", "Eta - Pt of generator muons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dElectronEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of electrons", "Eta - Pt of electrons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dMatchedElectronEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of matched electrons", "Eta - Pt of matched electrons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dGenElectronEtaPt[process] = TnP.make<TH2F>( "Eta-Pt of MC electrons", "Eta - Pt of generator electrons", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
  
    //TnP
    h2dTnPProbeSigEtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP signal probes", "Eta - Pt of TnP signal probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSigEtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP signal passed probes", "Eta - Pt of TnP signal passed probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPProbeSSSigEtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP signal SS probes", "Eta - Pt of TnP signal SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSSSigEtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP signal passed SS probes", "Eta - Pt of TnP signal passed SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPProbeSB1EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB1 probes", "Eta - Pt of TnP SB1 probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSB1EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB1 passed probes", "Eta - Pt of TnP SB1 passed probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPProbeSSSB1EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB1 SS probes", "Eta - Pt of TnP SB1 SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSSSB1EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB1 passed SS probes", "Eta - Pt of TnP SB1 passed SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPProbeSB2EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB2 probes", "Eta - Pt of TnP SB2 probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSB2EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB2 passed probes", "Eta - Pt of TnP SB2 passed probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPProbeSSSB2EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB2 SS probes", "Eta - Pt of TnP SB2 SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    h2dTnPPassSSSB2EtaPt[process] = TnP.make<TH2F>( "Eta-Pt TnP SB2 passed SS probes", "Eta - Pt of TnP SB2 passed SS probes", 1000, 0.0, 1000.0, 250, -2.5, 2.5); 
    hTnPProbeInvMass[process] = TnP.make<TH1F>( "Invariant mass of probes", "Invariant mass of probes", 300, 0, 300);
    hTnPPassInvMass[process] = TnP.make<TH1F>( "Invariant mass of passed probes", "Invariant mass of passed probes", 300, 0, 300);
    hTnPSSInvMass[process] = TnP.make<TH1F>( "Invariant mass of SS probes", "Invariant mass of SS probes", 300, 0, 300);
 
    //histograms for Missing ET
    TFileDirectory MET = theFile->mkdir("MET"); 
    hMissingET[process] = MET.make<TH1F>( "Missing transverse energy", "Missing transverse energy", 1000, 0.0, 1000.0);
    hMissingETmc[process] =  MET.make<TH1F>( "Missing transverse energy MC", "Missing transverse energy MC", 1000, 0.0, 1000.0);
    halphaT[process] = MET.make<TH1F>( "alphaT", "alphaT", 100, 0.0, 1.0);
    hEtSum[process] = MET.make<TH1F>( "ETsum", "Transverse energy sum ET", 2000, 0.0, 2000.0);
    hHT[process] = MET.make<TH1F>( "HT", "Transverse energy sum HT", 2000, 0.0, 2000.0);
    hMHT[process] = MET.make<TH1F>( "MHT", "Transverse energy sum HT + MET", 2000, 0.0, 2000.0);
    h2dMETEtSumJets[process] = MET.make<TH2F>( "MET - sum4Jets", "MET - sum4Jets", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    h2dMETHT[process] = MET.make<TH2F>( "MET - HT", "MET - HT", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    h2dMETMHT[process] = MET.make<TH2F>( "MET - MHT", "MET - MHT", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);

    //histograms for jets
    TFileDirectory Jets = theFile->mkdir("Jets"); 
    hJetEt[process] = Jets.make<TH1F>( "jet pt", "jet pt", 1500, 0.0, 1500.0);
    hJetEta[process] = Jets.make<TH1F>( "jet eta", "jet eta", 300, -3., 3.);
    hJetPhi[process] = Jets.make<TH1F>( "jet phi", "jet phi", 350, -3.5, 3.5);
    hEtJet1[process] = Jets.make<TH1F>( "Pt first jet", "Pt spectrum of the 1st jet", 1500, 0.0, 1500.0);
    hEtJet2[process] = Jets.make<TH1F>( "Pt second jet", "Pt spectrum of the 2nd jet", 1500, 0.0, 1500.0);
    hEtJet3[process] = Jets.make<TH1F>( "Pt third jet", "Pt spectrum of the 3rd jet", 1500, 0.0, 1500.0);
    hEtJet4[process] = Jets.make<TH1F>( "Pt fourth jet", "Pt spectrum of the 4th jet", 1500, 0.0, 1500.0);
 
    hTrigger[process] = theFile->make<TH1F>( "Trigger paths", "Trigger paths", 160, 0, 160);
    hWeight[process] = theFile->make<TH1F>( "Weigths", "Weights", 1000, 0, 1000);
}

//Destructor
SUSYDiLeptonAnalysis::~SUSYDiLeptonAnalysis()
{
    PrintStatistics();
    if (debug) cout << "************* Finished analysis" << endl;
} 

void SUSYDiLeptonAnalysis::ReadEfficiency(){
    std::string muon_fname = "NiklasMohr/SUSYDiLepton/data/efficiency_Muon.data"; 
    std::string electron_fname = "NiklasMohr/SUSYDiLepton/data/efficiency_Electron.data";
    edm::LogPrint("Summary")  << "Reading muon efficiencies from file " << muon_fname << "\n"; 
    std::ifstream muon_listfile(edm::FileInPath(muon_fname).fullPath().c_str());
    edm::LogPrint("Summary")  << "Reading electron efficiencies from file " << electron_fname << "\n";
    std::ifstream electron_listfile(edm::FileInPath(electron_fname).fullPath().c_str());
   
    float muon_eff=0.;
    float electron_eff=0.;
    int muon_nEta;
    int electron_nEta;
    int muon_nPt;
    int electron_nPt;
    muon_listfile >> muon_nEta;
    muon_listfile >> muon_nPt;
    electron_listfile >> electron_nEta;
    electron_listfile >> electron_nPt;
    if (muon_nEta != nEtaBins && electron_nEta != nEtaBins && muon_nPt != nPtBins && electron_nPt != nPtBins) { 
        std::cout << " *** ERROR -> Check data files : nEta bins " 
        << muon_nEta << " instead of " << nEtaBins << std::endl;
    }
    else{
        for (int i=0; i<nEtaBins; i++) {
            for (int j=0; j<nPtBins; j++) {
                muon_listfile >> muon_eff;
                electron_listfile >> electron_eff;
                Muon_Eff[i][j]=muon_eff;
                Electron_Eff[i][j]=electron_eff;
                if (debug) cout << "Cat: " << i << "," << j << " muon: " << Muon_Eff[i][j] << " electron: " << Electron_Eff[i][j] << endl; 
            }
        }
    }
}    

//Check if Muon is clean and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanMuon(const pat::Muon& muon)
{
    if(!muon.innerTrack().isNull()&&muon.isGlobalMuon()){
        if (muon.pt()>=cut_MuonPt &&
                abs(muon.eta())<=cut_MuonEta &&
                muon.globalTrack()->normalizedChi2()<cut_MuonChi2 &&
                muon.innerTrack()->numberOfValidHits()>=cut_MuonnHits &&
                abs(muon.innerTrack()->dxy(bs))<cut_Muond0 &&
                muon.hcalIsoDeposit()->candEnergy() < cut_MuonHCALIso &&
                muon.ecalIsoDeposit()->candEnergy() < cut_MuonECALIso 
                )
                {return true;} else return false;
    }
    else return false;
}

//Check if Electron is identified and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanElectron(const pat::Electron& electron)
{
    if(!electron.gsfTrack().isNull()){
        if (electron.electronID(ElectronID)==1.0 &&
            electron.pt()>=cut_ElectronPt && 
            abs(electron.eta())<=cut_ElectronEta &&
            abs(electron.gsfTrack()->dxy(bs))<cut_Electrond0
            ) 
            {return true;} else return false;
    }
    else return false;
}

//Check if Jet is not overlapping with electrons and in geometrical acceptance
bool SUSYDiLeptonAnalysis::IsCleanJet(const pat::Jet& jet)
{
    /*if (cut_JetEMF >= 0.){
        if (jet.emEnergyFraction()>cut_JetEMF&&jet.pt()>cut_JetPt&&abs(jet.eta())<cut_JetEta){return true;} else return false;
    }
    else*/
    if (pat::Flags::test(jet, pat::Flags::Overlap::Electrons) &&
            jet.emEnergyFraction()>cut_JetEMF &&
            jet.energyFractionHadronic()>cut_JetHEF &&
            jet.pt()>cut_JetPt &&
            abs(jet.eta())<cut_JetEta ){return true;} else return false;
}

//calculate the isolation of a pat lepton relative sum in cone: (tracker+ecal+hcal)/lept.pt
template < class T > 
double CalcIso(const T & lepton)
{
    //double cut_ConeSize = 0.3;
    //reco::isodeposit::Direction candDir(lepton.eta(), lepton.phi());
    double value = (lepton.trackIso()+lepton.ecalIso()+lepton.hcalIso())/lepton.pt();
    /*double value =  (lepton.trackerIsoDeposit()->depositWithin(cut_ConeSize)+
                    lepton.ecalIsoDeposit()->depositWithin(cut_ConeSize)+
                    lepton.hcalIsoDeposit()->depositWithin(cut_ConeSize))/lepton.pt();*/
    return value;
}

template < class T > 
int GetLeptKind(const T * lepton)
{
    int value = 0;
    if(lepton->genLepton()){
        if(lepton->genLepton()->status()==1 && lepton->genLepton()->numberOfMothers()==1){
            const reco::Candidate * mom = lepton->genLepton()->mother();
            //Check if lepton is promt (itself,tau,Z,W,SUSY)
            if(mom->pdgId()==lepton->genLepton()->pdgId()||abs(mom->pdgId())==15||abs(mom->pdgId())==23||abs(mom->pdgId())==24||abs(mom->pdgId())>1000000){
                value=5;
            }
            else {
                //LogPrint("Lepton") << "Non promt: " << mom->pdgId();
                value=6;
            }
        }
        else value = 6;
    }
    else value=4;
    return value;
}
//Find the combination of objects with the closest distance in pt
bool SUSYDiLeptonAnalysis::FindMinCombo(std::vector<const reco::Candidate*> objects,
				std::vector<unsigned int> & lista,
				std::vector<unsigned int> & listb) {
 
    unsigned int n = objects.size();
    if (n>Combinations::nmax) {LogPrint("alhpaCalc") << "FindMinDptCombo: " << n << " too big"; return false;}

    lista.clear(); listb.clear(); // clears the lists, just to be sure

    double mindiff = 1000000000., diff = 0.;
    // Strip indices from the relevant combination set
    std::vector<unsigned int> la; //!< Temporary list a for calculating best combo
    std::vector<unsigned int> lb; //!< Temporary list b for calculating best combo

    for (unsigned int r=0; r<Combinations::rmax[n]; r++) {
        for (unsigned int j=0; j<Combinations::jmax[n]; j++) {
        //cout << r << " " << j << endl;
        // populate list a
        for (unsigned int ia=0; ia<(n-1); ia++) {
	    if ((n==2) && (Combinations::la2[r][j][ia]!=-1) ) la.push_back( Combinations::la2[r][j][ia] );	
	    if ((n==3) && (Combinations::la3[r][j][ia]!=-1) ) la.push_back( Combinations::la3[r][j][ia] );
	    if ((n==4) && (Combinations::la4[r][j][ia]!=-1) ) la.push_back( Combinations::la4[r][j][ia] );
	    if ((n==5) && (Combinations::la5[r][j][ia]!=-1) ) la.push_back( Combinations::la5[r][j][ia] );
	    if ((n==6) && (Combinations::la6[r][j][ia]!=-1) ) la.push_back( Combinations::la6[r][j][ia] );
	    if ((n==7) && (Combinations::la7[r][j][ia]!=-1) ) la.push_back( Combinations::la7[r][j][ia] );
	//if (lista[ia] < 0) { lista.pop_back(); break; }
        }
        // populate list b
        for (unsigned int ib=0; ib<(n-1); ib++) {
	    if ((n==2) && (Combinations::lb2[r][j][ib]!=-1) ) lb.push_back( Combinations::lb2[r][j][ib] );
	    if ((n==3) && (Combinations::lb3[r][j][ib]!=-1) ) lb.push_back( Combinations::lb3[r][j][ib] );
	    if ((n==4) && (Combinations::lb4[r][j][ib]!=-1) ) lb.push_back( Combinations::lb4[r][j][ib] );
	    if ((n==5) && (Combinations::lb5[r][j][ib]!=-1) ) lb.push_back( Combinations::lb5[r][j][ib] );
	    if ((n==6) && (Combinations::lb6[r][j][ib]!=-1) ) lb.push_back( Combinations::lb6[r][j][ib] );
	    if ((n==7) && (Combinations::lb7[r][j][ib]!=-1) ) lb.push_back( Combinations::lb7[r][j][ib] );
	    //if (listb[ib] < 0) { listb.pop_back(); break; }
        }


        if ( la.size()==0 || lb.size()==0 ) break;

        double JetaEt = 0., JetbEt = 0.;
        for (std::vector<unsigned int>::iterator ia=la.begin();ia!=la.end();++ia) {
	//cout << (*ia) << " ";
	JetaEt += objects[ (*ia) ]->pt();
        }
        //cout << ", ";
        for (std::vector<unsigned int>::iterator ib=lb.begin();ib!=lb.end();++ib) {
	    JetbEt += objects[ (*ib) ]->pt();
	    //cout << (*ib) << " ";
        }
        //cout << endl;
        diff = fabs(JetaEt - JetbEt);
        //cout << "Difference in Et is " << diff << endl;
        if (diff < mindiff) { mindiff = diff; lista = la; listb = lb; }
        la.clear(); lb.clear();
        } // end of permutation list loop
    } // end of combination list loop
  
    return true;
}


double SUSYDiLeptonAnalysis::CalcalphaT(std::vector<const reco::Candidate*> objects,
			  std::vector<unsigned int> la,
			  std::vector<unsigned int> lb) {

    // Check we have enough jets for the lists supplied.
    if ((la[la.size()-1]>objects.size())||(lb[lb.size()-1]>objects.size())) {return 0.;}
    
    double jetaEt = 0., jetbEt = 0.; // sums for Et
    // Loop over jet list a to get pseudo-jet a
    for ( vector<unsigned int>::iterator ia = la.begin(); ia != la.end(); ++ia ) {
      jetaEt += objects[(*ia)]->pt();
    }
      
    // Loop over jet list b to get pseudo-jet b
    for ( vector<unsigned int>::iterator ib = lb.begin(); ib != lb.end(); ++ib ) {
      jetbEt += objects[(*ib)]->pt();
    }
    double ptSum = 0., pxSum = 0., pySum = 0.;
    for (unsigned int i=0; i < objects.size(); ++i){
        ptSum += objects[i]->pt();
        pxSum += objects[i]->px();
        pySum += objects[i]->py();
    }
    double MT = sqrt(ptSum*ptSum-pxSum*pxSum-pySum*pySum);
    double alphaT = min(jetaEt,jetbEt)/MT;
    return alphaT; 
}


//Check if Muon is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedMuon(const pat::Muon& muon)
{
    if (debug) cout << "Muon tracker isolation: "<< CalcIso(muon) << endl;
    if (CalcIso(muon)<cut_MuonIso){return true;} else return false;
}

//Check if Electron is isolated
bool SUSYDiLeptonAnalysis::IsIsolatedElectron(const pat::Electron& electron)
{
    if (debug) cout << "Electron tracker isolation: "<< CalcIso(electron) << endl;
    if (CalcIso(electron)<cut_ElectronIso){return true;} else return false;
}

//Calculate the muon effificiency
double SUSYDiLeptonAnalysis::getMuonWeight(const pat::Muon* muon)
{
    int catEta = 0;
    int catPt = 0;
    for(int i=0; i<nEtaBins; ++i){
        if(boundEta[i]<muon->eta() && muon->eta()<boundEta[i+1]){catEta=i;}
    }
    for(int j=0; j<nPtBins; ++j){
        if(boundPt[j]<muon->pt() && muon->pt()<boundPt[j+1]){catPt=j;}
    }
    float eff = Muon_Eff[catEta][catPt];
    if (debug) cout << "Muon eta: " << muon->eta() << " pt: " << muon->pt() << " category: " << catEta << "," << catPt << " Eff: " << eff << endl;
    if(eff!=0.) return 1/eff;
    else return 0;
}

//Calculate the electron effificiency
double SUSYDiLeptonAnalysis::getElectronWeight(const pat::Electron* electron)
{
    int catEta = 0;
    int catPt = 0;
    for(int i=0; i<nEtaBins; ++i){
        if(boundEta[i]<electron->eta() && electron->eta()<boundEta[i+1]){catEta=i;}
    }
    for(int j=0; j<nPtBins; ++j){
        if(boundPt[j]<electron->pt() && electron->pt()<boundPt[j+1]){catPt=j;}
    }
    float eff = Electron_Eff[catEta][catPt];
    if (debug) cout << "Electron eta: " << electron->eta() << " pt: " << electron->pt() << " category: " << catEta << "," << catPt << " Eff: " << eff << endl;
    if(eff!=0.) return 1/eff;
    else return 0;
}

//check if jets are above cuts
bool SUSYDiLeptonAnalysis::JetCut(const edm::Handle< std::vector<pat::Jet> >& jets){
    if (debug) cout << "Checking Jet cut: " << endl;
    int n_Jet=0;
    int n_bJet=0;
    double ptJets=0;
    bool firstJet=false;
    bool secondJet=false;
    bool thirdJet=false;
    bool fourthJet=false;
    bool bJets=false;
    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
        //Test only clean jets in acceptance
        if(IsCleanJet(*jet_i)){
	    ++n_Jet;
	    ptJets+=jet_i->pt()*(1+jet_Scale);
	    if(n_Jet==1 && jet_i->pt()*(1+jet_Scale)>cut_ptfirstJet){firstJet=true;}
	    if(n_Jet==2 && jet_i->pt()*(1+jet_Scale)>cut_ptsecondJet){secondJet=true;}
	    if(n_Jet==3 && jet_i->pt()*(1+jet_Scale)>cut_ptthirdJet){thirdJet=true;}
	    if(n_Jet==4 && jet_i->pt()*(1+jet_Scale)>cut_ptfourthJet){fourthJet=true;}
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
        if((ptJets>cut_ptfirstJet+cut_ptsecondJet+cut_ptthirdJet+cut_ptfourthJet)&&bJets) {return true;} else return false;
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

//check if MET is above cut
bool SUSYDiLeptonAnalysis::METCut(const edm::Handle< std::vector<pat::MET> >& met, const edm::Handle< std::vector<pat::Jet> >& jets){
    if (debug) cout << "Checking MET cut: " << endl;
    bool MET=false;
    double ptJets = 0.;
    if (jet_Scale!=0.){
	if (debug) cout << "Jet scaling" << endl;
        for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
            if(IsCleanJet(*jet_i)){
	        ptJets+=jet_i->pt();
            }
        }
	if (debug) cout << "ptJets = " << ptJets << endl;
    } 
    for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); ++met_i){      
        if (met_i->et()-jet_Scale*ptJets>cut_MET) {MET=true;} 
    }
    if (MET) {return true;} else return false;
}

//check if two leptons are the event
bool SUSYDiLeptonAnalysis::LeptonCut(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons){
    int n_Pos_Muons = 0, n_Pos_Electrons = 0, n_Neg_Muons = 0, n_Neg_Electrons = 0;
    if (debug) cout << "Checking Two Lepton cut: " << endl;
    for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
        if(IsCleanMuon(*mu_i)&&IsIsolatedMuon(*mu_i)){
            if(mu_i->charge()>0){++n_Pos_Muons;}
            if(mu_i->charge()<0){++n_Neg_Muons;}
        }
    }
    for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
        if(IsCleanElectron(*ele_i)&&IsIsolatedElectron(*ele_i)){
            if(ele_i->charge()>0){++n_Pos_Electrons;}
            if(ele_i->charge()<0){++n_Neg_Electrons;}
        }
    }
    if (rej_LeptonCut == "OS"){
        if (debug) cout << "Checking OS Lepton cut: " << endl;
        if( (n_Pos_Muons+n_Pos_Electrons)>=(cut_nPos_Muon+cut_nPos_Elec) &&
            (n_Neg_Muons+n_Neg_Electrons)>=(cut_nNeg_Muon+cut_nNeg_Elec) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "SS"){
        if (debug) cout << "Checking SS Lepton cut: " << endl;
        if( (n_Pos_Muons+n_Pos_Electrons)>=(cut_nPos_Muon+cut_nPos_Elec) ||
            (n_Neg_Muons+n_Neg_Electrons)>=(cut_nNeg_Muon+cut_nNeg_Elec) ){
                return true;
        } else return false;
    } 
    else if (rej_LeptonCut == "Single"){
        if (debug) cout << "Checking Single Lepton cut: " << endl;
        if( n_Pos_Muons>=cut_nPos_Muon ||
            n_Pos_Electrons >= cut_nPos_Elec ||
            n_Neg_Muons >=cut_nNeg_Muon ||
            n_Neg_Electrons >= cut_nNeg_Elec ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "OSEMU"){
        if (debug) cout << "Checking EMU Lepton cut: " << endl;
        if( (n_Pos_Muons == cut_nPos_Muon && n_Neg_Electrons == cut_nNeg_Elec && n_Neg_Muons == 0 && n_Pos_Electrons == 0) ||
            (n_Neg_Muons == cut_nNeg_Muon && n_Pos_Electrons == cut_nPos_Elec && n_Pos_Muons == 0 && n_Neg_Electrons == 0) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "OSEE"){
        if (debug) cout << "Checking EE Lepton cut: " << endl;
        if( (n_Pos_Electrons == cut_nPos_Elec && n_Neg_Electrons == cut_nNeg_Elec) &&
            (n_Neg_Muons == cut_nNeg_Muon && n_Pos_Muons == cut_nPos_Muon) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "OSMUMU"){
        if (debug) cout << "Checking  Lepton cut: " << endl;
        if( (n_Neg_Muons == cut_nNeg_Muon && n_Pos_Muons == cut_nPos_Muon) && 
            (n_Pos_Electrons == cut_nPos_Elec && n_Neg_Electrons == cut_nNeg_Elec) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "SSMUMU"){
        if (debug) cout << "Checking  Lepton cut: " << endl;
        if( (n_Neg_Electrons == 0 && n_Pos_Electrons == 0) && ( n_Neg_Muons == 2 || n_Pos_Muons == 2 ) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "SSEE"){
        if (debug) cout << "Checking  Lepton cut: " << endl;
        if( (n_Neg_Muons == 0 && n_Pos_Muons == 0) && ( n_Neg_Electrons == 2 || n_Pos_Electrons == 2 ) ){
                return true;
        } else return false;
    }
    else if (rej_LeptonCut == "SSEMU"){
        if (debug) cout << "Checking EMU Lepton cut: " << endl;
        if( (n_Pos_Muons == 1 && n_Pos_Electrons == 1 && n_Neg_Muons == 0 && n_Neg_Electrons == 0) ||
            (n_Neg_Muons == 1 && n_Neg_Electrons == 1 && n_Pos_Muons == 0 && n_Pos_Electrons == 0) ){
                return true;
        } else return false;
    }
    else if( n_Pos_Muons>=cut_nPos_Muon &&
             n_Pos_Electrons >= cut_nPos_Elec &&
             n_Neg_Muons >=cut_nNeg_Muon &&
             n_Neg_Electrons >= cut_nNeg_Elec ){
                return true;
    } else return false;
}

//Check all cuts
bool SUSYDiLeptonAnalysis::CheckCuts(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle< std::vector<pat::MET> >& met){
    if (rej_Cuts){
        bool passedJetCut=false;
        bool passedMETCut=false;
        bool passedLeptonCut=false;
        if (rej_JetCut){
	    if (JetCut(jets)){
	        if (debug) cout << "Jet cut passed" << endl;
		passedJetCut = true;
   	    } 
        } else passedJetCut = true;
   	if (rej_METCut){
	    if (METCut(met,jets)){
	        if (debug) cout << "MET cut passed" << endl;
        	passedMETCut = true;
	    } 
        } else passedMETCut = true;
   	if (!(rej_LeptonCut=="False")){
   	    if (LeptonCut(muons,electrons)){
	        if (debug) cout << "Lepton cut passed" << endl;
        	passedLeptonCut = true;
   		}
	} else passedLeptonCut = true; 
   	if (passedJetCut&&passedMETCut&&passedLeptonCut){return true;} else return false;
    }
    else return true;  
}

//Filling of all histograms and calculation of kinematics
//main part
void SUSYDiLeptonAnalysis::Analysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<pat::Jet> >& jets, const edm::Handle< std::vector<pat::MET> >& met, double weight, const int process){
    hWeight[process]->Fill(weight,weight);
 
    double MET=0;
    for (std::vector<pat::MET>::const_iterator met_i = met->begin(); met_i != met->end(); ++met_i){
        if (debug) cout <<"MET et = "<< met_i->et() << endl;
        MET = met_i->et();
    }

    hMissingET[process]->Fill(MET,weight);
  
    int n_Jet=0;
    int n_bJet=0;
    std::vector< pat::Jet > bJets;
    float et4Jets=0;
    float etFourthJet=0;
    float HT=0;

    std::vector<const reco::Candidate*> objects;

    for (std::vector<pat::Jet>::const_iterator jet_i = jets->begin(); jet_i != jets->end(); ++jet_i){       
        if (debug) cout <<"Jet et = "<< jet_i->pt() << endl;
        ++numTotJets;
        if(IsCleanJet(*jet_i)){
	    ++numTotCleanJets;
	    ++n_Jet;
            objects.push_back(static_cast<const reco::Candidate*>( &(*jet_i) ));
	    //Plots of leading Jets
	    if(n_Jet==1){
                hEtJet1[process]->Fill(jet_i->pt(),weight);
	        et4Jets+=jet_i->pt();}
	    if(n_Jet==2){
                hEtJet2[process]->Fill(jet_i->pt(),weight);
		et4Jets+=jet_i->pt();}
	    if(n_Jet==3){
                hEtJet3[process]->Fill(jet_i->pt(),weight);
		et4Jets+=jet_i->pt();}
	    if(n_Jet==4){
                hEtJet4[process]->Fill(jet_i->pt(),weight);
		et4Jets+=jet_i->pt();
		etFourthJet=jet_i->pt();}
	    //Jet base plots
	    if(jet_i->bDiscriminator(bJetAlgo)>cut_bTagDiscriminator){
	        ++n_bJet;
		bJets.push_back( *jet_i );
	    }
	    HT += jet_i->pt();
	    hJetEt[process]->Fill(jet_i->pt(),weight);
	    hJetEta[process]->Fill(jet_i->eta(),weight);
	    hJetPhi[process]->Fill(jet_i->phi(),weight);
	}
    }

    hJetMult[process]->Fill(n_Jet,weight);
    hbJetMult[process]->Fill(n_bJet,weight);
    if (debug) cout <<" Number of Jets = " << n_Jet << endl;

    //Sum of four leading jets
    hEtSum[process]->Fill(et4Jets,weight);
    h2dMETEtSumJets[process]->Fill(et4Jets,MET,weight);
    //Sum of jet et
    hHT[process]->Fill(HT,weight);
    h2dMETHT[process]->Fill(HT,MET,weight);
    //Sum of jet et + MET
    hMHT[process]->Fill(HT+MET,weight);
    h2dMETMHT[process]->Fill(HT+MET,MET,weight);
    
    //Inv mass variables
    invMOFOS = 0; 
    invMSFOS = 0; 
    invMMuon = 0; 
    invMElec = 0; 
    invweight = 1.; 
    double inv = 0;
    double weightcorr = 0;
    double iso = 0;
    //Muon histograms
    int n_Muons = 0;
    int n_CleanMuons = 0;
    int n_CleanIsoMuons = 0;
    for (std::vector<pat::Muon>::const_iterator mu_i = muons->begin(); mu_i != muons->end(); ++mu_i){
        if (debug) cout <<"mu eta = "<< mu_i->eta() << endl;
	//All muons
	++numTotMuons;
        ++n_Muons;
   	MuonMonitor(&(*mu_i),n_Muons,weight,all); 
	//Clean muons	
	if(IsCleanMuon(*mu_i)){
	    ++numTotCleanMuons;
            ++n_CleanMuons;
   	    MuonMonitor(&(*mu_i),n_CleanMuons,weight,clean);
            if(mcInfo){MuonMonitor(&(*mu_i),n_CleanMuons,weight,GetLeptKind(&(*mu_i)));}   
	    //Clean and isolated muons
	    if(IsIsolatedMuon(*mu_i)){
                objects.push_back(static_cast<const reco::Candidate*>( &(*mu_i) ));
	        ++numTotIsolatedMuons;
	        ++n_CleanIsoMuons;
   	        MuonMonitor(&(*mu_i),n_CleanIsoMuons,weight,general);   
   	        MuonMonitor(&(*mu_i),n_CleanIsoMuons,weight,effcor); 
            }
        }
	//Invariant mass plots 
	//Muon pairs
   	for (std::vector<pat::Muon>::const_iterator mu_j = muons->begin(); mu_j != muons->end(); ++mu_j){
       	    if( mu_i->charge()==+1 && mu_j->charge()==-1){
	        if (debug) cout <<"Invariant Mass mu+ mu-: " << (mu_i->p4()+mu_j->p4()).M() <<endl;
                inv = (mu_i->p4()+mu_j->p4()).M();
                iso = CalcIso(*mu_i)+CalcIso(*mu_j);
                weightcorr = getMuonWeight(&(*mu_i))*getMuonWeight(&(*mu_j))*weight;
                MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,all);
                if(IsCleanMuon(*mu_i)&&IsCleanMuon(*mu_j)){
                    MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,clean);
                    if(IsIsolatedMuon(*mu_i)&&IsIsolatedMuon(*mu_j)){
                        MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                        MuonInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                        invMSFOS = (mu_i->p4()+mu_j->p4()).M();
                        invMMuon = (mu_i->p4()+mu_j->p4()).M();
                        invweight = weightcorr;
                        treeSFOS->Fill();
                        treeMuon->Fill();
                    }
                }
	        if(n_bJet>=2){
		    for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                        for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
		            if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+mu_j->p4()).M(),weight);}
		        }
		    }
	        }
            }
       	    if( mu_i->charge()==mu_j->charge()&&mu_i<mu_j){hInvMMuonSS[process]->Fill((mu_i->p4()+mu_j->p4()).M(),weight);}
	}
        //Wrong pairings for different flavour subtraction
        for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
       	    if(( mu_i->charge()==+1 && ele_j->charge()==-1)|(mu_i->charge()==-1 && ele_j->charge()==+1)){
	        if (debug) cout <<"Invariant Mass mu+ e-, mu- e+: " << (mu_i->p4()+ele_j->p4()).M() <<endl;
                inv = (mu_i->p4()+ele_j->p4()).M();
                weightcorr = getMuonWeight(&(*mu_i))*getElectronWeight(&(*ele_j))*weight;
                iso = CalcIso(*mu_i)+CalcIso(*ele_j);
                OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,all);
                if(IsCleanMuon(*mu_i)&&IsCleanElectron(*ele_j)){
                    OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,clean);
                    if(IsIsolatedMuon(*mu_i)&&IsIsolatedElectron(*ele_j)){
                        OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                        OFOSInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                        invMOFOS = (mu_i->p4()+ele_j->p4()).M();
                        invweight = weightcorr;
                        treeOFOS->Fill();
                    }
                }
	        if(n_bJet>=2){
	            for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                        for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
		            if(bjet_i<bjet_j){hInvMbbllOFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+mu_i->p4()+ele_j->p4()).M(),weight);}
		        }
	            }
	        }
            }
        }            
    }

    //Loop over electrons 
    int n_Electrons = 0;
    int n_CleanElectrons = 0;
    int n_CleanIsoElectrons = 0;
    for (std::vector<pat::Electron>::const_iterator ele_i = electrons->begin(); ele_i != electrons->end(); ++ele_i){
        if (debug) cout <<"ele eta = "<< ele_i->eta() << endl;
        ++numTotElectrons;
	++n_Electrons;
   	ElectronMonitor(&(*ele_i),n_Electrons,weight,all); 
        //Clean electrons  
	if(IsCleanElectron(*ele_i)){
            ++numTotCleanElectrons;
	    ++n_CleanElectrons;
   	    ElectronMonitor(&(*ele_i),n_CleanElectrons,weight,clean);   
            if(mcInfo){ElectronMonitor(&(*ele_i),n_CleanElectrons,weight,GetLeptKind(&(*ele_i)));}  
            //Clean and isolated electrons 
	    if(IsIsolatedElectron(*ele_i)){
	        ++numTotIsolatedElectrons;
	        ++n_CleanIsoElectrons;
                objects.push_back(static_cast<const reco::Candidate*>( &(*ele_i) ));
   	        ElectronMonitor(&(*ele_i),n_CleanIsoElectrons,weight,general);   
   	        ElectronMonitor(&(*ele_i),n_CleanIsoElectrons,weight,effcor);
            }
        } 
        //Invariant mass plots
	//Electron pairs
   	for (std::vector<pat::Electron>::const_iterator ele_j = electrons->begin(); ele_j != electrons->end(); ++ele_j){
       	    if( ele_i->charge()==-1 && ele_j->charge()==+1){
	        if (debug) cout <<"Invariant Mass e+ e-: " << (ele_i->p4()+ele_j->p4()).M() <<endl;
                inv = (ele_i->p4()+ele_j->p4()).M();
                if(electron_Scale<0.00001||electron_Scale>0.00001){
                    reco::Particle::LorentzVector pb_i = reco::Particle::LorentzVector(ele_i->px(),ele_i->py(),ele_i->pz(),(ele_i->energy()+ele_i->energy()*electron_Scale));
                    
                    reco::Particle::LorentzVector pb_j = reco::Particle::LorentzVector(ele_j->px(),ele_j->py(),ele_j->pz(),(ele_j->energy()+ele_j->energy()*electron_Scale));
                    inv = (pb_i+pb_j).M();
                }
                weightcorr = getElectronWeight(&(*ele_i))*getElectronWeight(&(*ele_j))*weight;
                iso = CalcIso(*ele_i)+CalcIso(*ele_j);
                ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,all);
	        if(IsCleanElectron(*ele_i)&&IsCleanElectron(*ele_j)){
                    ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,clean);
                    if(IsIsolatedElectron(*ele_i)&&IsIsolatedElectron(*ele_j)){
                        ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weight,general);
                        ElectronInvMonitor(inv,MET,HT,et4Jets,etFourthJet,iso,weightcorr,effcor);
                        invMSFOS = inv;
                        invMElec = inv;
                        invweight = weightcorr;
                        treeSFOS->Fill();
                        treeElec->Fill();
                    }
                }
                if(n_bJet>=2){
		    for(std::vector<pat::Jet>::const_iterator bjet_i = bJets.begin(); bjet_i!=bJets.end(); ++bjet_i){
                	for(std::vector<pat::Jet>::const_iterator bjet_j = bJets.begin(); bjet_j!=bJets.end();++bjet_j){
			    if(bjet_i<bjet_j){hInvMbbllSFOS[process]->Fill((bjet_i->p4()+bjet_j->p4()+ele_i->p4()+ele_j->p4()).M(),weight);}
			}
                    }
		}		
	    }
       	    if( ele_i->charge()==ele_j->charge()&&ele_i<ele_j){hInvMElectronSS[process]->Fill((ele_i->p4()+ele_j->p4()).M(),weight);}
	}
    }
 
    //Global alphaT from all objects
    double alphaT = 0; 
    if(FindMinCombo(objects,mlaMinDpt,mlbMinDpt)){alphaT =  CalcalphaT(objects,mlaMinDpt,mlbMinDpt);}
    halphaT[process]->Fill(alphaT,weight);
    //Muon multiplicity 
    hMuonMult[process]->Fill(n_Muons,weight);
    //Electron multiplicity 
    hElectronMult[all]->Fill(n_Electrons,weight);
    hElectronMult[clean]->Fill(n_CleanElectrons,weight);
    hElectronMult[general]->Fill(n_CleanIsoElectrons,weight);
    //Lepton multiplicity
    hLeptonMult[general]->Fill(n_Muons+n_CleanIsoElectrons,weight);
}

//Fill all muon inv mass related quantities
void SUSYDiLeptonAnalysis::MuonInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,weight);
    hInvMMuon[process]->Fill(inv,weight);
    hInvMSFOS[process]->Fill(inv,weight);
    h2dMETInvMassSFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassSFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassSFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassSFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassSFOS[process]->Fill(inv,iso,weight);
}

//Fill all electron inv mass related quantities
void SUSYDiLeptonAnalysis::ElectronInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,weight);
    hInvMElectron[process]->Fill(inv,weight);
    hInvMSFOS[process]->Fill(inv,weight);
    h2dMETInvMassSFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassSFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassSFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassSFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassSFOS[process]->Fill(inv,iso,weight);
}

//Fill all OFOS inv mass related quantities
void SUSYDiLeptonAnalysis::OFOSInvMonitor(const double inv,const double MET,const double HT,
        const double et4Jets, const double etFourthJet, const double iso, double weight, const int process){
    hInvMass[process]->Fill(inv,-weight);
    hInvMOFOS[process]->Fill(inv,weight);
    h2dMETInvMassOFOS[process]->Fill(inv,MET,weight);
    h2dEtJetsInvMassOFOS[process]->Fill(inv,et4Jets,weight);
    h2dHTInvMassOFOS[process]->Fill(inv,HT,weight);
    h2dEtJet4InvMassOFOS[process]->Fill(inv,etFourthJet,weight);
    h2dIsoInvMassOFOS[process]->Fill(inv,iso,weight);
}

//TnP analysis muons
void SUSYDiLeptonAnalysis::MuonTnP(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<reco::Track> >& tracks, double weight, const int process){
    for (std::vector<pat::Muon>::const_iterator tag_i = muons->begin(); tag_i != muons->end(); ++tag_i){
        if(IsCleanMuon(*tag_i)&&IsIsolatedMuon(*tag_i)&&tag_i->pt()>cut_TnPTagPt&&abs(tag_i->eta()<cut_TnPTagEta)){
	    for (std::vector<reco::Track>::const_iterator pb_j = tracks->begin(); pb_j != tracks->end(); ++pb_j){
	    //for (std::vector<pat::Muon>::const_iterator pb_j = muons->begin(); pb_j != muons->end(); ++pb_j){
                if(pb_j->pt()>cut_MuonPt &&
                abs(pb_j->eta())<cut_MuonEta){
                //&&
                //pb_j->chi2()<cut_TnPChi &&
                //pb_j->numberOfValidHits()>=cut_TnPnHits){
                    int nMatch = 0;
                    double InvMass = 0;
                    for (std::vector<pat::Muon>::const_iterator tag_j = muons->begin(); tag_j != muons->end(); ++tag_j){
                        if(IsCleanMuon(*tag_j)&&IsIsolatedMuon(*tag_j)&&reco::deltaR(tag_j->eta(),tag_j->phi(),pb_j->eta(),pb_j->phi())<cut_TnPDr){++nMatch;}
                    }
                    reco::Particle::LorentzVector pb = reco::Particle::LorentzVector(pb_j->px(),pb_j->py(),pb_j->pz(),pb_j->p());
                    if(nMatch<=1){InvMass = (tag_i->p4()+pb).M();}
                    /*if(IsCleanMuon(*pb_j)&&IsIsolatedMuon(*pb_j)){nMatch=1;}
                    InvMass = (tag_i->p4()+pb_j->p4()).M();*/
                    //edm::LogPrint("TnP")  << "nMatch = " << nMatch << "\n"
                    //                      << "InvMass = " << InvMass << "\n";
                    if(InvMass>cut_TnPlow && InvMass<cut_TnPhigh ){
		        h2dTnPProbeSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
                            h2dTnPProbeSSSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}
                        }        
                    }
                    if(InvMass>cut_TnPlowSB1 && InvMass<cut_TnPhighSB1 ){
		        h2dTnPProbeSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
		            h2dTnPProbeSSSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
                        }
                    }
                    if(InvMass>cut_TnPlowSB2 && InvMass<cut_TnPhighSB2 ){
		        h2dTnPProbeSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
		            h2dTnPProbeSSSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}
                        }        
                    }
                    if(InvMass>cut_TnPlowSB1 && InvMass<cut_TnPhighSB2 ){
                        hTnPProbeInvMass[process]->Fill(InvMass,weight);
	                if(nMatch==1) {hTnPPassInvMass[process]->Fill(InvMass,weight);}
                        if(tag_i->charge()==pb_j->charge()){
                            hTnPSSInvMass[process]->Fill(InvMass,weight);
                        }        
                    }
                }
            }
        }
    }
}

//TnP analysis electrons
void SUSYDiLeptonAnalysis::ElectronTnP(const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<reco::Track> >& tracks, double weight, const int process){
    for (std::vector<pat::Electron>::const_iterator tag_i = electrons->begin(); tag_i != electrons->end(); ++tag_i){
        if(IsCleanElectron(*tag_i)&&IsIsolatedElectron(*tag_i)&&tag_i->pt()>cut_TnPTagPt&&abs(tag_i->eta()<cut_TnPTagEta)){
	    for (std::vector<reco::Track>::const_iterator pb_j = tracks->begin(); pb_j != tracks->end(); ++pb_j){
	    //for (std::vector<pat::Muon>::const_iterator pb_j = muons->begin(); pb_j != muons->end(); ++pb_j){
                if(pb_j->pt()>cut_ElectronPt &&
                abs(pb_j->eta())<cut_ElectronEta){
                //&&
                //pb_j->chi2()<cut_TnPChi &&
                //pb_j->numberOfValidHits()>=cut_TnPnHits){
                    int nMatch = 0;
                    double InvMass = 0;
                    for (std::vector<pat::Electron>::const_iterator tag_j = electrons->begin(); tag_j != electrons->end(); ++tag_j){
                        if(IsCleanElectron(*tag_j)&&IsIsolatedElectron(*tag_j)&&reco::deltaR(tag_j->eta(),tag_j->phi(),pb_j->eta(),pb_j->phi())<cut_TnPDr){++nMatch;}
                    }
                    reco::Particle::LorentzVector pb = reco::Particle::LorentzVector(pb_j->px(),pb_j->py(),pb_j->pz(),pb_j->p());
                    if(nMatch<=1){InvMass = (tag_i->p4()+pb).M();}
                    /*if(IsCleanMuon(*pb_j)&&IsIsolatedMuon(*pb_j)){nMatch=1;}
                    InvMass = (tag_i->p4()+pb_j->p4()).M();*/
                    //edm::LogPrint("TnP")  << "nMatch = " << nMatch << "\n"
                    //                      << "InvMass = " << InvMass << "\n";
                    if(InvMass>cut_TnPlow && InvMass<cut_TnPhigh ){
		        h2dTnPProbeSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
                            h2dTnPProbeSSSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSigEtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}
                        }        
                    }
                    if(InvMass>cut_TnPlowSB1 && InvMass<cut_TnPhighSB1 ){
		        h2dTnPProbeSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
		            h2dTnPProbeSSSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSB1EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
                        }
                    }
                    if(InvMass>cut_TnPlowSB2 && InvMass<cut_TnPhighSB2 ){
		        h2dTnPProbeSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                if(nMatch==1) {h2dTnPPassSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}	
		        if(tag_i->charge()==pb_j->charge()){
		            h2dTnPProbeSSSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);
	                    if(nMatch==1) {h2dTnPPassSSSB2EtaPt[process]->Fill(pb_j->pt(),pb_j->eta(),weight);}
                        }        
                    }
                    if(InvMass>cut_TnPlowSB1 && InvMass<cut_TnPhighSB2 ){
                        hTnPProbeInvMass[process]->Fill(InvMass,weight);
	                if(nMatch==1) {hTnPPassInvMass[process]->Fill(InvMass,weight);}
                        if(tag_i->charge()==pb_j->charge()){
                            hTnPSSInvMass[process]->Fill(InvMass,weight);
                        }        
                    }
                }
            }
        }
    }
}

//Fill all the trigger related quantities
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

//Fill all electron related quantities
void SUSYDiLeptonAnalysis::ElectronMonitor(const pat::Electron* electron,const int n_Electron, double weight, const int process){
    if(process==effcor){weight=getElectronWeight(electron);}
    //Electron base plot
    hElectronPt[process]->Fill(electron->pt(),weight);
    if(n_Electron == 1){hElectron1Pt[process]->Fill(electron->pt(),weight);}
    if(n_Electron == 2){hElectron2Pt[process]->Fill(electron->pt(),weight);}
    hElectronEta[process]->Fill(electron->eta(),weight);
    hElectronPhi[process]->Fill(electron->phi(),weight);
    //Electron isolation
    double IsoValue = CalcIso(*electron);
    hElectronIso[process]->Fill(IsoValue,weight);
    hElectronTrackIso[process]->Fill(electron->trackIso(),weight);
    hElectronCaloIso[process]->Fill(electron->caloIso(),weight);
    h2dElectronEtaPt[process]->Fill(electron->pt(),electron->eta(),weight);
  	
    double eOverP = electron->eSuperClusterOverP();
    //double eSeed = electron->superCluster()->seed()->energy();
    double pin  = electron->trackMomentumAtVtx().R();   
    //double eSeedOverPin = eSeed/pin; 
    double pout = electron->trackMomentumOut().R(); 
    double fBrem = (pin-pout)/pin;
    
      
    double hOverE = electron->hadronicOverEm();
    //double sigmaee = sqrt((shapeRef)->covEtaEta());
    double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
    double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  
    hElectronEoverP[process]->Fill(eOverP,weight); 
    hElectronfBrem[process]->Fill(fBrem,weight); 
    hElectronHoverE[process]->Fill(hOverE,weight); 
    hElectrondeltaPhiIn[process]->Fill(deltaPhiIn,weight); 
    hElectrondeltaEtaIn[process]->Fill(deltaEtaIn,weight);
    if(!electron->gsfTrack().isNull()){
        //impact parameter
        double d0tobs = electron->gsfTrack()->dxy(bs);
        hElectrond0[process]->Fill(d0tobs,weight);
        hElectronIsod0[process]->Fill(IsoValue,abs(d0tobs),weight);
    }
    if (mcInfo){
        if(electron->genLepton()){
            h2dMatchedElectronEtaPt[process]->Fill(electron->pt(),electron->eta(),weight);
        }
    } 
}


//Fill all muon related quantities
void SUSYDiLeptonAnalysis::MuonMonitor(const pat::Muon* muon,const int n_Muon, double weight, const int process){
    if(process==effcor){weight=getMuonWeight(muon);}
    h2dMuonEtaPt[process]->Fill(muon->pt(),muon->eta(),weight);
    //Muon base plots
    hMuonPt[process]->Fill(muon->pt(),weight);
	
    if(n_Muon == 1){hMuon1Pt[process]->Fill(muon->pt(),weight);}
    if(n_Muon == 2){hMuon2Pt[process]->Fill(muon->pt(),weight);}
    hMuonEta[process]->Fill(muon->eta(),weight);
    hMuonPhi[process]->Fill(muon->phi(),weight);
    
    //Muon isolation
    double IsoValue = CalcIso(*muon);
    hMuonIso[process]->Fill(IsoValue,weight);
    hMuonTrackIso[process]->Fill(muon->trackIso(),weight);
    hMuonCaloIso[process]->Fill(muon->caloIso(),weight);

    if(!muon->innerTrack().isNull()){
        double d0tobs = muon->innerTrack()->dxy(bs);
    	hMuonChi2[process]->Fill(muon->innerTrack()->normalizedChi2(),weight);
        hMuond0[process]->Fill(d0tobs,weight);
	hMuond0Sig[process]->Fill(abs(d0tobs)/muon->innerTrack()->dxyError(),weight);
	hMuonnHits[process]->Fill(muon->innerTrack()->numberOfValidHits(),weight);
        hMuonIsod0[process]->Fill(IsoValue,abs(d0tobs),weight);
    }
    if (mcInfo){
        if(muon->genLepton()){
            h2dMatchedMuonEtaPt[process]->Fill(muon->pt(),muon->eta(),weight);
        }
    }
}
 
//Event loop
//Gets all collections and calls Analysis
void SUSYDiLeptonAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    double weight = externalWeight;

    //Beam spot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel(beamSpotSrc, beamSpotHandle);
    bs = beamSpotHandle->position();
  
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
  
    //Trigger 
    edm::Handle< edm::TriggerResults > trigger;
    iEvent.getByLabel("TriggerResults", trigger);

    bool signal = false;
    ++numTotEvents;
    if(CheckCuts(muons, electrons, jets, met)){
        ++numTotEventsAfterCuts;
        if (mcInfo){
            //MC gen Particle
            edm::Handle< std::vector<reco::GenParticle> > genParticles;
            iEvent.getByLabel(mcSrc, genParticles);
            signal = MCAnalysis(muons,electrons,genParticles,weight,general);
        }
        if (signal&&Signal_Analysis){   
            TriggerMonitor(trigger,weight,general);
            Analysis(muons, electrons, jets, met, weight, general);
        }
        if(!Signal_Analysis){
            TriggerMonitor(trigger,weight,general);
            Analysis(muons, electrons, jets, met, weight, general);
        }
        if(methodTnP=="muons"){ 
            //Tracks
            edm::Handle< std::vector<reco::Track> > tracks;
            iEvent.getByLabel(trackSrc, tracks);
            MuonTnP(muons, tracks, weight, general);
        }
        if(methodTnP=="electrons"){ 
            //Tracks
            edm::Handle< std::vector<reco::Track> > tracks;
            iEvent.getByLabel(trackSrc, tracks);
            ElectronTnP(electrons, tracks, weight, general);
        }
    } 
}


//MC analysis of leptons
bool SUSYDiLeptonAnalysis::MCAnalysis(const edm::Handle< std::vector<pat::Muon> >& muons, const edm::Handle< std::vector<pat::Electron> >& electrons, const edm::Handle< std::vector<reco::GenParticle> >& genParticles, double weight, const int process){
    bool signal = false;
    int pid = 0;
    float metx = 0;
    float mety = 0;
    for (std::vector<reco::GenParticle>::const_iterator p_i = genParticles->begin(); p_i != genParticles->end(); ++p_i){
        pid = abs(p_i->pdgId());
        if (pid==1000023){
            //cout << "Found neutralino with 3 daughters" << endl;
            std::vector<const reco::Candidate *> mc_electrons;
            std::vector<const reco::Candidate *> mc_muons;
            std::vector<const reco::Candidate *> mc_chi1;
            for(size_t j = 0; j < p_i->numberOfDaughters(); ++j){
                const reco::Candidate * d = p_i->daughter(j);
                //cout << d->pdgId() << endl;
 	        if (d->pt()>cut_MuonPt&&abs(d->eta())<cut_MuonEta&& abs(d->pdgId())==11){mc_electrons.push_back(d);}     
 	        if (d->pt()>cut_MuonPt&&abs(d->eta())<cut_MuonEta&& abs(d->pdgId())==13){mc_muons.push_back(d);}     
                if( abs(d->pdgId()) == 1000022){mc_chi1.push_back(d);}
            }
            if (mc_chi1.size()==1 && mc_electrons.size()==2){
		hInvMassMC[process]->Fill((mc_electrons[0]->p4()+mc_electrons[1]->p4()).M(),weight);
                signal=true;
            }
            if (mc_chi1.size()==1 && mc_muons.size()==2){
		hInvMassMC[process]->Fill((mc_muons[0]->p4()+mc_muons[1]->p4()).M(),weight);
                signal=true;
            }
        }
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
 	    if ( pid == 12 || pid == 13 || pid == 14 || pid == 16 || 
		 pid == 1000022 || pid == 2000012 || pid == 2000014 ||
		 pid == 2000016 || pid == 1000039 || pid == 5000039 ||
		 pid == 4000012 || pid == 9900012 || pid == 9900014 ||
		 pid == 9900016 || pid == 39 ){
	            metx += p_i->px(); //TODO define met using et
		    mety += p_i->py();
	    }
        }
    }
    hMissingETmc[process]->Fill(sqrt(metx*metx+mety*mety),weight);

    return signal;
}

//square
double SUSYDiLeptonAnalysis::square(double input)
{ 
    return (input*input);
}

//SUSYDiLeptonAnalysis::GetInvMass
//-----------------------------------------------
//Computes invariant mass of 2 particles lept1 and lept2
/*double SUSYDiLeptonAnalysis::GetInvMass(double E1, double E2, double px1, double px2, double py1, double py2, double pz1, double pz4)
{
    double radicand = square(E1+E2) - (square(px1+px1) + square(py1+py2) + square(pz1+py2));
    if ( !(radicand < 0) ) return (sqrt(radicand));
    else return -(sqrt(-radicand));
}*/


//PrintStatics
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
}

DEFINE_FWK_MODULE(SUSYDiLeptonAnalysis);
