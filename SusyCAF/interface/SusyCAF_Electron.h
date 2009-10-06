#ifndef SUSY_CAF_ELECTRON
#define SUSY_CAF_ELECTRON

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <string>

template< typename T >
class SusyCAF_Electron : public edm::EDProducer {
 public: 
  explicit SusyCAF_Electron(const edm::ParameterSet&);
 private: 
  void produce(edm::Event &, const edm::EventSetup & );

  const edm::InputTag inputTag;
  const std::string Prefix,Suffix;
};

template< typename T >
SusyCAF_Electron<T>::SusyCAF_Electron(const edm::ParameterSet& iConfig) :
  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  Prefix(iConfig.getParameter<std::string>("Prefix")),
  Suffix(iConfig.getParameter<std::string>("Suffix"))
{
  produces <std::vector<float> > (  Prefix + "px" + Suffix);
  produces <std::vector<float> > (  Prefix + "py" + Suffix);
  produces <std::vector<float> > (  Prefix + "pz" + Suffix);
  produces <std::vector<float> > (  Prefix + "pt" + Suffix);
  produces <std::vector<float> > (  Prefix + "eta" + Suffix);
  produces <std::vector<float> > (  Prefix + "phi" + Suffix);
  produces <std::vector<float> > (  Prefix + "energy" + Suffix);
  produces <std::vector<float> > (  Prefix + "charge" + Suffix);
  produces <std::vector<float> > (  Prefix + "gsfTracknormalizedChi2" + Suffix);
  produces <std::vector<float> > (  Prefix + "gsfTracknumberOfValidHits" + Suffix);
  produces <std::vector<float> > (  Prefix + "gsfTrackdxy" + Suffix);
  produces <std::vector<float> > (  Prefix + "gsfTrackdxyError" + Suffix);
  produces <std::vector<float> > (  Prefix + "e1x5" + Suffix);
  produces <std::vector<float> > (  Prefix + "e5x5" + Suffix);
  produces <std::vector<float> > (  Prefix + "fbrem" + Suffix);
  produces <std::vector<float> > (  Prefix + "hcalOverEcal" + Suffix);
  produces <std::vector<float> > (  Prefix + "hcalDepth1OverEcal" + Suffix);
  produces <std::vector<float> > (  Prefix + "hcalDepth2OverEcal" + Suffix);
  produces <std::vector<float> > (  Prefix + "eEleClusterOverPout" + Suffix);
  produces <std::vector<float> > (  Prefix + "eSeedClusterOverPout" + Suffix);
  produces <std::vector<float> > (  Prefix + "eSeedClusterOverP" + Suffix);
  produces <std::vector<float> > (  Prefix + "eSuperClusterOverP" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaPhiSuperClusterTrackAtVtx" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaEtaSuperClusterTrackAtVtx" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaPhiSeedClusterTrackAtCalo" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaEtaSeedClusterTrackAtCalo" + Suffix);
  produces <std::vector<float> > (  Prefix + "sigmaEtaEta" + Suffix);
  produces <std::vector<float> > (  Prefix + "sigmaIetaIeta" + Suffix);
  produces <std::vector<float> > (  Prefix + "classification" + Suffix);
  produces <std::vector<float> > (  Prefix + "mva" + Suffix);
  produces <std::vector<float> > (  Prefix + "dr03TkSumPt" + Suffix);
  produces <std::vector<float> > (  Prefix + "dr03EcalRecHitSumEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "dr03HcalTowerSumEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "vx" + Suffix);
  produces <std::vector<float> > (  Prefix + "vy" + Suffix);
  produces <std::vector<float> > (  Prefix + "vz" + Suffix);
  produces <std::vector<float> > (  Prefix + "caloEnergy" + Suffix);
  produces <std::vector<float> > (  Prefix + "et" + Suffix);
  produces <std::vector<float> > (  Prefix + "mass" + Suffix);
  produces <std::vector<float> > (  Prefix + "mt" + Suffix);
  produces <std::vector<float> > (  Prefix + "vertexChi2" + Suffix);
  produces <std::vector<float> > (  Prefix + "vertexNdof" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaEtaEleClusterTrackAtCalo" + Suffix);
  produces <std::vector<float> > (  Prefix + "deltaPhiEleClusterTrackAtCalo" + Suffix);
  produces <std::vector<float> > (  Prefix + "dr03HcalDepth1TowerSumEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "dr03HcalDepth2TowerSumEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "e2x5Max" + Suffix);
  produces <std::vector<float> > (  Prefix + "ecalEnergy" + Suffix);
  produces <std::vector<float> > (  Prefix + "ecalEnergyError" + Suffix);
  produces <std::vector<float> > (  Prefix + "electronMomentumError" + Suffix);
  produces <std::vector<float> > (  Prefix + "numberOfTracks" + Suffix);
  produces <std::vector<int> >   (  Prefix + "numberOfBrems" + Suffix);
  produces <std::vector<float> > (  Prefix + "shFracInnerHits" + Suffix);
//  _hi
//  produces <std::vector<float> > (  Prefix + "test" + Suffix);
}

template< typename T >
void SusyCAF_Electron<T>::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<std::vector<float> >  px   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  py   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  pz   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  pt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eta   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  phi   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  energy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  charge   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  gsfTrack_normalizedChi2   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  gsfTrack_numberOfValidHits   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  gsfTrack_dxy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  gsfTrack_dxyError   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  e1x5   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  e5x5   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  fbrem   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  hcalOverEcal   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  hcalDepth1OverEcal   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  hcalDepth2OverEcal   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eEleClusterOverPout   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eSeedClusterOverPout   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eSeedClusterOverP   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eSuperClusterOverP   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaPhiSuperClusterTrackAtVtx   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaEtaSuperClusterTrackAtVtx   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaPhiSeedClusterTrackAtCalo   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaEtaSeedClusterTrackAtCalo   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  sigmaEtaEta   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  sigmaIetaIeta   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  classification   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  mva   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  dr03TkSumPt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  dr03EcalRecHitSumEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  dr03HcalTowerSumEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vx   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vz   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  caloEnergy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  et   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  mass   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  mt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vertexChi2   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vertexNdof   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaEtaEleClusterTrackAtCalo   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  deltaPhiEleClusterTrackAtCalo   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  dr03HcalDepth1TowerSumEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  dr03HcalDepth2TowerSumEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  e2x5Max   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  ecalEnergy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  ecalEnergyError   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  electronMomentumError   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  numberOfTracks   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<int> >    numberOfBrems   ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<float> >  shFracInnerHits   ( new std::vector<float>()  ) ;
  //_hi
  //std::auto_ptr<std::vector<float> >  test   ( new std::vector<float>()  ) ;
  
  edm::Handle<std::vector<T> > collection;
  iEvent.getByLabel(inputTag,collection);
  
  for(typename std::vector<T>::const_iterator it = collection->begin(); it!=collection->end(); it++) {
    px->push_back(it->px());
    py->push_back(it->py());
    pz->push_back(it->pz());
    pt->push_back(it->pt());
    eta->push_back(it->eta());
    phi->push_back(it->phi());
    energy->push_back(it->energy());
    charge->push_back(it->charge());
    gsfTrack_normalizedChi2->push_back(it->gsfTrack()->normalizedChi2());
    gsfTrack_numberOfValidHits->push_back(it->gsfTrack()->numberOfValidHits());
    gsfTrack_dxy->push_back(it->gsfTrack()->dxy());
    gsfTrack_dxyError->push_back(it->gsfTrack()->dxyError());
    e1x5->push_back(it->e1x5());
    e5x5->push_back(it->e5x5());
    fbrem->push_back(it->fbrem());
    hcalOverEcal->push_back(it->hcalOverEcal());
    hcalDepth1OverEcal->push_back(it->hcalDepth1OverEcal());
    hcalDepth2OverEcal->push_back(it->hcalDepth2OverEcal());
    eEleClusterOverPout->push_back(it->eEleClusterOverPout());
    eSeedClusterOverPout->push_back(it->eSeedClusterOverPout());
    eSeedClusterOverP->push_back(it->eSeedClusterOverP());
    eSuperClusterOverP->push_back(it->eSuperClusterOverP());
    deltaPhiSuperClusterTrackAtVtx->push_back(it->deltaPhiSuperClusterTrackAtVtx());
    deltaEtaSuperClusterTrackAtVtx->push_back(it->deltaEtaSuperClusterTrackAtVtx());
    deltaPhiSeedClusterTrackAtCalo->push_back(it->deltaPhiSeedClusterTrackAtCalo());
    deltaEtaSeedClusterTrackAtCalo->push_back(it->deltaEtaSeedClusterTrackAtCalo());
    sigmaEtaEta->push_back(it->sigmaEtaEta());
    sigmaIetaIeta->push_back(it->sigmaIetaIeta());
    classification->push_back(it->classification());
    mva->push_back(it->mva());
    dr03TkSumPt->push_back(it->dr03TkSumPt());
    dr03EcalRecHitSumEt->push_back(it->dr03EcalRecHitSumEt());
    dr03HcalTowerSumEt->push_back(it->dr03HcalTowerSumEt());
    vx->push_back(it->vx());
    vy->push_back(it->vy());
    vz->push_back(it->vz());
    caloEnergy->push_back(it->caloEnergy());
    et->push_back(it->et());
    mass->push_back(it->mass());
    mt->push_back(it->mt());
    vertexChi2->push_back(it->vertexChi2());
    vertexNdof->push_back(it->vertexNdof());
    deltaEtaEleClusterTrackAtCalo->push_back(it->deltaEtaEleClusterTrackAtCalo());
    deltaPhiEleClusterTrackAtCalo->push_back(it->deltaPhiEleClusterTrackAtCalo());
    dr03HcalDepth1TowerSumEt->push_back(it->dr03HcalDepth1TowerSumEt());
    dr03HcalDepth2TowerSumEt->push_back(it->dr03HcalDepth2TowerSumEt());
    e2x5Max->push_back(it->e2x5Max());
    ecalEnergy->push_back(it->ecalEnergy());
    ecalEnergyError->push_back(it->ecalEnergyError());
    electronMomentumError->push_back(it->electronMomentumError());
    numberOfTracks->push_back(it->numberOfTracks());
    numberOfBrems->push_back(it->numberOfBrems());
    shFracInnerHits->push_back(it->shFracInnerHits());
//  _hi
//  test->push_back(it->test());
  }
  
  iEvent.put( px,  Prefix + "px" + Suffix );
  iEvent.put( py,  Prefix + "py" + Suffix );
  iEvent.put( pz,  Prefix + "pz" + Suffix );
  iEvent.put( pt,  Prefix + "pt" + Suffix );
  iEvent.put( eta,  Prefix + "eta" + Suffix );
  iEvent.put( phi,  Prefix + "phi" + Suffix );
  iEvent.put( energy,  Prefix + "energy" + Suffix );
  iEvent.put( charge,  Prefix + "charge" + Suffix );
  iEvent.put( gsfTrack_normalizedChi2,  Prefix + "gsfTracknormalizedChi2" + Suffix );
  iEvent.put( gsfTrack_numberOfValidHits,  Prefix + "gsfTracknumberOfValidHits" + Suffix );
  iEvent.put( gsfTrack_dxy,  Prefix + "gsfTrackdxy" + Suffix );
  iEvent.put( gsfTrack_dxyError,  Prefix + "gsfTrackdxyError" + Suffix );
  iEvent.put( e1x5,  Prefix + "e1x5" + Suffix );
  iEvent.put( e5x5,  Prefix + "e5x5" + Suffix );
  iEvent.put( fbrem,  Prefix + "fbrem" + Suffix );
  iEvent.put( hcalOverEcal,  Prefix + "hcalOverEcal" + Suffix );
  iEvent.put( hcalDepth1OverEcal,  Prefix + "hcalDepth1OverEcal" + Suffix );
  iEvent.put( hcalDepth2OverEcal,  Prefix + "hcalDepth2OverEcal" + Suffix );
  iEvent.put( eEleClusterOverPout,  Prefix + "eEleClusterOverPout" + Suffix );
  iEvent.put( eSeedClusterOverPout,  Prefix + "eSeedClusterOverPout" + Suffix );
  iEvent.put( eSeedClusterOverP,  Prefix + "eSeedClusterOverP" + Suffix );
  iEvent.put( eSuperClusterOverP,  Prefix + "eSuperClusterOverP" + Suffix );
  iEvent.put( deltaPhiSuperClusterTrackAtVtx,  Prefix + "deltaPhiSuperClusterTrackAtVtx" + Suffix );
  iEvent.put( deltaEtaSuperClusterTrackAtVtx,  Prefix + "deltaEtaSuperClusterTrackAtVtx" + Suffix );
  iEvent.put( deltaPhiSeedClusterTrackAtCalo,  Prefix + "deltaPhiSeedClusterTrackAtCalo" + Suffix );
  iEvent.put( deltaEtaSeedClusterTrackAtCalo,  Prefix + "deltaEtaSeedClusterTrackAtCalo" + Suffix );
  iEvent.put( sigmaEtaEta,  Prefix + "sigmaEtaEta" + Suffix );
  iEvent.put( sigmaIetaIeta,  Prefix + "sigmaIetaIeta" + Suffix );
  iEvent.put( classification,  Prefix + "classification" + Suffix );
  iEvent.put( mva,  Prefix + "mva" + Suffix );
  iEvent.put( dr03TkSumPt,  Prefix + "dr03TkSumPt" + Suffix );
  iEvent.put( dr03EcalRecHitSumEt,  Prefix + "dr03EcalRecHitSumEt" + Suffix );
  iEvent.put( dr03HcalTowerSumEt,  Prefix + "dr03HcalTowerSumEt" + Suffix );
  iEvent.put( vx,  Prefix + "vx" + Suffix );
  iEvent.put( vy,  Prefix + "vy" + Suffix );
  iEvent.put( vz,  Prefix + "vz" + Suffix );
  iEvent.put( caloEnergy,  Prefix + "caloEnergy" + Suffix );
  iEvent.put( et,  Prefix + "et" + Suffix );
  iEvent.put( mass,  Prefix + "mass" + Suffix );
  iEvent.put( mt,  Prefix + "mt" + Suffix );
  iEvent.put( vertexChi2,  Prefix + "vertexChi2" + Suffix );
  iEvent.put( vertexNdof,  Prefix + "vertexNdof" + Suffix );
  iEvent.put( deltaEtaEleClusterTrackAtCalo,  Prefix + "deltaEtaEleClusterTrackAtCalo" + Suffix );
  iEvent.put( deltaPhiEleClusterTrackAtCalo,  Prefix + "deltaPhiEleClusterTrackAtCalo" + Suffix );
  iEvent.put( dr03HcalDepth1TowerSumEt,  Prefix + "dr03HcalDepth1TowerSumEt" + Suffix );
  iEvent.put( dr03HcalDepth2TowerSumEt,  Prefix + "dr03HcalDepth2TowerSumEt" + Suffix );
  iEvent.put( e2x5Max,  Prefix + "e2x5Max" + Suffix );
  iEvent.put( ecalEnergy,  Prefix + "ecalEnergy" + Suffix );
  iEvent.put( ecalEnergyError,  Prefix + "ecalEnergyError" + Suffix );
  iEvent.put( electronMomentumError,  Prefix + "electronMomentumError" + Suffix );
  iEvent.put( numberOfTracks,  Prefix + "numberOfTracks" + Suffix );
  iEvent.put( numberOfBrems,  Prefix + "numberOfBrems" + Suffix );
  iEvent.put( shFracInnerHits,  Prefix + "shFracInnerHits" + Suffix );
//  _hi
//  iEvent.put( test,  Prefix + "test" + Suffix );
}

#endif
