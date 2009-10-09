#ifndef SUSY_CAF_MUON
#define SUSY_CAF_MUON

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <string>

template< typename T >
class SusyCAF_Muon : public edm::EDProducer {
 public: 
  explicit SusyCAF_Muon(const edm::ParameterSet&);
 private: 
  void produce(edm::Event &, const edm::EventSetup & );

  const edm::InputTag inputTag;
  const std::string Prefix,Suffix;
};

template< typename T >
SusyCAF_Muon<T>::SusyCAF_Muon(const edm::ParameterSet& iConfig) :
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
  produces <std::vector<float> > (  Prefix + "et" + Suffix);
  produces <std::vector<float> > (  Prefix + "mass" + Suffix);
  produces <std::vector<float> > (  Prefix + "mt" + Suffix);
  produces <std::vector<float> > (  Prefix + "globalTracknormalizedChi2" + Suffix);
  produces <std::vector<float> > (  Prefix + "globalTracknumberOfValidHits" + Suffix);
  produces <std::vector<float> > (  Prefix + "globalTrackdxy" + Suffix);
  produces <std::vector<float> > (  Prefix + "globalTrackdxyError" + Suffix);
  produces <std::vector<float> > (  Prefix + "isolationR03sumPt" + Suffix);
  produces <std::vector<float> > (  Prefix + "isolationR03emEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "isolationR03hadEt" + Suffix);
  produces <std::vector<float> > (  Prefix + "vx" + Suffix);
  produces <std::vector<float> > (  Prefix + "vy" + Suffix);
  produces <std::vector<float> > (  Prefix + "vz" + Suffix);
  produces <std::vector<float> > (  Prefix + "vertexChi2" + Suffix);
  produces <std::vector<float> > (  Prefix + "vertexNdof" + Suffix);
//  _hi
//  produces <std::vector<float> > (  Prefix + "test" + Suffix);
}

template< typename T >
void SusyCAF_Muon<T>::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<std::vector<float> >  px   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  py   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  pz   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  pt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  eta   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  phi   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  energy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  charge   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  et   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  mass   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  mt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  globalTrack_normalizedChi2   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  globalTrack_numberOfValidHits   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  globalTrack_dxy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  globalTrack_dxyError   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  isolationR03sumPt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  isolationR03emEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  isolationR03hadEt   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vx   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vz   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vertexChi2   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >  vertexNdof   ( new std::vector<float>()  ) ;
  //_hi
  //std::auto_ptr<std::vector<float> >  test   ( new std::vector<float>()  ) ;
  
  edm::Handle<std::vector<T> > collection;
  iEvent.getByLabel(inputTag,collection);
  
  for(typename std::vector<T>::const_iterator it = collection->begin(); it!=collection->end(); it++) {
    if( !it->isGlobalMuon()) continue;
    px->push_back(it->px());
    py->push_back(it->py());
    pz->push_back(it->pz());
    pt->push_back(it->pt());
    eta->push_back(it->eta());
    phi->push_back(it->phi());
    energy->push_back(it->energy());
    charge->push_back(it->charge());
    et->push_back(it->et());
    mass->push_back(it->mass());
    mt->push_back(it->mt());
    globalTrack_normalizedChi2->push_back(it->globalTrack()->normalizedChi2());
    globalTrack_numberOfValidHits->push_back(it->globalTrack()->numberOfValidHits());
    globalTrack_dxy->push_back(it->globalTrack()->dxy());
    globalTrack_dxyError->push_back(it->globalTrack()->dxyError());
    isolationR03sumPt->push_back(it->isolationR03().sumPt);
    isolationR03emEt->push_back(it->isolationR03().emEt);
    isolationR03hadEt->push_back(it->isolationR03().hadEt);
    vx->push_back(it->vx());
    vy->push_back(it->vy());
    vz->push_back(it->vz());
    vertexChi2->push_back(it->vertexChi2());
    vertexNdof->push_back(it->vertexNdof());
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
  iEvent.put( et,  Prefix + "et" + Suffix );
  iEvent.put( mass,  Prefix + "mass" + Suffix );
  iEvent.put( mt,  Prefix + "mt" + Suffix );
  iEvent.put( globalTrack_normalizedChi2,  Prefix + "globalTracknormalizedChi2" + Suffix );
  iEvent.put( globalTrack_numberOfValidHits,  Prefix + "globalTracknumberOfValidHits" + Suffix );
  iEvent.put( globalTrack_dxy,  Prefix + "globalTrackdxy" + Suffix );
  iEvent.put( globalTrack_dxyError,  Prefix + "globalTrackdxyError" + Suffix );
  iEvent.put( isolationR03sumPt,  Prefix + "isolationR03sumPt" + Suffix );
  iEvent.put( isolationR03emEt,  Prefix + "isolationR03emEt" + Suffix );
  iEvent.put( isolationR03hadEt,  Prefix + "isolationR03hadEt" + Suffix );
  iEvent.put( vx,  Prefix + "vx" + Suffix );
  iEvent.put( vy,  Prefix + "vy" + Suffix );
  iEvent.put( vz,  Prefix + "vz" + Suffix );
  iEvent.put( vertexChi2,  Prefix + "vertexChi2" + Suffix );
  iEvent.put( vertexNdof,  Prefix + "vertexNdof" + Suffix );
//  _hi
//  iEvent.put( test,  Prefix + "test" + Suffix );
}

#endif
