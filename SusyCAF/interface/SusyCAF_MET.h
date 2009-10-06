#ifndef SUSY_CAF_MET
#define SUSY_CAF_MET

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

template< typename T >
class SusyCAF_MET : public edm::EDProducer {
 public: 
  explicit SusyCAF_MET(const edm::ParameterSet&);
 private: 
  void produce(edm::Event &, const edm::EventSetup & );

  const edm::InputTag inputTag;
  const std::string Prefix,Suffix;
};

template< typename T >
SusyCAF_MET<T>::SusyCAF_MET(const edm::ParameterSet& iConfig) :
  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  Prefix(iConfig.getParameter<std::string>("Prefix")),
  Suffix(iConfig.getParameter<std::string>("Suffix"))
{
  produces <std::vector<double> > ( Prefix + "px" + Suffix );
  produces <std::vector<double> > ( Prefix + "py" + Suffix );
}

template< typename T >
void SusyCAF_MET<T>::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<std::vector<double> >  px   ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  py   ( new std::vector<double>()  ) ;
  
  edm::Handle<std::vector<T> > metcollection;
  iEvent.getByLabel(inputTag, metcollection);

  for(typename std::vector<T>::const_iterator it = metcollection->begin(); it != metcollection->end(); ++it) {
    px->push_back(it->px());
    py->push_back(it->py());
  }
  
  iEvent.put( px,  Prefix + "px" + Suffix );
  iEvent.put( py,  Prefix + "py" + Suffix );
}

#endif
