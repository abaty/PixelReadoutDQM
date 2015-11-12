//#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

//
// class declaration
//

class HLTPixelActivityHFSumEnergyFilter : public edm::stream::EDFilter<> {
public:
  explicit HLTPixelActivityHFSumEnergyFilter(const edm::ParameterSet&);
  ~HLTPixelActivityHFSumEnergyFilter();
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override; //, trigger::TriggerFilterObjectWithRefs & filterproduct) const override;
 virtual void endStream() override;
  
  edm::InputTag inputTag_;          // input tag identifying product containing pixel digis
  //edm::EDGetTokenT<edm::DetSetVector<PixelDigi> > inputToken_;
  edm::EDGetTokenT<edm::DetSetVector<SiPixelCluster> > inputToken_;
  edm::EDGetTokenT<HFRecHitCollection> HFHitsToken_;
  edm::InputTag HFHits_;
  double eCut_HF_;
  double eMin_HF_;
  double offset_;
  double slope_;

};

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

//
// constructors and destructor
//

HLTPixelActivityHFSumEnergyFilter::HLTPixelActivityHFSumEnergyFilter(const edm::ParameterSet& config)
{
  inputTag_      = config.getParameter<edm::InputTag>("inputTag");
  HFHits_        = config.getParameter<edm::InputTag>("HFHitCollection");
  eCut_HF_       = config.getParameter<double>("eCut_HF");
  eMin_HF_       = config.getParameter<double>("eMin_HF");
  offset_        = config.getParameter<double>("offset");
  slope_         = config.getParameter<double>("slope");

  inputToken_ = consumes<edm::DetSetVector<SiPixelCluster> >(inputTag_);
  //inputToken_ = consumes<edm::DetSetVector<PixelDigi> >(inputTag_);
  HFHitsToken_ = consumes<HFRecHitCollection>(HFHits_); 
  
  LogDebug("") << "Using the " << inputTag_ << " input collection";
  //  LogDebug("") << "Requesting at least " << min_clusters_ << " clusters";
  // if(max_clusters_ > 0)
  // LogDebug("") << "...but no more than " << max_clusters_ << " clusters";
}

HLTPixelActivityHFSumEnergyFilter::~HLTPixelActivityHFSumEnergyFilter()
{
}

void
HLTPixelActivityHFSumEnergyFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //makeEDFilterDescription(desc);
  desc.add<edm::InputTag>("inputTag",edm::InputTag("hltSiPixelClusters"));
  desc.add<edm::InputTag>("inputHFTag",edm::InputTag("hfHitCollection"));
  desc.add<double>("eCut_HF",0);
  desc.add<double>("eMin_HF",0);
  desc.add<double>("offset",0);
  desc.add<double>("slope",0);
  descriptions.add("hltPixelActivityHFSumEnergyFilter",desc);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool HLTPixelActivityHFSumEnergyFilter::filter(edm::Event& event, const edm::EventSetup& iSetup) { //, trigger::TriggerFilterObjectWithRefs & filterproduct) const
  //{

  Printf("running filter");
  
  using namespace edm;

  // The filter object
  //  if (saveTags()) filterproduct.addCollectionTag(inputTag_);

  // get hold of products from Event
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > clusterColl;
  event.getByToken(inputToken_, clusterColl);

  unsigned int clusterSize = clusterColl->dataSize();
  LogDebug("") << "Number of clusters accepted: " << clusterSize;

  edm::Handle<HFRecHitCollection> HFRecHitsH;
  event.getByToken(HFHitsToken_,HFRecHitsH);

  double sumE = 0.;

  for (HFRecHitCollection::const_iterator it=HFRecHitsH->begin(); it!=HFRecHitsH->end(); it++) {
    if (it->energy()>eCut_HF_) {
      sumE += it->energy();
    }
  }

  bool accept = kFALSE;

  double thres = offset_ + slope_ * clusterSize;
  if(sumE>eMin_HF_ && sumE<thres) accept = kTRUE;

  Printf("clusterSize: %d  sumE: %f thres: %f ",clusterSize,sumE,thres);
  
  // return with final filter decision
  return accept;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HLTPixelActivityHFSumEnergyFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HLTPixelActivityHFSumEnergyFilter::endStream() {
}

// define as a framework module

DEFINE_FWK_MODULE(HLTPixelActivityHFSumEnergyFilter);
