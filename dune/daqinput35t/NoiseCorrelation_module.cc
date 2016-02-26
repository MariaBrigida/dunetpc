////////////////////////////////////////////////////////////////////////
// Class:       NoiseCorrelation
// Module Type: analyser
// File:        NoiseCorrelation_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), Feb 2016
// 
// Looks at the correlation between waveforms on all combinations
// of channels.
// Runs over an artdaq-formatted file and produces a tree or histogram,
// filled for each channel combination and the respective correlation.
////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// lbne-artdaq
#include "lbne-raw-data/Overlays/TpcMilliSliceFragment.hh"
#include "artdaq-core/Data/Fragments.hh"

// larsoft
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "tpcFragmentToRawDigits.h"
#include "utilities/UnpackFragment.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "dune/RunHistory/DetPedestalDUNE.h"
#include "cetlib/getenv.h"

// c++
#include <memory>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>

// ROOT
#include "TMath.h"
#include "TTree.h"
#include "TH2F.h"
#include "TStyle.h"

namespace DAQToOffline {
  class NoiseCorrelation;
}

class DAQToOffline::NoiseCorrelation : public art::EDAnalyzer {
public:

  explicit NoiseCorrelation(fhicl::ParameterSet const & pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NoiseCorrelation(NoiseCorrelation const &) = delete;
  NoiseCorrelation(NoiseCorrelation &&) = delete;
  NoiseCorrelation & operator = (NoiseCorrelation const &) = delete;
  NoiseCorrelation & operator = (NoiseCorrelation &&) = delete;

  void analyze(art::Event const& evt) override;
  void reconfigure(const fhicl::ParameterSet &pset);

private:

  std::string fFragType;
  std::string fRawDataLabel;

  // Correlation tree
  TTree* fCorrelationTree;
  TH2F* fCorrelationHist;
  float fCorrelation;
  int fChannel1;
  int fChannel2;

  std::map<int,int> fChannelMap;

  const lariov::DetPedestalProvider& fPedestalRetrievalAlg = *(lar::providerFrom<lariov::DetPedestalService>());

};

DAQToOffline::NoiseCorrelation::NoiseCorrelation(fhicl::ParameterSet const & pset) : art::EDAnalyzer(pset) {

  this->reconfigure(pset);
  gStyle->SetOptStat(0);

  art::ServiceHandle<art::TFileService> tfs;
  fCorrelationTree = tfs->make<TTree>("Correlations","Correlations");
  fCorrelationHist = tfs->make<TH2F>("Correlation","Correlation",2048,0,2048,2048,0,2048);
  fCorrelationTree->Branch("Channel1",&fChannel1);
  fCorrelationTree->Branch("Channel2",&fChannel2);
  fCorrelationTree->Branch("Correlation",&fCorrelation);

  BuildTPCChannelMap("rce_channel_map_dune35t.txt",fChannelMap);

}

void DAQToOffline::NoiseCorrelation::reconfigure(fhicl::ParameterSet const& pset) {
  fFragType = pset.get<std::string>("FragType");
  fRawDataLabel = pset.get<std::string>("RawDataLabel");
}

void DAQToOffline::NoiseCorrelation::analyze(art::Event const& evt) {

  art::Handle<artdaq::Fragments> rawFragments;
  evt.getByLabel(fRawDataLabel, fFragType, rawFragments);

  art::EventNumber_t eventNumber = evt.event();

  // Check that the data are valid
  if(!rawFragments.isValid()){
    std::cerr << "Run: " << evt.run()
	      << ", SubRun: " << evt.subRun()
	      << ", Event: " << eventNumber
	      << " is NOT VALID" << std::endl;
    throw cet::exception("rawFragments NOT VALID");
  }

  // Create a map containing (fragmentID, fragIndex) for the event, will be used to check if each channel is present
  std::map < unsigned int, unsigned int > mapFragID;
  for(size_t fragIndex = 0; fragIndex < rawFragments->size(); fragIndex++){
    const artdaq::Fragment &singleFragment = (*rawFragments)[fragIndex];
    unsigned int fragmentID = singleFragment.fragmentID();
    mapFragID.insert(std::pair<unsigned int, unsigned int>(fragmentID,fragIndex));
  }

  // Create a 2D vector to save the correlations
  std::vector<std::vector<float> > correlationArray(2048,std::vector<float>(2048,0));

  art::ServiceHandle<geo::Geometry> geometry;
  size_t numChannels = geometry->Nchannels();

  // Get the ADC vector for each channel
  std::vector<std::vector<short> > adcVectors;

  for (size_t channel = 0; channel < numChannels; ++channel) {

    // Follow the steps written by J Davies in TpcDAQToOffline
    unsigned int fragmentID = UnpackFragment::getFragIDForChan(channel);
    unsigned int sample = UnpackFragment::getNanoSliceSampleForChan(channel);
    std::vector<short> adcvec;
    if (mapFragID.find(fragmentID) != mapFragID.end()) {
      unsigned int fragIndex = mapFragID[fragmentID];
      const artdaq::Fragment &singleFragment = (*rawFragments)[fragIndex];
      lbne::TpcMilliSliceFragment millisliceFragment(singleFragment);
      auto numMicroSlices = millisliceFragment.microSliceCount();
      for(unsigned int i_micro = 0; i_micro < numMicroSlices; i_micro++) {
	std::unique_ptr <const lbne::TpcMicroSlice> microSlice = millisliceFragment.microSlice(i_micro);
	auto numNanoSlices = microSlice->nanoSliceCount();
	for(uint32_t i_nano = 0; i_nano < numNanoSlices; i_nano++){
	  uint16_t val = std::numeric_limits<uint16_t>::max();
	  bool success = microSlice->nanosliceSampleValue(i_nano, sample, val);
	  if (success) {
	    float const pedestal = fPedestalRetrievalAlg.PedMean(fChannelMap.at(channel));
	    short adc = short(val - pedestal);
	    if ((adc & 0x3F) == 0x0 || (adc & 0x3F) == 0x3F)
	      adcvec.push_back(short(-999));
	    else
	      adcvec.push_back(adc);
	  }
	}
      }
    }

    adcVectors.push_back(adcvec);

  }

  // Look at all possible channel combinations
  for (size_t chan1 = 0; chan1 < numChannels; ++chan1) {
    for (size_t chan2 = chan1; chan2 < numChannels; ++chan2) {

      fChannel1 = chan1;
      fChannel2 = chan2;
      fCorrelation = -999;

      std::vector<short> adcvec1 = adcVectors[chan1];
      std::vector<short> adcvec2 = adcVectors[chan2];

      if (fChannel1 % 50 == 0 and fChannel2 % 1000 == 0)
	std::cout << "Looking at correlations between channel " << fChannel1 << " and " << fChannel2 << std::endl;

      // Now we have the ADC vectors, we can look at the correlations
      // Use the Pearson Correlation Coefficient

      double n1 = adcvec1.size(), n2 = adcvec2.size();
      if (n1 != n2)
	continue;

      double sumxy = 0, sumx = 0, sumy = 0, sumx2 = 0, sumy2 = 0;
      for (size_t tick = 0; tick < n1; ++tick) {
	short adc1 = adcvec1[tick], adc2 = adcvec2[tick];
	if (adc1 == -999 or adc2 == -999)
	  continue;
	sumx += adc1;
	sumy += adc2;
	sumxy += adc1*adc2;
	sumx2 += adc1*adc1;
	sumy2 += adc2*adc2;
      }

      double denom = ( (n1*sumx2) - (sumx*sumx) ) * ( (n1*sumy2) - (sumy*sumy) );
      if (n1 > 0 and denom > 0) {
	fCorrelation = ( (n1*sumxy) - (sumx*sumy) ) / ( TMath::Sqrt( denom ) );
	correlationArray[fChannel1][fChannel2] = fCorrelation;
	correlationArray[fChannel2][fChannel1] = fCorrelation;
      }
      //fCorrelationTree->Fill();

    }
  } // channel loops

  for (unsigned int channel1 = 0; channel1 < correlationArray.size(); ++channel1)
    for (unsigned int channel2 = 0; channel2 < correlationArray.size(); ++channel2)
	fCorrelationHist->SetBinContent(channel1, channel2, correlationArray[channel1][channel2]);
  fCorrelationHist->GetZaxis()->SetRangeUser(-1,1);

  return;

}

DEFINE_ART_MODULE(DAQToOffline::NoiseCorrelation)
