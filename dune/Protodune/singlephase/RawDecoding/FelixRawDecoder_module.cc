// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// artdaq and dune-raw-data includes
#include "dune-raw-data/Overlays/FelixFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

// larsoft includes
#include "lardataobj/RawData/RawDigit.h"

// ROOT includes
#include "TH1.h"

// C++ Includes
#include <memory>
#include <iostream>

namespace dune {
  class FelixRawDecoder;
}

class dune::FelixRawDecoder : public art::EDProducer {
public:
  explicit FelixRawDecoder(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resoufelix use.

  // Plugins should not be copied or assigned.
  FelixRawDecoder(FelixRawDecoder const &) = delete;
  FelixRawDecoder(FelixRawDecoder &&) = delete;
  FelixRawDecoder & operator = (FelixRawDecoder const &) = delete;
  FelixRawDecoder & operator = (FelixRawDecoder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void reconfigure(const fhicl::ParameterSet &pset);
  void beginJob() override;

  void setRootObjects();

private:
  typedef std::vector<raw::RawDigit> RawDigits;

  bool _process(
          const artdaq::Fragment& frag, 
          RawDigits& raw_digits);

  // Declare member data here.
  std::string _input_label; 
  std::string _output_label;
  bool _expect_container_fragments;
  TH1D* _h_nframes;

  std::vector<uint16_t> _buffer;
};


dune::FelixRawDecoder::FelixRawDecoder(fhicl::ParameterSet const & pset)
  : EDProducer{pset}
// Initialize member data here.
{
  art::ServiceHandle<art::TFileService> fs;
  //fs->registerFileSwitchCallback(this, &FelixRawDecoder::setRootObjects);
  setRootObjects();

  reconfigure(pset);
  // Call appropriate produces<>() functions here.
  //produces< std::vector<raw::RawDigit> > ( _output_label );  
  produces<RawDigits>( _output_label );  
}

void dune::FelixRawDecoder::reconfigure(fhicl::ParameterSet const& pset) {

  _input_label = pset.get<std::string>("RawDataLabel");
  _output_label = pset.get<std::string>("OutputDataLabel");
  _expect_container_fragments = pset.get<bool>("ExpectContainerFragments", true);
}

void dune::FelixRawDecoder::setRootObjects(){
  art::ServiceHandle<art::TFileService> file_srv;
  _h_nframes = file_srv->make<TH1D>("felix_NFrames","FELIX: Number of frames",  1000, 0, 1000);
}
void dune::FelixRawDecoder::beginJob(){
}

void dune::FelixRawDecoder::produce(art::Event & evt){
    // TODO Use MF_LOG_DEBUG
    MF_LOG_INFO("FelixRawDecoder")
      << "-------------------- FELIX RawDecoder -------------------";

  RawDigits raw_digits;
  unsigned int n_felix_frags = 0;  

  if (_expect_container_fragments) {
    art::InputTag itag1(_input_label, "ContainerFELIX");
    auto cont_frags = evt.getHandle<artdaq::Fragments>(itag1);
    art::EventNumber_t eventNumber = evt.event();
    // Check if there is Timing data in this event
    // Don't crash code if not present, just don't save anything
    try { cont_frags->size(); }
    catch(std::exception const&) {
      std::cout << "WARNING: Container FELIX data not found in event " << eventNumber << std::endl;
      std::vector<raw::RawDigit> digits;
      evt.put(std::make_unique<std::vector<raw::RawDigit>>(std::move(digits)), _output_label);
      return;
    }
    //Check that the data is valid
    if(!cont_frags){
      MF_LOG_ERROR("FelixRawDecoder")
          << "Run: " << evt.run()
		  << ", SubRun: " << evt.subRun()
          << ", Event: " << evt.event()
		  << " Container Fragments is NOT VALID";
    }
    
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	{
	  //artdaq::Fragment frag;
          //size_t frag_size = cont_frag.fragSize(ii);
	  //frag.resizeBytes(frag_size);
	  //memcpy(frag.headerAddress(), cont_frag.at(ii), frag_size);
          //if (_process(frag, raw_digits)) ++n_felix_frags;
          if (_process(*cont_frag[ii], raw_digits)) ++n_felix_frags;
	}
    }
  }
  else
  {
    art::InputTag itag2(_input_label, "FELIX");
    auto frags = evt.getHandle<artdaq::Fragments>(itag2);
    // Check if there is Timing data in this event
    // Don't crash code if not present, just don't save anything
    art::EventNumber_t eventNumber = evt.event();
    try { frags->size(); }
    catch(std::exception const&) {
      std::cout << "WARNING: Raw FELIX data not found in event " << eventNumber << std::endl;
      std::vector<raw::RawDigit> digits;
      evt.put(std::make_unique<std::vector<raw::RawDigit>>(std::move(digits)), _output_label);
      return;
    }

    //Check that the data is valid
    if(!frags.isValid()){
      MF_LOG_ERROR("FelixRawDecoder")
          << "Run: " << evt.run()
		  << ", SubRun: " << evt.subRun()
          << ", Event: " << evt.event()
		  << " Fragments is NOT VALID";
    }

    for(auto const& frag: *frags)
    {
      if (_process(frag, raw_digits)) ++n_felix_frags;
    }
  }

  MF_LOG_INFO("FelixRawDecoder")
      << " Processed " << n_felix_frags
      << " FELIX Fragments, "
      << raw_digits.size()
      << " RawDigits.";

  evt.put(std::make_unique<decltype(raw_digits)>(std::move(raw_digits)),
          _output_label);
}

//FIXME: need to re-write this function using Felix overlay class
bool dune::FelixRawDecoder::_process(
        const artdaq::Fragment& frag, 
        RawDigits& raw_digits
        )
{
  // FIXME: Remove hard-coded fragment type
  //if((unsigned)frag.type() != 2) return false;

  MF_LOG_INFO("FelixRawDecoder")
      << "   SequenceID = " << frag.sequenceID()
      << "   fragmentID = " << frag.fragmentID()
      << "   fragmentType = " << (unsigned)frag.type()
      << "   Timestamp =  " << frag.timestamp();
  art::ServiceHandle<dune::PdspChannelMapService> channelMap;
  //Load overlay class.
  dune::FelixFragment felix(frag);
  //Get detector elemen number
  uint8_t crate = felix.crate_no(0);
  uint8_t slot = felix.slot_no(0);
  uint8_t fiber = felix.fiber_no(0); // two numbers? 
  const unsigned n_frames = felix.total_frames(); // One frame contains 20 ticks.
  std::cout<<" Nframes = "<<n_frames<<std::endl;
  _h_nframes->Fill(n_frames);
  const unsigned n_channels = dune::FelixFrame::num_ch_per_frame;// 256
  raw::RawDigit::ADCvector_t v_adc;
  //v_adc.reserve(n_frames*n_channels);
  // Fill the adc vector.
  //typedef std::tuple<uint8_t, uint8_t, uint8_t, unsigned> WireInfo_tuple; // unused
  for(unsigned ch = 0; ch < n_channels; ++ch) {
    v_adc.clear();
    std::cout<<"crate:slot:fiber = "<<crate<<", "<<slot<<", "<<fiber<<std::endl;
    std::vector<dune::adc_t> waveform( felix.get_ADCs_by_channel(ch) );
    for(unsigned int nframe=0;nframe<waveform.size();nframe++){
      if(ch==0 && nframe<100) {
        if(nframe==0) std::cout<<"Print the first 100 ADCs of Channel#1"<<std::endl;  
        std::cout<<waveform.at(nframe)<<"  ";
        if(nframe==99) std::cout<<std::endl;
      }
      v_adc.push_back(waveform.at(nframe));  
    }
    int offlineChannel = -1;
    offlineChannel = channelMap->GetOfflineNumberFromDetectorElements(crate, slot, fiber, ch,dune::PdspChannelMapService::kFELIX); // FIXME
    // Push to raw_digits.
    raw::RawDigit raw_digit(offlineChannel, n_frames, v_adc);
    raw_digits.push_back(raw_digit);
  }
  return true;
}

DEFINE_ART_MODULE(dune::FelixRawDecoder)
