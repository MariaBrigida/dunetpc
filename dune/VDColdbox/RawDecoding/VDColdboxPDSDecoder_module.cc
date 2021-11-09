////////////////////////////////////////////////////////////////////////
// Class:       VDColdboxPDSDecoder
// Plugin Type: producer (Unknown Unknown)
// File:        VDColdboxPDSDecoder_module.cc
//
// Generated at Tue Nov  9 16:28:27 2021 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <memory>

namespace dune {
  class VDColdboxPDSDecoder;
}


class dune::VDColdboxPDSDecoder : public art::EDProducer {
public:
  explicit VDColdboxPDSDecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VDColdboxPDSDecoder(VDColdboxPDSDecoder const&) = delete;
  VDColdboxPDSDecoder(VDColdboxPDSDecoder&&) = delete;
  VDColdboxPDSDecoder& operator=(VDColdboxPDSDecoder const&) = delete;
  VDColdboxPDSDecoder& operator=(VDColdboxPDSDecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string fOutputDataLabel;

  // Declare member data here.

};


dune::VDColdboxPDSDecoder::VDColdboxPDSDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fOutputDataLabel{p.get<std::string>("OutputDataLabel")} {
  produces<std::vector<raw::OpDetWaveform>>(fOutputDataLabel);
  produces<std::vector<recob::OpHit>>(fOutputDataLabel);
}

void dune::VDColdboxPDSDecoder::produce(art::Event& e) {
  //To-do: put in sizes here?
  std::unique_ptr<std::vector<raw::OpDetWaveform>> output_wfs
      = std::make_unique<std::vector<raw::OpDetWaveform>>();
  std::unique_ptr<std::vector<recob::OpHit>> output_hits
      = std::make_unique<std::vector<recob::OpHit>>();

  e.put(std::move(output_wfs), fOutputDataLabel);
  e.put(std::move(output_hits), fOutputDataLabel);
}

DEFINE_ART_MODULE(dune::VDColdboxPDSDecoder)
