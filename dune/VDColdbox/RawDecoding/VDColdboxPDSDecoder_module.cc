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
#include <hdf5.h>
#include "dune/DuneObj/DUNEHDF5FileInfo.h"
#include "VDColdboxHDF5Utils.h"

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

  hid_t fPrevStoredHandle = -1;
  hid_t fHDFFile = -1;
  std::string fOutputDataLabel;
  std::string fFileInfoLabel;
  bool fForceOpen;

  // Declare member data here.

};


dune::VDColdboxPDSDecoder::VDColdboxPDSDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fOutputDataLabel{p.get<std::string>("OutputDataLabel")},
    fFileInfoLabel{p.get<std::string>("FileInfoLabel", "daq")},
    fForceOpen(p.get<bool>("ForceOpen", false)) {
  produces<std::vector<raw::OpDetWaveform>>(fOutputDataLabel);
  produces<std::vector<recob::OpHit>>(fOutputDataLabel);
}

void dune::VDColdboxPDSDecoder::produce(art::Event& e) {

  using namespace dune::VDColdboxHDF5Utils;

  //To-do: put in sizes here?
  std::unique_ptr<std::vector<raw::OpDetWaveform>> output_wfs
      = std::make_unique<std::vector<raw::OpDetWaveform>>();
  std::unique_ptr<std::vector<recob::OpHit>> output_hits
      = std::make_unique<std::vector<recob::OpHit>>();


  auto infoHandle = e.getHandle<raw::DUNEHDF5FileInfo>(fFileInfoLabel);
  const std::string & group_name = infoHandle->GetEventGroupName();
  const std::string & file_name = infoHandle->GetFileName();
  hid_t file_id = infoHandle->GetHDF5FileHandle();
  
  std::cout << "HDF5 FileName: " << file_name << std::endl;
  std::cout << "Top-Level Group Name: " << group_name << std::endl;
  
  //If the fcl file said to force open the file
  //(i.e. because one is just running DataPrep), then open
  //but only if we are on a new file -- identified by if the handle
  //stored in the event is different
  if (fForceOpen && (file_id != fPrevStoredHandle)) {
    std::cout << "Opening" << std::endl;
    fHDFFile = H5Fopen(file_name.data(), H5F_ACC_RDONLY, H5P_DEFAULT);
  }//If the handle is the same, fHDFFile won't change
  else if (!fForceOpen) {
    fHDFFile = file_id;
  }
  fPrevStoredHandle = file_id;

  hid_t PDS_group = getGroupFromPath(fHDFFile, group_name + "/PDS");
  std::vector<std::string> region_names = readMidLevelGroupNames(PDS_group);
  std::cout << "Got " << region_names.size() << " regions" << std::endl;
  for (const auto & n : region_names) {
    std::cout << n << std::endl;

    hid_t region_group = getGroupFromPath(PDS_group, n);
    std::vector<std::string> element_names = readMidLevelGroupNames(region_group);
    std::cout << "Got " << element_names.size() << " elements" << std::endl;   
    for (const auto & element_name : element_names) {
      std::cout << element_name << std::endl;

      hid_t dataset = H5Dopen(region_group, element_name.data(), H5P_DEFAULT);
      hsize_t ds_size = H5Dget_storage_size(dataset);
      std::cout << "\tDataset size: " << ds_size << std::endl;

      std::vector<char> ds_data(ds_size);
      H5Dread(dataset, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              ds_data.data());
      H5Dclose(dataset);

      std::cout << "\tRead " << ds_data.size() << " bytes" << std::endl;
    }
    H5Gclose(region_group);
  }
  H5Gclose(PDS_group);

  e.put(std::move(output_wfs), fOutputDataLabel);
  e.put(std::move(output_hits), fOutputDataLabel);
}

DEFINE_ART_MODULE(dune::VDColdboxPDSDecoder)
