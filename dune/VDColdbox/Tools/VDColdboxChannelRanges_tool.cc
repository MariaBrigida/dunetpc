#include "VDColdboxChannelRanges.h"
#include "dune/ArtSupport/DuneToolManager.h"

VDColdboxChannelRanges::VDColdboxChannelRanges(fhicl::ParameterSet const& ps)
  : m_LogLevel(ps.get<int>("LogLevel", 0)) { }


//Expects a name of the form "apaX[_dccY]"
//X can be nonzero, but nothing will be returned
//
//_dccY is optional, and governs whether or not 
//ghost wires are requested, either for starting
//at 3456 or 0
IndexRange VDColdboxChannelRanges::get(std::string range_name) const {
  //The coldbox only has apa/cru0
  //If other is requested, return nothing
  if (atoi(&range_name[3]) != 1) {
    std::cout << "VDColdboxChannelRanges::get: " << range_name <<
                 " requested. Expect only APA/CRU0" << std::endl;
    return IndexRange();
  }

  //If ghost channels are requested
  if (range_name.find("dcc") != std::string::npos) {
    if (m_LogLevel > 0)
      std::cout << "VDColdboxChannelRanges::get: " <<
                   "Returning ghost channels starting at ";
    //If those channels start at 3456
    if (range_name.find("3456") != std::string::npos) {
      if (m_LogLevel > 0)
        std::cout << "3456" << std::endl;
      return IndexRange(3456, 3648);
    }
    //If they start at 0
    else {
      if (m_LogLevel > 0)
        std::cout << "0" << std::endl;
      return IndexRange(0, 192);
    }
  }

  //Return normal range for the coldbox
  if (m_LogLevel > 0)
    std::cout << "VDColdboxChannelRanges::get: " <<
                 "Returning normal range: 1600 to 3200" << std::endl;
  return IndexRange(1600, 3200);
}

DEFINE_ART_CLASS_TOOL(VDColdboxChannelRanges)
