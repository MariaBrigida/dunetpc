#include "VDColdboxChannelGroups.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"

VDColdboxChannelGroups::VDColdboxChannelGroups(fhicl::ParameterSet const& ps)
  : m_LogLevel(ps.get<int>("LogLevel", 0)) { }

IndexRangeGroup VDColdboxChannelGroups::get(std::string group_name) const {
  std::cout << "VDColdboxChannelGroups::get: " <<
               "attempting to return range group for " << group_name <<
               std::endl;
  std::vector<IndexRange> ranges;

  if (atoi(&group_name[3]) == 1) {
    std::cout << "VDColdboxChannelGroups::get: " <<
                 "Apa 0" << std::endl;
    ranges.push_back(IndexRange("apa1", 1600, 3200));
  }

  if (group_name.find("dcc") != std::string::npos) {
    if (group_name.find("3456") != std::string::npos) {
      std::cout << "VDColdboxChannelGroups::get: " <<
                   "Ghost channels 3456" << std::endl;
      ranges.push_back(IndexRange("apa1_dcc3456", 3456, 3648));
    }
    else {
      std::cout << "VDColdboxChannelGroups::get: " <<
                   "Ghost channels 3456" << std::endl;
      ranges.push_back(IndexRange("apa1_dcc0", 0, 192));
    }
  }
  
  std::cout << "VDColdboxChannelGroups::get: " <<
               "Group size: " << ranges.size() << std::endl;
  return IndexRangeGroup(group_name, ranges);
}

DEFINE_ART_CLASS_TOOL(VDColdboxChannelGroups)
