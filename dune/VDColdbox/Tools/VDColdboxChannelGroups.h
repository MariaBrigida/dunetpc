#ifndef VDColdboxChannelGroups_H
#define VDColdboxChannelGroups_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeGroupTool.h"
#include <map>

class IndexRangeTool;

class VDColdboxChannelGroups : public IndexRangeGroupTool {

public:

  /*
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = IndexRangeGroup::Index;
  using IndexVector = std::vector<Index>;
  using GroupMap = std::map<Name, NameVector>;
  */

  // Ctor.
  VDColdboxChannelGroups(fhicl::ParameterSet const& ps);

  // Dtor.
  ~VDColdboxChannelGroups() override =default;

  // Return a range.
  IndexRangeGroup get(std::string nam) const override;

private:

  // Configuration parameters.
  int m_LogLevel;
  /*
  IndexVector m_ApaNumbers;
  Name m_IndexRangeTool;

  // Derived from configuration.
  const IndexRangeTool* m_pIndexRangeTool =nullptr;
  GroupMap m_groups;
  GroupMap m_labels;*/

};


#endif
