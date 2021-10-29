#ifndef VDColdboxChannelRanges_H
#define VDColdboxChannelRanges_H

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/DuneInterface/Tool/IndexRangeTool.h"
#include <map>

class VDColdboxChannelRanges : public IndexRangeTool {

public:

  /*
  using Name = std::string;
  using NameVector = std::vector<Name>;
  using Index = IndexRange::Index;
  using IndexVector = std::vector<Index>;
  using IndexRangeMap = std::map<Name, IndexRange>;
  */

  // Ctor.
  VDColdboxChannelRanges(fhicl::ParameterSet const& ps);

  // Dtor.
  ~VDColdboxChannelRanges() override =default;

  // Return a range.
  IndexRange get(std::string range_name) const override;

private:

  // Configuration parameters.
  int m_LogLevel;
  //IndexVector m_ApaNumbers;
  //NameVector m_ApaLocationNames;
  //Name m_ExtraRanges;

  //IndexRangeMap m_Ranges;
  //const IndexRangeTool* m_pExtraRanges =nullptr;

  // Add an entry to the range map.
  //void insertLen(Name nam, Index begin, Index len, Name lab, Name lab1 ="", Name lab2 ="");

};
#endif
