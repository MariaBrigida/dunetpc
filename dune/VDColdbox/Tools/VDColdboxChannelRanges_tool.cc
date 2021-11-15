#include "VDColdboxChannelRanges.h"
#include "dune/ArtSupport/DuneToolManager.h"

using std::string;
using std::cout;
using std::endl;

using Name = VDColdboxChannelRanges::Name;
using Index = VDColdboxChannelRanges::Index;

VDColdboxChannelRanges::VDColdboxChannelRanges(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel", 0)),
  m_GhostRange(ps.get<IndexVector>("GhostRange")),
  m_glo(0), m_ghi(0) {
  const string myname = "VDColdboxChannelRanges::ctor: ";
  const Index nu = 384;
  const Index ny = 640;
  const Index ntot = 3200;
  const Index nhaf = ntot/2;
  if ( m_GhostRange.size() == 2 && m_GhostRange[1] >= m_GhostRange[0] ) {
    m_glo = m_GhostRange[0];
    m_ghi = m_GhostRange[1] + 1;
  } else if ( m_GhostRange.size() ) {
    cout << myname << "WARNING: " << "Ignoring invalid ghost range." << endl;
  }
  if ( m_LogLevel >= 1 ) {
    cout << myname << "     LogLevel: " << m_LogLevel << endl;
    cout << myname << "  Ghost range: [";
    if ( m_ghi > m_glo ) cout << m_glo << ", " << m_ghi-1;
    cout << "]" << endl;
  }
  // Build range map.
  insert("cru",          0,       ntot, "CRU");
  insert("crt",          0,       nhaf, "CRT");
  insert("crb",       1600,       ntot, "CRB");
  insert("crtu",         0,         nu, "CRTu");
  insert("crty",        nu,      nu+ny, "CRTy");
  insert("crtz",     nu+ny,       nhaf, "CRTz");
}


IndexRange VDColdboxChannelRanges::get(string sran) const {
  const string myname = "VDColdboxChannelRanges::get: ";
  const Index nu = 384;
  const Index ny = 640;
  const Index ntot = 3200;
  const Index nhaf = ntot/2;
  RangeMap::const_iterator iran = m_rans.find(sran);
  if ( iran  != m_rans.end() ) return iran->second;
  //if ( sran == "cru" )  return IndexRange(sran,          0,       ntot, "CRU");
  if ( sran == "crt" )  return IndexRange(sran,          0,       nhaf, "CRT");
  if ( sran == "crb" )  return IndexRange(sran,       1600,       ntot, "CRB");
  if ( sran == "crtu" ) return IndexRange(sran,          0,         nu, "CRTu");
  if ( sran == "crty" ) return IndexRange(sran,         nu,      nu+ny, "CRTy");
  if ( sran == "crtz" ) return IndexRange(sran,      nu+ny,       nhaf, "CRTz");
  if ( sran == "crbu" ) return IndexRange(sran,       nhaf,    nhaf+nu, "CRBu");
  if ( sran == "crby" ) return IndexRange(sran,    nhaf+nu, nhaf+nu+ny, "CRBy");
  if ( sran == "crbz" ) return IndexRange(sran, nhaf+nu+ny,       ntot, "CRBz");
  if ( sran == "crbg" ) return IndexRange(sran,      m_glo,      m_ghi, "CRBghost");
  if ( m_LogLevel >= 2 ) {
    cout << myname << "Invalid channel range name: " << sran << endl;
  }
  return IndexRange(0, 0);
}

void VDColdboxChannelRanges::insert(Name sran, Index ich1, Index ich2, Name slab1) {
  m_rans[sran] = IndexRange(sran, ich1, ich2, slab1);
}

DEFINE_ART_CLASS_TOOL(VDColdboxChannelRanges)
