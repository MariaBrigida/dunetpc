#ifndef PROTODUNETIMESTAMP_H
#define PROTODUNETIMESTAMP_H

#include "Rtypes.h"
// I want to do this to get the TimingCommand enum from dune-raw-data:
//
// #include "dune-raw-data/Overlays/TimingFragment.hh" // for TimingCommand enum
// 
// but TimingFragment #includes Fragment.hh from artdaq-core, and
// large parts of artdaq::Fragment appear to be hidden from ROOT with
// #ifdefs, so I get compile errors from rootcling, eg:
// 
// TimingFragment.hh:115:77: error: no member named 'dataBeginBytes' in 'artdaq::Fragment'
//
// Life's too short, so I'm copying the enum here with a different name

namespace dune {

    enum class ProtoDUNETimingCommand {
        // From https://twiki.cern.ch/twiki/bin/view/CENF/TimingSystemAdvancedOp retrieved on 2018-09-07
        // The 'sync' bus has the following commands at the moment:
        TimeSync   = 0x0,
        Echo       = 0x1,
        SpillStart = 0x2,
        SpillStop  = 0x3,
        RunStart   = 0x4,
        RunStop    = 0x5,
        WibCalib   = 0x6,
        SSPCalib   = 0x7,
        FakeTrig0  = 0x8,
        FakeTrig1  = 0x9,
        FakeTrig2  = 0xa,
        FakeTrig3  = 0xb,
        BeamTrig   = 0xc,
        NoBeamTrig = 0xd,
        ExtFakeTrig= 0xe
  };

    class ProtoDUNETimeStamp
    {
    public:
        ProtoDUNETimeStamp();
        
        UInt_t getCookie()  const { return m_cookie; } 
        ProtoDUNETimingCommand getTriggerType()  const { return m_triggerType; } 
        UInt_t getReservedBits()  const { return m_reservedBits; } 
        ULong64_t getTimeStamp()  const { return m_timeStamp; } 
        UInt_t getEventCounter()  const { return m_eventCounter; } 
        bool isChecksumGood()  const { return m_checksumGood; }
        ULong64_t getLastRunStart()  const { return m_lastRunStart; } 
        ULong64_t getLastSpillStart()  const { return m_lastSpillStart; } 
        ULong64_t getLastSpillEnd()  const { return m_lastSpillEnd; } 
        UInt_t getVersion() const { return m_version; }

        void setCookie(UInt_t arg) { m_cookie = arg; } 
        void setTriggerType(dune::ProtoDUNETimingCommand arg) { m_triggerType = arg; } 
        void setReservedBits(UInt_t arg) { m_reservedBits = arg; } 
        void setTimeStamp(ULong64_t arg) { m_timeStamp = arg; } 
        void setEventCounter(UInt_t arg) { m_eventCounter = arg; } 
        void setChecksumGood(bool arg) {m_checksumGood = arg;}
        void setLastRunStart(ULong64_t arg) { m_lastRunStart = arg; } 
        void setLastSpillStart(ULong64_t arg) { m_lastSpillStart = arg; } 
        void setLastSpillEnd(ULong64_t arg) { m_lastSpillEnd = arg; } 
        void setVersion(UInt_t arg) { m_version = arg; }

    protected:
        /// The cookie which identifies the event as coming from the timing board
        UInt_t m_cookie; 
        /// The type of trigger command
        ProtoDUNETimingCommand m_triggerType; 
        /// Reserved bits from the trigger command word (ought to be zero)
        UInt_t m_reservedBits; 
        /// The 50 MHz timestamp of the trigger
        ULong64_t m_timeStamp; 
        /// The event counter from the timing board reader
        UInt_t m_eventCounter; 
        /// Whether the checksum is good (always true)
        bool m_checksumGood;
        /// The timestamp of the last run start seen by the timing board reader
        ULong64_t m_lastRunStart; 
        /// The timestamp of the last spill start seen by the timing board reader
        ULong64_t m_lastSpillStart; 
        /// The timestamp of the last spill end seen by the timing board reader
        ULong64_t m_lastSpillEnd;
        /// The version of the artdaq timing fragment that this object is made from
        UInt_t m_version;
    };

}

#endif
