// test_StandardRawDigitPrepService.cxx
//
// David Adams
// May 2016
//
// Test StandardRawDigitPrepService.

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/RawDigitPrepService.h"
#include "dune/DuneInterface/WiredAdcChannelDataMap.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::fixed;
using std::vector;
using std::map;
using art::ServiceHandle;
using raw::RawDigit;

//**********************************************************************

bool sigequal(AdcSignal sig1, AdcSignal sig2) {
  AdcSignal sigdiff = sig2 - sig1;
  if ( sigdiff < -0.5 || sigdiff >  0.5 ) {
    cout << "sigequal: " << sig1 << " != " << sig2 << endl;
    return false;
  }
  return true;
}

bool flagequal(AdcFlag flg1, AdcFlag flg2) {
  if ( flg1 == flg2 ) return true;
  if ( flg1 == AdcInterpolated && flg2 == AdcStuckOff ) return true;
  if ( flg1 == AdcInterpolated && flg2 == AdcStuckOn ) return true;
  if ( flg1 == AdcSetFixed ) return true;
  cout << "flagequal: " << flg1 << " != " << flg2 << endl;
  return false;
}

//**********************************************************************

// If usePedestalAdjustment is true, then extra pedestals are added and
// then removed using the DoPedestalAdjustment option.

int test_StandardRawDigitPrepService(bool useExistingFcl =false, bool useFclFile =false) {
  const string myname = "test_StandardRawDigitPrepService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  string fclfile = "test_StandardRawDigitPrepService.fcl";
  bool usePedestalAdjustment = false;
  bool checkExtracted = false;
  bool checkMitigated = false;
  bool checkNoiseRemoved = false;
  if ( useExistingFcl ) {
  } else if ( useFclFile ) {
    ofstream fout(fclfile.c_str());
    fout << "#include \"services_dune.fcl\"" << endl;
    fout << "services: @local::dune35tdata_reco_services" << endl;
    fout << "services.RawDigitPrepService.DoMitigation: false" << endl;
    fout << "services.RawDigitPrepService.DoNoiseRemoval: false" << endl;
    fout << "services.RawDigitPrepService.DoDeconvolution: false" << endl;
    fout.close();
  } else {
    ofstream fout(fclfile.c_str());
    fout << "#include \"services_dune.fcl\"" << endl;
    fout << "services.user: @local::dune35t_services" << endl;
    fout << "services.RawDigitExtractService: {" << endl;
    fout << "  service_provider: StandardRawDigitExtractService" << endl;
    fout << "  LogLevel:        1" << endl;
    fout << "  PedestalOption:  1" << endl;
    fout << "  FlagStuckOff: true" << endl;
    fout << "  FlagStuckOn:  true" << endl;
    fout << "}" << endl;
    fout << "services.AdcMitigationService: {" << endl;
    fout << "  service_provider: InterpolatingAdcMitigationService" << endl;
    fout << "  LogLevel:              1" << endl;
    fout << "  SkipUnderflows:     true" << endl;
    fout << "  SkipOverflows:      true" << endl;
    fout << "  MaxConsecutiveSamples: 3" << endl;
    fout << "  MaxConsecutiveFlag:    0" << endl;
    fout << "}" << endl;
    fout << "services.PedestalEvaluationService: {" << endl;
    fout << "  service_provider: MedianPedestalService" << endl;
    fout << "  LogLevel:           1" << endl;
    fout << "  SkipFlaggedSamples: true" << endl;
    fout << "  SkipSignals:        true" << endl;
    fout << "}" << endl;
    fout << "services.AdcChannelNoiseRemovalService: {" << endl;
    fout << "  service_provider: ThresholdNoiseRemovalService" << endl;
    fout << "  Threshold:  10.0" <<  endl;
    fout << "  LogLevel:      1" << endl;
    fout << "}" << endl;
    fout << "services.AdcNoiseRemovalService: {" << endl;
    fout << "  service_provider: MultiChannelNoiseRemovalService" << endl;
    fout << "  LogLevel:       1" << endl;
    fout << "}" << endl;
    fout << "services.AdcChannelDataCopyService: {" << endl;
    fout << "  service_provider: ConfigurableAdcChannelDataCopyService" << endl;
    fout << "  LogLevel: 1" << endl;
    fout << "  CopyChannel:     true" << endl;
    fout << "  CopyPedestal:    true" << endl;
    fout << "  CopySamples:     true" << endl;
    fout << "  CopyFlags:       true" << endl;
    fout << "}" << endl;
    fout << "services.AdcRoiBuildingService: {" << endl;
    fout << "  service_provider: KeepAllRoiBuildingService" << endl;
    fout << "  LogLevel:       1" << endl;
    fout << "}" << endl;
    fout << "services.AdcWireBuildingService: {" << endl;
    fout << "  service_provider: StandardAdcWireBuildingService" << endl;
    fout << "  LogLevel:       1" << endl;
    fout << "}" << endl;
    fout << "services.RawDigitPrepService: {" << endl;
    fout << "  service_provider: StandardRawDigitPrepService" << endl;
    fout << "  LogLevel:                 1" << endl;
    fout << "  SkipBad:              false" << endl;
    fout << "  SkipNoisy:            false" << endl;
    fout << "  DoMitigation:          true" << endl;
    fout << "  DoEarlySignalFinding: false" << endl;
    fout << "  DoNoiseRemoval:        true" << endl;
    fout << "  DoDeconvolution:      false" << endl;
    fout << "  DoROI:                 true" << endl;
    fout << "  DoWires:               true" << endl;
    fout << "  DoPedestalAdjustment: false" << endl;
    fout << "  IntermediateStates:      [\"extracted\", \"mitigated\", \"noiseRemoved\"]" << endl;
    fout << "}" << endl;
    fout.close();
    usePedestalAdjustment = false;
    checkExtracted = true;
    checkMitigated = true;
    checkNoiseRemoved = true;
  }

  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Add services." << endl;
  assert( ash.addServices(fclfile, true) == 0 );
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  assert( ash.loadServices() == 1 );
  ash.print();

  AdcCount lowbits = 63;
  AdcCount highbits = 4095 - 63;

  cout << myname << line << endl;
  cout << myname << "Create raw digits." << endl;
  AdcChannel nchan = 8;
  unsigned int nsig = 64;
  AdcSignalVectorVector sigsin(nchan);
  float fac = 250.0;
  unsigned int isig_stucklo = 15;
  unsigned int isig_stuckhi = 25;
  AdcSignal peds[nchan] = {2000.2, 2010.1, 2020.3, 1990.4, 1979.6, 1979.2, 1995.0, 2001.3};
  AdcSignal xpeds[nchan] = {0, 0, 0, 0, 0, 0, 0, 0};  // Need pedestal adju
  if ( usePedestalAdjustment ) {
    xpeds[4] = 100.0;
    xpeds[5] = 100.0;
    xpeds[6] = 100.0;
    xpeds[7] = 100.0;
  }
  vector<RawDigit> digs;
  map<AdcChannel, AdcCountVector> adcsmap;
  map<AdcChannel, AdcFlagVector> expflagsmap;
  // Each channel is filled with a bipolar signal offset by one from preceding channel.
  for ( AdcChannel chan=0; chan<nchan; ++chan ) {
    unsigned int isig1 = 10 + chan;
    for ( unsigned int isig=0; isig<isig1; ++isig ) sigsin[chan].push_back(0);
    for ( unsigned int i=0; i<10; ++i )             sigsin[chan].push_back(fac*i);
    for ( unsigned int i=10; i<1000; --i )            sigsin[chan].push_back(fac*i);
    for ( unsigned int i=19; i<1000; --i )            sigsin[chan].push_back(-sigsin[chan][i+isig1]);
    for ( unsigned int isig=sigsin[chan].size();
                               isig<nsig; ++isig )  sigsin[chan].push_back(0);
    assert(sigsin[chan].size() == nsig);
    AdcCountVector adcsin;
    for ( unsigned int isig=0; isig<nsig; ++isig) {
      AdcSignal sig = sigsin[chan][isig] + peds[chan] + xpeds[chan];
      AdcCount adc = 0;
      if ( sig > 0.0 ) adc = int(sig+0.5);
      if ( adc > 4095 ) adc = 4095;
      AdcCount adchigh = adc & highbits;
      if ( isig == isig_stucklo ) adc = adchigh;           // Stuck low bits to zero.
      if ( isig == isig_stuckhi ) adc = adchigh + lowbits; // Stuck low bits to one.
      adcsin.push_back(adc);
    }
    assert(adcsin.size() == nsig);
    adcsmap[chan] = adcsin;
    // Create RawDigit.
    RawDigit& dig = *digs.emplace(digs.end(), chan, nsig, adcsin, raw::kNone);
    dig.SetPedestal(peds[chan]);
    cout << myname << "    Compressed size: " << dig.NADC() << endl;
    cout << myname << "  Uncompressed size: " << dig.Samples() << endl;
    cout << myname << "           Pedestal: " << dig.GetPedestal() << endl;
    cout << myname << "            Channel: " << dig.Channel() << endl;
    assert(dig.Samples() == nsig);
    assert(dig.Channel() == chan);
    assert(dig.GetPedestal() == peds[chan]);
  }
  assert( adcsmap.size() == nchan );

  cout << myname << line << endl;
  cout << myname << "Create the expected flag vector." << endl;
  for ( AdcChannel chan=0; chan<nchan; ++chan ) {
    AdcFlagVector expflags(nsig, AdcGood);
    for ( unsigned int isig=0; isig<nsig; ++isig) {
      AdcCount adc = adcsmap[chan][isig];
      AdcCount adclow = adc & lowbits;
      if ( adc <= 0 ) expflags[isig] = AdcUnderflow;
      else if ( adc >= 4095 ) expflags[isig] = AdcOverflow;
      else if ( adclow == 0 ) expflags[isig] = AdcStuckOff;
      else if ( adclow == lowbits ) expflags[isig] = AdcStuckOn;
    }
    expflagsmap[chan] = expflags;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch raw digit prep service." << endl;
  ServiceHandle<RawDigitPrepService> hrdp;
  hrdp->print();
  ash.print();

  cout << myname << line << endl;
  cout << myname << "Prep data from digits." << endl;
  AdcSignalVector sigs;
  AdcFlagVector flags;
  AdcChannelDataMap prepdigs;
  std::vector<recob::Wire> wires;
  wires.reserve(nchan);
  WiredAdcChannelDataMap intStates;
  assert( hrdp->prepare(digs, prepdigs, &wires, &intStates) == 0 );
  cout << myname << "      # prepared digits: " << prepdigs.size() << endl;
  cout << myname << "                # wires: " << wires.size() << endl;
  cout << myname << "  # intermediate states: " << intStates.dataMaps.size() << endl;
  for ( const auto& namedadm : intStates.dataMaps ) {
    string sname = namedadm.first;
    const AdcChannelDataMap& adm = namedadm.second;
    auto iwco = intStates.wires.find(sname);
    if ( iwco == intStates.wires.end() ) {
      cout << myname << "Wires not found for intermediate state " << sname << "." << endl;
      assert(false);
    }
    const vector<recob::Wire>* pwires = iwco->second;
    assert( pwires != nullptr );
    cout << myname << "  State " << sname << " has " << adm.size() << " ADC channels";
    cout << " and " << pwires->size() << " wires";
    cout <<"." << endl;
  }
  cout << myname << "   # intermediate wires: " << intStates.wires.size() << endl;
  for ( AdcChannelDataMap::const_iterator ichdat=prepdigs.begin(); ichdat!=prepdigs.end(); ++ichdat ) {
    AdcChannel chan = ichdat->first;
    const AdcChannelData& acd = ichdat->second;
    const AdcSignalVector& sigs = acd.samples;
    const AdcFlagVector& flags = acd.flags;
    const raw::RawDigit* pdig = acd.digit;
    AdcSignal ped = acd.pedestal;
    cout << myname << "----- Channel " << chan << endl;
    cout << myname << "Output vector size: " << sigs.size() << endl;
    cout << myname << " Output flags size: " << flags.size() << endl;
    cout << myname << "          Pedestal: " << ped << endl;
    cout << myname << "        samples[0]: " << sigs[0] << endl;
    cout << myname << "Check final data." << endl;
    assert( sigs.size() == nsig );
    assert( flags.size() == nsig );
    assert( pdig != nullptr );
    assert( pdig == &digs[chan] );
    assert( pdig->Channel() == chan );
    assert( ichdat->second.digitIndex == chan );
    const AdcFlagVector& expflags = expflagsmap[chan];
    assert( expflagsmap[chan].size() == nsig );
    // Fetch intermediate data.
    cout << myname << "Fetch intermediate data." << endl;
    vector<const AdcSignalVector*> intSigs;
    vector<const AdcFlagVector*> intFlags;
    cout << myname << "  ... bad" << endl;
    auto iacd = intStates.dataMaps.find("bad");
    assert( iacd == intStates.dataMaps.end() );
    if  ( checkExtracted ) {
      cout << myname << "  ... extracted" << endl;
      iacd = intStates.dataMaps.find("extracted");
      assert( iacd != intStates.dataMaps.end() );
      intSigs.push_back(&iacd->second[chan].samples);
      intFlags.push_back(&iacd->second[chan].flags);
    }
    if ( checkMitigated ) {
      cout << myname << "  ... mitigated" << endl;
      iacd = intStates.dataMaps.find("mitigated");
      assert( iacd != intStates.dataMaps.end() );
      intSigs.push_back(&iacd->second[chan].samples);
      intFlags.push_back(&iacd->second[chan].flags);
    }
    if ( checkNoiseRemoved ) {
      cout << myname << "  ... noiseRemoved" << endl;
      iacd = intStates.dataMaps.find("noiseRemoved");
      assert( iacd != intStates.dataMaps.end() );
      intSigs.push_back(&iacd->second[chan].samples);
      intFlags.push_back(&iacd->second[chan].flags);
    }
    // Display results.
    cout << myname << "Display intermediate and final samples." << endl;
    for ( unsigned int isig=0; isig<nsig; ++isig ) {
      cout << setw(4) << chan << "-" 
           << setw(2) << isig << ": " << setw(4) << adcsmap[chan][isig];
      for ( unsigned int ista=0; ista<intSigs.size(); ++ista ) {
        cout << fixed << setprecision(1) << setw(8) << intSigs[ista]->at(isig);
        cout << " [" << intFlags[ista]->at(isig) << "]";
      }
      cout << fixed << setprecision(1) << setw(8) << sigs[isig]
           << " [" << flags[isig] << "]" << endl;
      assert( adcsmap[chan][isig] == acd.raw[isig] );
      if ( flags[isig] == AdcGood ) assert( sigequal(sigs[isig], sigsin[chan][isig]) );
      assert( flagequal(flags[isig], expflags[isig]) );
    }
  }

  cout << myname << line << endl;
  cout << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool useFclFile = true;   // If true ware are also test 35t dune reco fcl
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [UseExisting] [UseFclFile]" << endl;
      cout << "  If UseExisting = true, existing FCL file is used [false]." << endl;
      cout << "  If UseFclFile = true, FCL file dataprep_dune [true]." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    string sarg(argv[2]);
    useFclFile = sarg == "true" || sarg == "1";
  }
  return test_StandardRawDigitPrepService(useExistingFcl, useFclFile);
}

//**********************************************************************
