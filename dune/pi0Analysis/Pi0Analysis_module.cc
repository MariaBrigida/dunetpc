////////////////////////////////////////////////////////////////////////
// Class:       Pi0Analysis
// Plugin Type: analyzer (art v3_06_03)
// File:        Pi0Analysis_module.cc
//
// Generated at Mon Sep 27 11:01:20 2021 by Maria Brigida Brunetti using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace test {
  class Pi0Analysis;
}


class test::Pi0Analysis : public art::EDAnalyzer {
public:
  explicit Pi0Analysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0Analysis(Pi0Analysis const&) = delete;
  Pi0Analysis(Pi0Analysis&&) = delete;
  Pi0Analysis& operator=(Pi0Analysis const&) = delete;
  Pi0Analysis& operator=(Pi0Analysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


test::Pi0Analysis::Pi0Analysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::Pi0Analysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void test::Pi0Analysis::beginJob()
{
  // Implementation of optional member function here.
}

void test::Pi0Analysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::Pi0Analysis)
