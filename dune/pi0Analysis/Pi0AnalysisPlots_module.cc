////////////////////////////////////////////////////////////////////////
// Class:       Pi0AnalysisPlots
// Plugin Type: analyzer (art v3_06_03)
// File:        Pi0AnalysisPlots_module.cc
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"
#include "art_root_io/TFileService.h"

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
  class Pi0AnalysisPlots;
}


class test::Pi0AnalysisPlots : public art::EDAnalyzer {
public:
  explicit Pi0AnalysisPlots(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Pi0AnalysisPlots(Pi0AnalysisPlots const&) = delete;
  Pi0AnalysisPlots(Pi0AnalysisPlots&&) = delete;
  Pi0AnalysisPlots& operator=(Pi0AnalysisPlots const&) = delete;
  Pi0AnalysisPlots& operator=(Pi0AnalysisPlots&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  TTree *fTestTree;
  std::string fTestFileName;
  TFile *fTestFile;
  std::vector<double> fMcParticleVertexX_test;
};


test::Pi0AnalysisPlots::Pi0AnalysisPlots(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fTestFileName = p.get<std::string>("TestFileName");
  std::cout << "deb0" << std::endl;
  const char * fileNameChar = fTestFileName.c_str();
  std::cout << "deb1" << std::endl;
  fTestFile = TFile::Open(fileNameChar);
  std::cout << "deb2" << std::endl;
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::Pi0AnalysisPlots::analyze(art::Event const& e)
{
}

void test::Pi0AnalysisPlots::beginJob()
{
  
  art::ServiceHandle<art::TFileService> tfs;
  // Implementation of optional member function here.
  fTestTree = tfs->make<TTree>("Test","Test Tree");


  TTree *pi0Tree = (TTree*)fTestFile->Get("ana/Pi0");
  //TTree *eventTree = (TTree*)fTestFile->Get("ana/Event");

  //pi0Tree->Print();

  std::cout << "deb3" << std::endl;

  //Read pi0 branches
  unsigned int pi0Daughter1TrackID, pi0Daughter2TrackID, pi0GlobalEventID;
  unsigned int pi0Daughter1Position, pi0Daughter2Position;
  pi0Tree->SetBranchAddress("pi0GlobalEventID",&pi0GlobalEventID);
  pi0Tree->SetBranchAddress("pi0Daughter1TrackID",&pi0Daughter1TrackID);
  pi0Tree->SetBranchAddress("pi0Daughter2TrackID",&pi0Daughter2TrackID);
  pi0Tree->SetBranchAddress("pi0Daughter1Position", &pi0Daughter1Position);
  pi0Tree->SetBranchAddress("pi0Daughter2Position", &pi0Daughter2Position);
  for(int iEntry=0; iEntry<pi0Tree->GetEntries(); iEntry++){
    pi0Tree->GetEntry(iEntry);
    std::cout << "pi0GlobalEventID = " << pi0GlobalEventID << " pi0Daughter1TrackID = " << pi0Daughter1TrackID << " pi0Daughter2TrackID = " << pi0Daughter2TrackID << std::endl;
    std::cout << "pi0Daughter1Position = " << pi0Daughter1Position << " pi0Daughter2Position = " << pi0Daughter2Position << std::endl;
    
  }


/*
  //TestReadFile
  //fTestTree->Branch("mcParticleVertexX_test", &fMcParticleVertexX_test);

  std::cout << "deb5" << std::endl;
  std::vector<double> *mcParticleVertexX = 0;
  std::cout << "deb6" << std::endl;
  //tree->ls();
  TBranch *b_mcParticleVertexX;
  tree->SetBranchAddress("mcParticleVertexX",&mcParticleVertexX, &b_mcParticleVertexX);
  std::cout << "debug 61" << std::endl;
  //tree->Print();
  std::cout << "deb7. tree->GetEntries() = " << tree->GetEntries() << std::endl;
  for(int iEntry=0; iEntry<tree->GetEntries(); iEntry++){
    std::cout << "deb 70" << std::endl;
    tree->GetEntry(iEntry);
    std::cout << "debug mcParticleVertexX size = " << mcParticleVertexX->size() << std::endl;
    for(std::vector<double>::size_type iVect=0; iVect<mcParticleVertexX->size(); iVect++){
      fMcParticleVertexX_test.push_back(mcParticleVertexX->at(iVect));
    }
    fTestTree->Fill();
    fMcParticleVertexX_test.clear();
  }
  std::cout << "deb8" << std::endl;
*/
}

void test::Pi0AnalysisPlots::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::Pi0AnalysisPlots)
