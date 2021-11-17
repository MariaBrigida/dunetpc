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
  const char * fileNameChar = fTestFileName.c_str();
  fTestFile = TFile::Open(fileNameChar);
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
  TH1D* hTruePi0InvMass;
  hTruePi0InvMass = tfs->make<TH1D>("hTruePi0InvMass","hTruePi0InvMass",100,0,200);

  TTree *eventTree = (TTree*)fTestFile->Get("ana/Event");
  //eventTree->Print();
  //std::vector<double> *mcParticleVertexX = 0;
  //TBranch *b_mcParticleVertexX;
  //tree->SetBranchAddress("mcParticleVertexX",&mcParticleVertexX, &b_mcParticleVertexX);
 
  unsigned int globalEventID; 
  std::vector<int> *mcParticleTrackID = 0;
  std::vector<int> *mcParticlePdgCode = 0;
  std::vector<int> *mcParticleMotherPdgCode = 0;
  std::vector<int> *mcParticleMotherPosition = 0;
  std::vector<double> *mcParticleStartMomentumX = 0;
  std::vector<double> *mcParticleStartMomentumY = 0;
  std::vector<double> *mcParticleStartMomentumZ = 0;
  std::vector<double> *mcParticleEnergy = 0;
  
  //Define Event tree branches

  TBranch *b_mcParticleTrackID, *b_mcParticlePdgCode, *b_mcParticleMotherPdgCode, *b_mcParticleMotherPosition;
  TBranch *b_mcParticleEnergy;
  TBranch *b_mcParticleStartMomentumX, *b_mcParticleStartMomentumY, *b_mcParticleStartMomentumZ;


  //Set Event tree branch addresses
  eventTree->SetBranchAddress("globalEventID",&globalEventID);
  eventTree->SetBranchAddress("mcParticleTrackID",&mcParticleTrackID, &b_mcParticleTrackID);
  eventTree->SetBranchAddress("mcParticlePdgCode",&mcParticlePdgCode, &b_mcParticlePdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPdgCode",&mcParticleMotherPdgCode, &b_mcParticleMotherPdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPosition",&mcParticleMotherPosition, &b_mcParticleMotherPosition);
  eventTree->SetBranchAddress("mcParticleStartMomentumX",&mcParticleStartMomentumX, &b_mcParticleStartMomentumX);
  eventTree->SetBranchAddress("mcParticleStartMomentumY",&mcParticleStartMomentumY, &b_mcParticleStartMomentumY);
  eventTree->SetBranchAddress("mcParticleStartMomentumZ",&mcParticleStartMomentumZ, &b_mcParticleStartMomentumZ);
  eventTree->SetBranchAddress("mcParticleEnergy",&mcParticleEnergy, &b_mcParticleEnergy);


  //Pi0 vectors
  std::vector<double> pi0Daughter1TrueEnergy, pi0Daughter1TrueDirectionX, pi0Daughter1TrueDirectionY, pi0Daughter1TrueDirectionZ;
  std::vector<double> pi0Daughter2TrueEnergy, pi0Daughter2TrueDirectionX, pi0Daughter2TrueDirectionY, pi0Daughter2TrueDirectionZ;

  //Loop over eventTree entries
  for(int iEntry=0; iEntry<eventTree->GetEntries(); iEntry++){
    eventTree->GetEntry(iEntry);
    std::cout << "mcParticleTrackID size = " << mcParticleTrackID->size() << std::endl;

    int pi0Daughter1MotherPos(-999), pi0Daughter2MotherPos(-999);
    bool daughter1Found(false), daughter2Found(false);
    for(std::vector<int>::size_type iMCPart = 0; iMCPart<mcParticleTrackID->size(); iMCPart++){
      //std::cout << "global event ID = " << globalEventID << " iMCPart = " << iMCPart << " mcParticleTrackID = " << mcParticleTrackID->at(iMCPart) << " mcParticlePdgCode = " << mcParticlePdgCode->at(iMCPart) << " mcParticleMotherPdgCode = " << mcParticleMotherPdgCode->at(iMCPart) << std::endl;
      if(mcParticlePdgCode->at(iMCPart)==22 && mcParticleMotherPdgCode->at(iMCPart)==111){
        //std::cout << "Found a pi0 at position = " << iMCPart << std::endl;
        if(!daughter1Found && !daughter2Found) {
          pi0Daughter1MotherPos = mcParticleMotherPosition->at(iMCPart);
          pi0Daughter1TrueEnergy.push_back(mcParticleEnergy->at(iMCPart));
          pi0Daughter1TrueDirectionX.push_back(mcParticleStartMomentumX->at(iMCPart)); 
          pi0Daughter1TrueDirectionY.push_back(mcParticleStartMomentumY->at(iMCPart)); 
          pi0Daughter1TrueDirectionZ.push_back(mcParticleStartMomentumZ->at(iMCPart));
          daughter1Found=true;
          //std::cout << "Found first daughter photon at position = " << iMCPart << ". its mother position is " << mcParticleMotherPosition->at(iMCPart) << " energy = " << pi0Daughter1TrueEnergy.back() << " mom = " << pi0Daughter1TrueDirectionX.back() << " " << pi0Daughter1TrueDirectionY.back() << " " << pi0Daughter1TrueDirectionZ.back() << std::endl; 
            //std::cout << "debug = mcParticleStartMomentumX->at(iMCPart) = " << mcParticleStartMomentumX->at(iMCPart) << " " << mcParticleStartMomentumY->at(iMCPart) << " " << mcParticleStartMomentumZ->at(iMCPart) << std::endl; 
        }
        else if(daughter1Found && !daughter2Found) {
          pi0Daughter2MotherPos = mcParticleMotherPosition->at(iMCPart);
          if(pi0Daughter2MotherPos!=pi0Daughter1MotherPos){
            pi0Daughter1TrueEnergy.pop_back();
            pi0Daughter1TrueDirectionX.pop_back();
            pi0Daughter1TrueDirectionY.pop_back();
            pi0Daughter1TrueDirectionZ.pop_back();
            pi0Daughter1TrueEnergy.push_back(mcParticleEnergy->at(iMCPart));
            pi0Daughter1TrueDirectionX.push_back(mcParticleStartMomentumX->at(iMCPart));
            pi0Daughter1TrueDirectionY.push_back(mcParticleStartMomentumY->at(iMCPart));
            pi0Daughter1TrueDirectionZ.push_back(mcParticleStartMomentumZ->at(iMCPart));
            pi0Daughter1MotherPos = mcParticleMotherPosition->at(iMCPart);
            //std::cout << "Previous particle was from dalitz. Replaced. Found first daughter photon at position = " << iMCPart << ". its mother position is " << mcParticleMotherPosition->at(iMCPart) << " energy = " << pi0Daughter1TrueEnergy.back() << " mom = " << pi0Daughter1TrueDirectionX.back() << " " << pi0Daughter1TrueDirectionY.back() << " " << pi0Daughter1TrueDirectionZ.back() << std::endl;
            //std::cout << "debug = mcParticleStartMomentumX->at(iMCPart) = " << mcParticleStartMomentumX->at(iMCPart) << " " << mcParticleStartMomentumY->at(iMCPart) << " " << mcParticleStartMomentumZ->at(iMCPart) << std::endl; 
          }
          else{
            pi0Daughter2TrueEnergy.push_back(mcParticleEnergy->at(iMCPart));
            pi0Daughter2TrueDirectionX.push_back(mcParticleStartMomentumX->at(iMCPart)); 
            pi0Daughter2TrueDirectionY.push_back(mcParticleStartMomentumY->at(iMCPart)); 
            pi0Daughter2TrueDirectionZ.push_back(mcParticleStartMomentumZ->at(iMCPart)); 
            daughter2Found=true; 
            pi0Daughter2MotherPos = mcParticleMotherPosition->at(iMCPart);
            //std::cout << "Found second daughter photon at position = " << iMCPart << ". its mother position is " << mcParticleMotherPosition->at(iMCPart) << " energy = " << pi0Daughter2TrueEnergy.back() << " mom = " << pi0Daughter2TrueDirectionX.back() << " " << pi0Daughter2TrueDirectionY.back() << " " << pi0Daughter2TrueDirectionZ.back() << std::endl; 
            //std::cout << "debug = mcParticleStartMomentumX->at(iMCPart) = " << mcParticleStartMomentumX->at(iMCPart) << " " << mcParticleStartMomentumY->at(iMCPart) << " " << mcParticleStartMomentumZ->at(iMCPart) << std::endl; 
         }
        }
        std::cout << "WARNING: daughter1Found = " << daughter1Found << " daughter2Found = " << daughter2Found << std::endl;
        if(daughter1Found && daughter2Found){daughter1Found=false; daughter2Found=false;}

      }

    }
  }
  
  std::cout << "sizes. pi0Daughter1TrueEnergy = " << pi0Daughter1TrueEnergy.size() << " pi0Daughter1TrueDirection: " << pi0Daughter1TrueDirectionX.size() << " " << pi0Daughter1TrueDirectionY.size() << " " << pi0Daughter1TrueDirectionZ.size() << " pi0Daughter2TrueEnergy = " << pi0Daughter2TrueEnergy.size() << " pi0Daughter2TrueDirection: " << pi0Daughter2TrueDirectionX.size() << " " << pi0Daughter2TrueDirectionY.size() << " " << pi0Daughter2TrueDirectionZ.size() << std::endl;

  //Make pi0 mass plots
  for(long unsigned int iPi0=0; iPi0<pi0Daughter1TrueEnergy.size(); iPi0++){

    TVector3 phot1mom(pi0Daughter1TrueDirectionX.at(iPi0),pi0Daughter1TrueDirectionY.at(iPi0),pi0Daughter1TrueDirectionZ.at(iPi0));
    TVector3 phot2mom(pi0Daughter2TrueDirectionX.at(iPi0),pi0Daughter2TrueDirectionY.at(iPi0),pi0Daughter2TrueDirectionZ.at(iPi0));
    double angle = phot1mom.Angle(phot2mom);
    double invMass = TMath::Sqrt(2*pi0Daughter1TrueEnergy.at(iPi0)*pi0Daughter2TrueEnergy.at(iPi0)*(1-TMath::Cos(angle))); 
    std::cout << "Pi0 n. " << iPi0 << " phot 1: E Mom: " << pi0Daughter1TrueEnergy.at(iPi0) << " " << pi0Daughter1TrueDirectionX.at(iPi0) << " " << pi0Daughter1TrueDirectionY.at(iPi0) << " " << pi0Daughter1TrueDirectionZ.at(iPi0) << std::endl;
    std::cout << "Pi0 n. " << iPi0 << " phot 2: E Mom: " << pi0Daughter2TrueEnergy.at(iPi0) << " " << pi0Daughter2TrueDirectionX.at(iPi0) << " " << pi0Daughter2TrueDirectionY.at(iPi0) << " " << pi0Daughter2TrueDirectionZ.at(iPi0) << std::endl;
    std::cout << "pi0 inv mass = " << invMass << std::endl;
    hTruePi0InvMass->Fill(invMass);
  }


  //pi0Tree->Print();
/*
  TTree *pi0Tree = (TTree*)fTestFile->Get("ana/Pi0");
  std::cout << "deb3" << std::endl;

  //ROOT::Experimental::TDataFrame d("ana/Event", fTestFile);
  //Read pi0 branches
  unsigned int pi0Daughter1TrackID, pi0Daughter2TrackID, pi0GlobalEventID;
  unsigned int pi0Daughter1Position, pi0Daughter2Position;
  pi0Tree->SetBranchAddress("pi0GlobalEventID",&pi0GlobalEventID);
  pi0Tree->SetBranchAddress("pi0Daughter1TrackID",&pi0Daughter1TrackID);
  pi0Tree->SetBranchAddress("pi0Daughter2TrackID",&pi0Daughter2TrackID);
  pi0Tree->SetBranchAddress("pi0Daughter1Position", &pi0Daughter1Position);
  pi0Tree->SetBranchAddress("pi0Daughter2Position", &pi0Daughter2Position);

  //Read Event branches
  unsigned int globalEventID;
  eventTree->SetBranchAddress("globalEventID",&globalEventID);

  for(int iEntry=0; iEntry<pi0Tree->GetEntries(); iEntry++){
    pi0Tree->GetEntry(iEntry);
    std::cout << "pi0GlobalEventID = " << pi0GlobalEventID << " pi0Daughter1TrackID = " << pi0Daughter1TrackID << " pi0Daughter2TrackID = " << pi0Daughter2TrackID << std::endl;
    std::cout << "pi0Daughter1Position = " << pi0Daughter1Position << " pi0Daughter2Position = " << pi0Daughter2Position << std::endl;

   //Lambda to cut on global event number
    //auto globalEvCut = [](double x) { return x == pi0GlobalEventID; };
    
   
  }
*/

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
