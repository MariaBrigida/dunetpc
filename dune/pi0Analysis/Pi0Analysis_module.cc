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
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


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
  TVector3 centerOfClosestApproachPoint(TVector3 startPoint1, TVector3 direction1,TVector3 startPoint2, TVector3 direction2);
  bool photonIsReconstructable(TruthMatchUtils::G4ID g4Id);
  void FillMcToHitMaps(art::Event const& e);
private:

  TTree *fEventTree;
  TTree *fPi0Tree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  unsigned int fNMCParticles;
  unsigned int fNPFParticles;
  unsigned int fNReconstructedPhotons;
  double fPhotonPairDistCenterToMCStart;
  std::vector<art::Ptr<recob::Hit> > fAllHits;
  std::map<int,std::vector<art::Ptr<recob::Hit>>> fTrueParticleHits, fTrueParticleHitsView0, fTrueParticleHitsView1, fTrueParticleHitsView2;
  art::Handle<std::vector<recob::Hit>> fHitHandle;


  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;
  const geo::Geometry* fGeom;
  bool fRollUpUnsavedIDs;
  long unsigned int fNMinimumShowerHits;
  // Declare member data here.

};


test::Pi0Analysis::Pi0Analysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitLabel = p.get<std::string>("HitLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs"); 
  fNMinimumShowerHits = p.get<long unsigned int>("NMinimumShowerHits"); 

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::Pi0Analysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  std::cout << "fEventID = " << fEventID << std::endl;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  //Access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    return;
  }
  
  //Access the truth information
  if(!e.isRealData()){
    test::Pi0Analysis::FillMcToHitMaps(e);
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
     fNMCParticles=mcParticles->size();
     //Loop over MC particles and identify two photons coming from a pi0
      bool isMCPrimary(false);
      for(unsigned int iMc=0; iMc< mcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        int trueParticlePdgCode = trueParticle.PdgCode();
        trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
	  if(trueParticlePdgCode==111) {
            std::cout << "--------------->Found a pi0<---------------" << std::endl;
            fNReconstructedPhotons=0; std::vector<TVector3> photon1Direction, photon2Direction, photon1Start, photon2Start;
            for(int iDaughter=0; iDaughter<trueParticle.NumberDaughters(); iDaughter++){
              int daughterID = trueParticle.Daughter(iDaughter);
              std::cout << "daughter n. " << iDaughter << " = " << daughterID <<std::endl;
    	      std::vector<art::Ptr<recob::Hit>> pfpHits;
              for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
                if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp,e,fPFParticleLabel,fShowerLabel)){
                  art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
                  pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
                //Find if there is a reco particle that matches best this photon
    		  std::vector<art::Ptr<recob::Hit>> pfpHits;
      		  pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
        	  TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpHits,fRollUpUnsavedIDs));
		  if(g4ID==daughterID) {
                    std::cout << "Found a reco particle that matches photon with ID = " << daughterID << std::endl;
                    if (test::Pi0Analysis::photonIsReconstructable(g4ID)){
		    //if(pfpHits.size()>fNMinimumShowerHits){
                      if(daughterID==0){photon1Direction.push_back(shower->Direction()); photon1Start.push_back(shower->ShowerStart());}
                      if(daughterID==1){photon2Direction.push_back(shower->Direction()); photon2Start.push_back(shower->ShowerStart());}
                    }
                    fNReconstructedPhotons++;
                  }
		}
	      }

            }
	      //Loop over pairs of shower directions for this pi0
	      for(long unsigned int iShower1=0; iShower1<photon1Direction.size(); iShower1++){
	        for(long unsigned int iShower2=0; iShower2<photon2Direction.size(); iShower2++){
		  TVector3 closestApproachPoint =  test::Pi0Analysis::centerOfClosestApproachPoint(photon1Start.at(iShower1),photon1Direction.at(iShower1),photon2Start.at(iShower2),photon2Direction.at(iShower2));
	          std::cout << "TEST: middle point of closest approach between shower directions = " << closestApproachPoint.X() << " " << closestApproachPoint.Y() << " " << closestApproachPoint.Z()  << std::endl;
                  fPhotonPairDistCenterToMCStart= (closestApproachPoint-trueParticle.Position().Vect()).Mag();
		  std::cout << "fPhotonPairDistCenterToMCStart = " << fPhotonPairDistCenterToMCStart << std::endl;
	        }
	      }

	    if(fNReconstructedPhotons==2){
            }
	    fPi0Tree->Fill();
          }
	}
	/*if(trueParticlePdgCode==111) {
          double pi0VertexX = trueParticle.Position().X();
          double pi0VertexY = trueParticle.Position().Y();
          double pi0VertexZ = trueParticle.Position().Z();
          //Check if there is any particle with vertex close to pi0 vertex
          for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
            if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp,e,fPFParticleLabel,fShowerLabel)){
              art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
              //double showerStartX = shower->ShowerStart().X();
              //double showerStartY = shower->ShowerStart().Y();
              //double showerStartZ = shower->ShowerStart().Z();
              (trueParticle.Position().Vect()-shower->ShowerStart().Vect()).Mag();
            }
	  }
          std::cout << "--------------->Found a pi0<---------------" << std::endl;
          std::cout << "pi0VertexX = " << pi0VertexX << "pi0VertexY = " << pi0VertexY << " pi0VertexZ = " << pi0VertexZ << std::endl;
          std::cout << "trueParticle.PdgCode() = " << trueParticle.PdgCode() << " is primary? " << isMCPrimary << std::endl;
          std::cout << "Start momentum: " << trueParticle.Momentum().X() << " " << trueParticle.Momentum().Y() << " " << trueParticle.Momentum().Z() << " " << std::endl;
          //std::cout << "End momentum: " << trueParticle.EndMomentum().X() << " " << trueParticle.EndMomentum().Y() << " " << trueParticle.EndMomentum().Z() << " " << std::endl;
          
        }*/


    }
  }
  fEventTree->Fill();
}

void test::Pi0Analysis::beginJob()
{
  // Implementation of optional member function here.
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event","Event Tree");
  fPi0Tree = tfs->make<TTree>("Pi0","Pi0 Tree");

  //Event branches
  fEventTree->Branch("eventID",&fEventID,"eventID/i");
  fEventTree->Branch("runID",&fRunID,"runID/i");
  fEventTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fEventTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fEventTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  fEventTree->Branch("nReconstructedPhotons",&fNReconstructedPhotons,"nReconstructedPhotons/i");
  fPi0Tree->Branch("nPhotonPairDistCenterToMCStart",&fPhotonPairDistCenterToMCStart,"nReconstructedPhotons/D");
}

void test::Pi0Analysis::endJob()
{
  // Implementation of optional member function here.
}

TVector3 test::Pi0Analysis::centerOfClosestApproachPoint(TVector3 startPoint1, TVector3 direction1, TVector3 startPoint2, TVector3 direction2){

  /*std::cout << "DEBUG centerOfClosestApproach" << std::endl;
  std::cout << "direction1: " << direction1.X() << " " << direction1.Y() << " " << direction1.Z() << std::endl;
  std::cout << "direction2: " << direction2.X() << " " << direction2.Y() << " " << direction2.Z() << std::endl;
  std::cout << "startPoint1: " << startPoint1.X() << " " << startPoint1.Y() << " " << startPoint1.Z() << std::endl;
  std::cout << "startPoint2: " << startPoint2.X() << " " << startPoint2.Y() << " " << startPoint2.Z() << std::endl;*/

  TVector3 perpDir = direction1.Cross(direction2);
  TVector3 perpDir1 = direction1.Cross(perpDir);
  TVector3 perpDir2 = direction2.Cross(perpDir);
  double factor1 = perpDir2.Dot((startPoint2 - startPoint1))/perpDir2.Dot(direction1);
  double factor2 = perpDir1.Dot((startPoint1 - startPoint2))/perpDir1.Dot(direction2);
  TVector3 clostestApproachPoint1 = startPoint1 + factor1*direction1;
  TVector3 clostestApproachPoint2 = startPoint2 + factor2*direction2;
  TVector3 middlePoint = (clostestApproachPoint1 + clostestApproachPoint2)*0.5;
  return middlePoint;
}


bool test::Pi0Analysis::photonIsReconstructable(TruthMatchUtils::G4ID g4Id){

//  const simb::MCParticle trueParticle = mcParticles->at(g4Id);
  if(fTrueParticleHits[g4Id].size()<fNMinimumShowerHits) return false;
  return true;
}

void test::Pi0Analysis::FillMcToHitMaps(art::Event const& e){
    //Get all hits for all MCparticles in event
    if(e.getByLabel(fHitLabel,fHitHandle))
    {art::fill_ptr_vector(fAllHits, fHitHandle);}

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    for (const auto& hit: fAllHits){
      TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
      if (TruthMatchUtils::Valid(g4ID)){
          fTrueParticleHits[g4ID].push_back(hit);
          if(hit->View()==0)fTrueParticleHitsView0[g4ID].push_back(hit);
          if(hit->View()==1)fTrueParticleHitsView1[g4ID].push_back(hit);
          if(hit->View()==2)fTrueParticleHitsView2[g4ID].push_back(hit);
      }
    }
}

DEFINE_ART_MODULE(test::Pi0Analysis)
