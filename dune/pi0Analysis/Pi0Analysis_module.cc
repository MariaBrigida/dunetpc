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
  void FillMcParticleVectors();
  int GetNHitsMatchingG4Id(std::vector<art::Ptr<recob::Hit>> pfpHits, art::Event const& e, int matchedMcParticleG4ID);
  double GetPfpCompleteness(int iShowerPfp);
  double GetPfpPurity(int iShowerPfp);
  void ClearAllVectors();
  int GetVectorPositionForTrackID(int trackID, std::vector<int> trackIDvector);

private:

  TTree *fEventTree;
  TTree *fPi0Tree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  //MC particles vectors
  std::vector<int> fMCParticleNHits, fMCParticleTrackID, fMCParticlePdgCode;
  std::vector<double> fMCParticleMass, fMCParticleEnergy, fMCParticleStartMomentumX, fMCParticleStartMomentumY, fMCParticleStartMomentumZ;
  std::vector<double> fMCParticleEndMomentumX, fMCParticleEndMomentumY, fMCParticleEndMomentumZ;
  std::vector<double> fMCParticleVertexX, fMCParticleVertexY, fMCParticleVertexZ;
  std::vector<double> fMCParticleEndX, fMCParticleEndY, fMCParticleEndZ;


  unsigned int fNMCParticles;
  unsigned int fNMCPhotons;
  unsigned int fNMCPi0s;

  unsigned int fNPFParticles;
  unsigned int fNShowerPFParticles;
  unsigned int fNReconstructedPhotons;
  unsigned int fNReconstructedPi0s;

  double fPhotonPairDistCenterToMCStart;
  std::vector<art::Ptr<recob::Hit> > fAllHits;
  std::map<int,std::vector<art::Ptr<recob::Hit>>> fTrueParticleHits, fTrueParticleHitsView0, fTrueParticleHitsView1, fTrueParticleHitsView2;
  art::Handle<std::vector<recob::Hit>> fHitHandle;
  art::Handle<std::vector<simb::MCParticle>> fMcParticles;

  std::vector<double> fPFPCompleteness, fPFPPurity, fShowerPFPCompleteness, fShowerPFPPurity;
  std::vector<double> fPFPShowerLength, fPFPShowerOpeningAngle, fPFPShowerNHits;
  std::vector<int> fPFPShowerTrueParticleMatchedId, fPFPShowerTrueParticleMatchedPosition, fPFPShowerTrueParticleMatchedPdgCode, fPFPShowerNSharedTrueParticleHits;
  
  //Pi0 vectors
  std::vector<int> fPi0Daughter1TrackID, fPi0Daughter2TrackID, fPi0Daughter1Position, fPi0Daughter2Position;
  std::vector<int> fPi0Daughter1ShowerPfpMatchedPosition, fPi0Daughter2ShowerPfpMatchedPosition;
 

  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;
  const geo::Geometry* fGeom;
  bool fRollUpUnsavedIDs;
  long unsigned int fNMinimumShowerHits;
  double fMinimumShowerEnergy;
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
  fMinimumShowerEnergy = p.get<long unsigned int>("MinimumShowerEnergy"); 

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::Pi0Analysis::analyze(art::Event const& e)
{
  ///////////////////Initialise/////////
  //Clear vectors
  test::Pi0Analysis::ClearAllVectors();
  //
  // Read event info
  fEventID = e.id().event();
  //std::cout << "fEventID = " << fEventID << std::endl;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  //
  //Fill vector of art::Ptr to access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    return;
  }
  //
  if(!e.isRealData()) {
     //Get MC Particles
     fMcParticles = e.getHandle<std::vector<simb::MCParticle>>(fTruthLabel);
     fNMCParticles = fMcParticles->size();
     //Fill mc to hit maps
     test::Pi0Analysis::FillMcToHitMaps(e);
     //Fill Mc particle vectors
     test::Pi0Analysis::FillMcParticleVectors();
   }
   /////////////////////////////////////

  //Loop over pfps
  int iPfp(0), iShowerPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
    //If it's a shower, fill shower plots////////////////////////////////
    if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp,e,fPFParticleLabel,fShowerLabel)){
      
      art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
      fPFPShowerLength.push_back(shower->Length());
      fPFPShowerOpeningAngle.push_back(shower->OpenAngle());
      std::vector<art::Ptr<recob::Hit>> pfpShowerHits;
      pfpShowerHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
      fPFPShowerNHits.push_back(pfpShowerHits.size());
      
      //Fill MC quantities such as best match true particle G4ID, purity, completeness 
      if(!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpShowerHits,fRollUpUnsavedIDs));
        fPFPShowerTrueParticleMatchedId.push_back(g4ID);
	//int pos(999999); for(unsigned int ipos=0; ipos<fNMCParticles; ipos++) {
	//  if(fMCParticleTrackID[ipos]==g4ID){pos=ipos; /*std::cout << "found pos = " << pos << std::endl;*/} 
          //std::cout << "fMCParticleTrackID[ipos] = " << fMCParticleTrackID[ipos] << std::endl;
 	//}
        //fPFPShowerTrueParticleMatchedPosition.push_back(pos);
        int pos = GetVectorPositionForTrackID(g4ID,fMCParticleTrackID);
        fPFPShowerTrueParticleMatchedPosition.push_back(pos);
        const simb::MCParticle trueParticle = fMcParticles->at(pos);
        int pfpShowerPdgCode = trueParticle.PdgCode();
        fPFPShowerTrueParticleMatchedPdgCode.push_back(pfpShowerPdgCode);
        int nSharedPfpShowerHits = test::Pi0Analysis::GetNHitsMatchingG4Id(pfpShowerHits,e,g4ID);
        fPFPShowerNSharedTrueParticleHits.push_back(nSharedPfpShowerHits);
        double comp = test::Pi0Analysis::GetPfpCompleteness(iShowerPfp);
        fShowerPFPCompleteness.push_back(comp);
        double purity = test::Pi0Analysis::GetPfpPurity(iShowerPfp);
        fShowerPFPPurity.push_back(purity);
        //std::cout << "iPfp = " << iPfp << " iShowerPfp = " << iShowerPfp << " g4ID = " << g4ID << " pos = " << pos << " pdg code = " << pfpShowerPdgCode << " shared hits = " << nSharedPfpShowerHits << " completeness = " << comp << " purity = " << purity << " length = " << shower->Length() << " opening angle = " << shower->OpenAngle() << std::endl;
      }
      iShowerPfp++;
    }
    ////////////////////////////////////////////////////////////////////
    iPfp++;
  }
  fNShowerPFParticles=iShowerPfp;
 
  //Access the truth information
  int nPi0s(0);
  if(!e.isRealData()){
    if(fMcParticles.isValid()){

     //Loop over MC particles and identify two photons coming from a pi0
      //bool isMCPrimary(false);
      //Loop over all MCparticles
      for(unsigned int iMc=0; iMc< fMcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = fMcParticles->at(iMc);
        int trueParticlePdgCode = trueParticle.PdgCode();
        //trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
        //If it's a photon, fill photon plots////////////////////////////////
         

          //If a pi0 is found
	  if(trueParticlePdgCode==111) {
            nPi0s++;
            std::cout << "--------------->Found a pi0<---------------" << std::endl;
            //fNReconstructedPhotons=0; std::vector<TVector3> photon1Direction, photon2Direction, photon1Start, photon2Start;
            std::cout << "N pi0 daughters = " << trueParticle.NumberDaughters() << std::endl;
            for(int iDaughter=0; iDaughter<trueParticle.NumberDaughters(); iDaughter++){
              int daughterID = trueParticle.Daughter(iDaughter);
              int daughterPosition = GetVectorPositionForTrackID(daughterID,fMCParticleTrackID);
              int daughterMatchedPfpPosition = GetVectorPositionForTrackID(daughterID,fPFPShowerTrueParticleMatchedId);
              if(iDaughter==0) {
                fPi0Daughter1TrackID.push_back(daughterID);
                fPi0Daughter1Position.push_back(daughterPosition);
                fPi0Daughter1ShowerPfpMatchedPosition.push_back(daughterMatchedPfpPosition);
                std::cout << "debug pi0 daughter 0: trackID/position/pfp position: " << daughterID << "/" << daughterPosition << "/" << daughterMatchedPfpPosition << std::endl;
              }
              else if(iDaughter==1){
                fPi0Daughter2TrackID.push_back(daughterID);
                fPi0Daughter2Position.push_back(daughterPosition);
                fPi0Daughter2ShowerPfpMatchedPosition.push_back(daughterMatchedPfpPosition);
                std::cout << "debug pi0 daughter 1: trackID/position/pfp position: " << daughterID << "/" << daughterPosition << "/" << daughterMatchedPfpPosition << std::endl;
              }
              //int daughterID = trueParticle.Daughter(iDaughter);
              //std::cout << "daughter n. " << iDaughter << " = " << daughterID <<std::endl;
    	    /*  std::vector<art::Ptr<recob::Hit>> pfpHits;
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
         */   
          }
	  fPi0Tree->Fill();
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
  }
  //std::cout << "fNMCPi0s = " << fNMCPi0s << std::endl;
  fNMCPi0s=nPi0s;
  fEventTree->Fill();
}

void test::Pi0Analysis::beginJob()
{
  // Implementation of optional member function here.
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event","Event Tree");
  fPi0Tree = tfs->make<TTree>("Pi0","Pi0 Tree");

  //Event tree
  fEventTree->Branch("eventID",&fEventID,"eventID/i");
  fEventTree->Branch("runID",&fRunID,"runID/i");
  fEventTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fEventTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fEventTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  ////MC particle branches
  fEventTree->Branch("mcParticleMass",&fMCParticleMass);
  fEventTree->Branch("mcParticleEnergy",&fMCParticleEnergy);
  fEventTree->Branch("mcParticleStartMomentumX",&fMCParticleStartMomentumX);
  fEventTree->Branch("mcParticleStartMomentumY",&fMCParticleStartMomentumY);
  fEventTree->Branch("mcParticleStartMomentumZ",&fMCParticleStartMomentumZ);
  fEventTree->Branch("mcParticleEndMomentumX",&fMCParticleEndMomentumX);
  fEventTree->Branch("mcParticleEndMomentumY",&fMCParticleEndMomentumY);
  fEventTree->Branch("mcParticleEndMomentumZ",&fMCParticleEndMomentumZ);
  fEventTree->Branch("mcParticleVertexX",&fMCParticleVertexX);
  fEventTree->Branch("mcParticleVertexY",&fMCParticleVertexY);
  fEventTree->Branch("mcParticleVertexZ",&fMCParticleVertexZ);
  fEventTree->Branch("mcParticleEndX",&fMCParticleEndX);
  fEventTree->Branch("mcParticleEndY",&fMCParticleEndY);
  fEventTree->Branch("mcParticleEndZ",&fMCParticleEndZ);

  fEventTree->Branch("nShowerPFParticles",&fNShowerPFParticles,"nShowerPFParticles/i");
  //fEventTree->Branch("pfpCompleteness",&fPFPCompleteness);
  fEventTree->Branch("showerPfpCompleteness",&fShowerPFPCompleteness);
  //fEventTree->Branch("pfpPurity",&fPFPPurity);
  fEventTree->Branch("showerPfpPurity",&fShowerPFPPurity);

  fEventTree->Branch("pfpShowerLength",&fPFPShowerLength);
  fEventTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpeningAngle);
  //fEventTree->Branch("pfpShowerTrueParticleMatchedId",&fPFPShowerTrueParticleMatchedId);
  //fEventTree->Branch("pfpShowerTrueParticleMatchedPdgCode",&fPFPShowerTrueParticleMatchedPdgCode);

  fEventTree->Branch("nMCPi0s",&fNMCPi0s,"nMCPi0s/i");
  fEventTree->Branch("nMCPhotons",&fNMCPhotons,"nMCPhotons/i");
  fEventTree->Branch("nReconstructedPhotons",&fNReconstructedPhotons,"nReconstructedPhotons/i");
  fEventTree->Branch("nReconstructedPi0s",&fNReconstructedPi0s,"nReconstructedPi0s/i");

  //Pi0 tree
  fPi0Tree->Branch("nPhotonPairDistCenterToMCStart",&fPhotonPairDistCenterToMCStart,"nReconstructedPhotons/D");
  fPi0Tree->Branch("pi0Daughter1TrackID", &fPi0Daughter1TrackID);
  fPi0Tree->Branch("pi0Daughter2TrackID", &fPi0Daughter2TrackID);
  fPi0Tree->Branch("pi0Daughter1Position", &fPi0Daughter1Position);
  fPi0Tree->Branch("pi0Daughter2Position", &fPi0Daughter2Position);
  fPi0Tree->Branch("pi0Daughter1ShowerPfpMatchedPosition", &fPi0Daughter1ShowerPfpMatchedPosition);
  fPi0Tree->Branch("pi0Daughter2ShowerPfpMatchedPosition", &fPi0Daughter2ShowerPfpMatchedPosition);

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

  if(fTrueParticleHits[g4Id].size()<fNMinimumShowerHits) return false;

  const simb::MCParticle trueParticle = fMcParticles->at(g4Id);
  if(trueParticle.E()<fMinimumShowerEnergy) return false;
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

void test::Pi0Analysis::FillMcParticleVectors(){
     for(unsigned int iMc=0; iMc<fMcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = fMcParticles->at(iMc);
        fMCParticleTrackID.push_back(trueParticle.TrackId());
        fMCParticleNHits.push_back(fTrueParticleHits[trueParticle.TrackId()].size());
        fMCParticlePdgCode.push_back(trueParticle.PdgCode());
        fMCParticleMass.push_back(trueParticle.Mass());
        fMCParticleEnergy.push_back(trueParticle.E());
        fMCParticleStartMomentumX.push_back(trueParticle.Px());
        fMCParticleStartMomentumY.push_back(trueParticle.Py());
        fMCParticleStartMomentumZ.push_back(trueParticle.Pz());
        fMCParticleEndMomentumX.push_back(trueParticle.EndPx());
        fMCParticleEndMomentumY.push_back(trueParticle.EndPy());
        fMCParticleEndMomentumZ.push_back(trueParticle.EndPz());
        fMCParticleVertexX.push_back(trueParticle.Vx());
        fMCParticleVertexY.push_back(trueParticle.Vy());
        fMCParticleVertexZ.push_back(trueParticle.Vz());
        fMCParticleEndX.push_back(trueParticle.EndX());
        fMCParticleEndY.push_back(trueParticle.EndY());
        fMCParticleEndZ.push_back(trueParticle.EndZ());
     }
} 


int test::Pi0Analysis::GetNHitsMatchingG4Id(std::vector<art::Ptr<recob::Hit>> pfpHits, art::Event const& e, int matchedMcParticleG4ID){
      int nPfpSharedTrueParticleHits(0);
      for (const auto& hit: pfpHits){
          auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);	//at some point try to avoid repeating this
          TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(clockData, hit, fRollUpUnsavedIDs));
          if (TruthMatchUtils::Valid(g4ID)){
              if(g4ID==matchedMcParticleG4ID)++nPfpSharedTrueParticleHits;
//              if(g4ID==fPFPShowerTrueParticleMatchedId[iPfp] && hit->View()==0)++fPFPShowerNSharedTrueParticleHitsView[iPfp][0];
//              if(g4ID==fPFPShowerTrueParticleMatchedId[iPfp] && hit->View()==1)++fPFPShowerNSharedTrueParticleHitsView[iPfp][1];
//             if(g4ID==fPFPShowerTrueParticleMatchedId[iPfp] && hit->View()==2)++fPFPShowerNSharedTrueParticleHitsView[iPfp][2];            
          }
      }
      return nPfpSharedTrueParticleHits;
}

double test::Pi0Analysis::GetPfpCompleteness(int iShowerPfp){
     if(fMCParticleNHits.at(fPFPShowerTrueParticleMatchedPosition[iShowerPfp]) > 0){
       return (float)fPFPShowerNSharedTrueParticleHits[iShowerPfp] / fMCParticleNHits[fPFPShowerTrueParticleMatchedPosition[iShowerPfp]];
     }
     return 0;
}

double test::Pi0Analysis::GetPfpPurity(int iPfp){
      if(fPFPShowerNHits.at(iPfp) > 0) return (float)fPFPShowerNSharedTrueParticleHits[iPfp] / fPFPShowerNHits.at(iPfp);
    return 0;
}

void test::Pi0Analysis::ClearAllVectors(){
      //MC particle vectors
      fMCParticleNHits.clear();
      fMCParticleTrackID.clear();
      fMCParticleMass.clear();
      fMCParticleEnergy.clear();
      fMCParticleStartMomentumX.clear();
      fMCParticleStartMomentumY.clear();
      fMCParticleStartMomentumZ.clear();
      fMCParticleEndMomentumX.clear();
      fMCParticleEndMomentumY.clear();
      fMCParticleEndMomentumZ.clear();
      fMCParticleVertexX.clear();
      fMCParticleVertexY.clear();
      fMCParticleVertexZ.clear();
      fMCParticleEndX.clear();
      fMCParticleEndY.clear();
      fMCParticleEndZ.clear();

      fAllHits.clear();
      fTrueParticleHits.clear();
      fTrueParticleHitsView0.clear();
      fTrueParticleHitsView1.clear();
      fTrueParticleHitsView2.clear();
      fPFPCompleteness.clear();
      fShowerPFPCompleteness.clear();
      fPFPPurity.clear();
      fShowerPFPPurity.clear();
      fPFPShowerLength.clear();
      fPFPShowerOpeningAngle.clear();
      fPFPShowerNHits.clear();
      fPFPShowerTrueParticleMatchedId.clear();
      fPFPShowerTrueParticleMatchedPosition.clear();
      fPFPShowerTrueParticleMatchedPdgCode.clear();
      fPFPShowerNSharedTrueParticleHits.clear();
      //Pi0 vectors
      fPi0Daughter1TrackID.clear();
      fPi0Daughter2TrackID.clear();
      fPi0Daughter1Position.clear();
      fPi0Daughter2Position.clear();
      fPi0Daughter1ShowerPfpMatchedPosition.clear();
      fPi0Daughter2ShowerPfpMatchedPosition.clear();

}

int test::Pi0Analysis::GetVectorPositionForTrackID(int trackID, std::vector<int> trackIDvector){
      int pos(999999);
      for(unsigned int ipos=0; ipos<trackIDvector.size(); ipos++) {
        if(trackIDvector.at(ipos)==trackID){pos=ipos;} 
      }
      return pos;
}

DEFINE_ART_MODULE(test::Pi0Analysis)
