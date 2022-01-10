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

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;
  unsigned int fGlobalEventID;

  //MC particles vectors
  std::vector<int> fMCParticleNHits, fMCParticleNHitsU, fMCParticleNHitsV, fMCParticleNHitsW;
  std::vector<int> fMCParticleTrackID, fMCParticleMotherID, fMCParticleMotherPosition, fMCParticleMotherPdgCode, fMCParticlePdgCode;
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

  std::vector<double> fPFPCompleteness, fPFPPurity, fPFPShowerCompleteness, fPFPShowerPurity;
  std::vector<double> fPFPShowerLength, fPFPShowerOpeningAngle, fPFPShowerNHits, fPFPShowerNHitsU, fPFPShowerNHitsV, fPFPShowerNHitsW, fPFPShowerCollectionPlaneEnergy;
  std::vector<double> fPFPShowerDirectionX, fPFPShowerDirectionY, fPFPShowerDirectionZ; //direction cosines at start of shower
  std::vector<int> fPFPShowerTrueParticleMatchedId, fPFPShowerTrueParticleMatchedPosition, fPFPShowerTrueParticleMatchedPdgCode, fPFPShowerNSharedTrueParticleHits;
  std::vector<int> fPFPShowerTrueParticleNHits, fPFPShowerTrueParticleNHitsU, fPFPShowerTrueParticleNHitsV, fPFPShowerTrueParticleNHitsW;
  
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
  //std::cout << "------------------------------------------------------" << std::endl;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  //
  //Fill vector of art::Ptr to access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    //std::cout << "No PFParticles found!" << std::endl;
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
      //std::cout << "shower energy size = " << shower->Energy().size() << std::endl;
      if(shower->Energy().size()==3)fPFPShowerCollectionPlaneEnergy.push_back(shower->Energy().at(2));  //collection plane energy
      fPFPShowerDirectionX.push_back(shower->Direction().X());
      fPFPShowerDirectionY.push_back(shower->Direction().Y());
      fPFPShowerDirectionZ.push_back(shower->Direction().Z());
      fPFPShowerOpeningAngle.push_back(shower->OpenAngle());
      std::vector<art::Ptr<recob::Hit>> pfpShowerHits;
      int nHitsViewU(0), nHitsViewV(0), nHitsViewW(0);
      pfpShowerHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
      for(auto hit:pfpShowerHits){
        if(hit->View()==geo::kU)nHitsViewU++;
        else if(hit->View()==geo::kV)nHitsViewV++;
        else if(hit->View()==geo::kW)nHitsViewW++;
      }
      fPFPShowerNHitsU.push_back(nHitsViewU);
      fPFPShowerNHitsV.push_back(nHitsViewV);
      fPFPShowerNHitsW.push_back(nHitsViewW);
      fPFPShowerNHits.push_back(pfpShowerHits.size());
      
      //Fill MC quantities such as best match true particle G4ID, purity, completeness 
      if (!e.isRealData()) {
        TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,pfpShowerHits,fRollUpUnsavedIDs));
        std::cout << "iPfp = " << iPfp << "/" << pfparticleVect.size() << "  shower TrueParticleIDFromTotalRecoHits = " << g4ID << std::endl;
	//int pos(999999); for(unsigned int ipos=0; ipos<fNMCParticles; ipos++) {
	//  if(fMCParticleTrackID[ipos]==g4ID){pos=ipos; /*std::cout << "found pos = " << pos << std::endl;*/} 
          //std::cout << "fMCParticleTrackID[ipos] = " << fMCParticleTrackID[ipos] << std::endl;
 	//}
        //fPFPShowerTrueParticleMatchedPosition.push_back(pos);
        int pos = GetVectorPositionForTrackID(g4ID,fMCParticleTrackID);
        if(pos<900000) {
            fPFPShowerTrueParticleMatchedPosition.push_back(pos);
            fPFPShowerTrueParticleMatchedId.push_back(g4ID);
            const simb::MCParticle trueParticle = fMcParticles->at(pos);
            int pfpShowerPdgCode = trueParticle.PdgCode();
            fPFPShowerTrueParticleMatchedPdgCode.push_back(pfpShowerPdgCode);
            fPFPShowerTrueParticleNHits.push_back(fMCParticleNHits.at(pos));
            fPFPShowerTrueParticleNHitsU.push_back(fMCParticleNHitsU.at(pos));
            fPFPShowerTrueParticleNHitsV.push_back(fMCParticleNHitsV.at(pos));
            fPFPShowerTrueParticleNHitsW.push_back(fMCParticleNHitsW.at(pos));
            int nSharedPfpShowerHits = test::Pi0Analysis::GetNHitsMatchingG4Id(pfpShowerHits,e,g4ID);
            fPFPShowerNSharedTrueParticleHits.push_back(nSharedPfpShowerHits);
            double comp = test::Pi0Analysis::GetPfpCompleteness(iShowerPfp);
            fPFPShowerCompleteness.push_back(comp);
            double purity = test::Pi0Analysis::GetPfpPurity(iShowerPfp);
            fPFPShowerPurity.push_back(purity);
        }
        else{
            fPFPShowerTrueParticleMatchedPosition.push_back(999999);
            fPFPShowerTrueParticleMatchedId.push_back(g4ID);
            fPFPShowerTrueParticleMatchedPdgCode.push_back(999999);
            fPFPShowerTrueParticleNHits.push_back(999999);
            fPFPShowerTrueParticleNHitsU.push_back(999999);
            fPFPShowerTrueParticleNHitsV.push_back(999999);
            fPFPShowerTrueParticleNHitsW.push_back(999999);
            fPFPShowerNSharedTrueParticleHits.push_back(999999);
            fPFPShowerCompleteness.push_back(999999);
            fPFPShowerPurity.push_back(999999);
        }
        //std::cout << "deb: iShowerPfp = " << fPFPShowerLength.back() << " matched MC position = " << pos << " matched MC PDG code = " << fPFPShowerTrueParticleMatchedPdgCode.back() << " nshared hits = " << fPFPShowerNSharedTrueParticleHits.back() << " completeness = " << fPFPShowerCompleteness.back() << " purity = " << fPFPShowerPurity.back() << std::endl;
        //std::cout << "iPfp = " << iPfp << " iShowerPfp = " << iShowerPfp << " g4ID = " << g4ID << " pos = " << pos << " pdg code = " << pfpShowerPdgCode << " shared hits = " << nSharedPfpShowerHits << " completeness = " << comp << " purity = " << purity << " length = " << shower->Length() << " opening angle = " << shower->OpenAngle() << std::endl;
      }
      else{
            fPFPShowerTrueParticleMatchedPosition.push_back(999999);
            fPFPShowerTrueParticleMatchedId.push_back(999999);
            fPFPShowerTrueParticleMatchedPdgCode.push_back(999999);
            fPFPShowerTrueParticleNHits.push_back(999999);
            fPFPShowerTrueParticleNHitsU.push_back(999999);
            fPFPShowerTrueParticleNHitsV.push_back(999999);
            fPFPShowerTrueParticleNHitsW.push_back(999999);
            fPFPShowerNSharedTrueParticleHits.push_back(999999);
            fPFPShowerCompleteness.push_back(999999);
            fPFPShowerPurity.push_back(999999);
      }
      iShowerPfp++;
    }
    ////////////////////////////////////////////////////////////////////
    iPfp++;
  }
  fNShowerPFParticles=iShowerPfp;
  std::cout << "DEBUG: fGlobalEventID = " << fGlobalEventID << " fPFPShowerTrueParticleMatchedPosition.size() = " << fPFPShowerTrueParticleMatchedPosition.size() << std::endl; 
  fEventTree->Fill();
  fGlobalEventID++;
}

void test::Pi0Analysis::beginJob()
{
  fGlobalEventID = 0;
  // Implementation of optional member function here.
  
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event","Event Tree");
  //fPi0Tree = tfs->make<TTree>("Pi0","Pi0 Tree");

  //Event tree
  fEventTree->Branch("globalEventID",&fGlobalEventID,"globalEventID/i");
  fEventTree->Branch("eventID",&fEventID,"eventID/i");
  fEventTree->Branch("runID",&fRunID,"runID/i");
  fEventTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fEventTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fEventTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  ////MC particle branches
  fEventTree->Branch("mcParticleTrackID", &fMCParticleTrackID);
  fEventTree->Branch("mcParticleMotherID", &fMCParticleMotherID);
  fEventTree->Branch("mcParticleMotherPosition", &fMCParticleMotherPosition);
  fEventTree->Branch("mcParticleMotherPdgCode", &fMCParticleMotherPdgCode);
  fEventTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode);
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
  fEventTree->Branch("showerPfpCompleteness",&fPFPShowerCompleteness);
  //fEventTree->Branch("pfpPurity",&fPFPPurity);
  fEventTree->Branch("showerPfpPurity",&fPFPShowerPurity);

  fEventTree->Branch("pfpShowerLength",&fPFPShowerLength);
  fEventTree->Branch("pfpShowerCollectionPlaneEnergy",&fPFPShowerCollectionPlaneEnergy);
  fEventTree->Branch("pfpShowerDirectionX",&fPFPShowerDirectionX);
  fEventTree->Branch("pfpShowerDirectionY",&fPFPShowerDirectionY);
  fEventTree->Branch("pfpShowerDirectionZ",&fPFPShowerDirectionZ);
  fEventTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpeningAngle);
  fEventTree->Branch("pfpShowerTrueParticleMatchedId",&fPFPShowerTrueParticleMatchedId);
  fEventTree->Branch("pfpShowerTrueParticleMatchedPdgCode",&fPFPShowerTrueParticleMatchedPdgCode);
  fEventTree->Branch("pfpShowerTrueParticleMatchedPosition",&fPFPShowerTrueParticleMatchedPosition);
  fEventTree->Branch("pfpShowerTrueParticleNHits",&fPFPShowerTrueParticleNHits);
  fEventTree->Branch("pfpShowerTrueParticleNHitsU",&fPFPShowerTrueParticleNHitsU);
  fEventTree->Branch("pfpShowerTrueParticleNHitsV",&fPFPShowerTrueParticleNHitsV);
  fEventTree->Branch("pfpShowerTrueParticleNHitsW",&fPFPShowerTrueParticleNHitsW);
  fEventTree->Branch("pfpShowerCompleteness",&fPFPShowerCompleteness);
  fEventTree->Branch("pfpShowerPurity",&fPFPShowerPurity);
  fEventTree->Branch("pfpShowerNHits",&fPFPShowerNHits);
  fEventTree->Branch("pfpShowerNHitsU",&fPFPShowerNHitsU);
  fEventTree->Branch("pfpShowerNHitsV",&fPFPShowerNHitsV);
  fEventTree->Branch("pfpShowerNHitsW",&fPFPShowerNHitsW);

  fEventTree->Branch("nMCPi0s",&fNMCPi0s,"nMCPi0s/i");
  fEventTree->Branch("nMCPhotons",&fNMCPhotons,"nMCPhotons/i");
  fEventTree->Branch("nReconstructedPhotons",&fNReconstructedPhotons,"nReconstructedPhotons/i");
  fEventTree->Branch("nReconstructedPi0s",&fNReconstructedPi0s,"nReconstructedPi0s/i");

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
        fMCParticleMotherID.push_back(trueParticle.Mother());
        //if(fMCParticleMotherID.back()==0)std::cout << "debug iMc = " << iMc << " track ID = " << fMCParticleTrackID.back() << " mother ID = " << fMCParticleMotherID.back() << std::endl;
        //if(fMCParticleMotherID.back()==0)std::cout << "mother position is = " << motherPosition << std::endl;
        fMCParticleNHits.push_back(fTrueParticleHits[trueParticle.TrackId()].size());
        int nHitsU(0), nHitsV(0), nHitsW(0);
        for(auto hit : fTrueParticleHits[trueParticle.TrackId()]){
	  if(hit->View()==geo::kU) nHitsU++; 
	  else if(hit->View()==geo::kV) nHitsV++; 
	  else if(hit->View()==geo::kW) nHitsW++; 
        }
        fMCParticleNHitsU.push_back(nHitsU);
        fMCParticleNHitsV.push_back(nHitsV);
        fMCParticleNHitsW.push_back(nHitsW);
        fMCParticlePdgCode.push_back(trueParticle.PdgCode());
        fMCParticleMass.push_back(trueParticle.Mass());
        fMCParticleEnergy.push_back(trueParticle.E());
        //std::cout << "trueParticle.PdgCode() = " << trueParticle.PdgCode() << " trueParticle.E() = " << trueParticle.E() << std::endl;
        if(trueParticle.PdgCode()==111) std::cout << "pi0! energy = " << trueParticle.E() << " global ev num = " << fGlobalEventID << " ev num = " << fEventID << " fRunID = " << fRunID << " fSubRunID = " << fSubRunID << std::endl;
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
        auto it = find(fMCParticleTrackID.begin(), fMCParticleTrackID.end(), trueParticle.Mother());
        //if(it != fMCParticleTrackID.end()) std::cout << "mother position is = " << it - fMCParticleTrackID.begin() << std::endl;
        //if (it == fMCParticleTrackID.end() && trueParticle.Mother()!=0) std::cout << "cannot find mother" << std::endl;
        int motherPosition(999999), motherPdgCode(999999);
        if(it != fMCParticleTrackID.end())  {
          motherPosition = it - fMCParticleTrackID.begin();
          motherPdgCode = fMCParticlePdgCode.at(motherPosition);
        }
        fMCParticleMotherPosition.push_back(motherPosition);
        fMCParticleMotherPdgCode.push_back(motherPdgCode);
        if(trueParticle.PdgCode()==22 && fMCParticleMotherPdgCode.back()==111) std::cout <<"photon from pi! energy = " << trueParticle.E() << std::endl;
        //std::cout << "pdg code = " << fMCParticlePdgCode.back() << " motherPosition = " << fMCParticleMotherPosition.back() << " mother pdg code = " << fMCParticleMotherPdgCode.back() << std::endl;
        //if(fMCParticleMotherPdgCode.back()==111)std::cout << "iMc = " << iMc << " pdg code = " << fMCParticlePdgCode.back() << " motherPosition = " << fMCParticleMotherPosition.back() << std::endl;
        //if(fMCParticleMotherPdgCode.back()==111) std::cout << "pdg: " << trueParticle.PdgCode() << " energy: " << trueParticle.E() << "mom: " << trueParticle.Px() << " " << trueParticle.Py() << " " << trueParticle.Pz() << std::endl;

     }
     //Fill mother position and pdg code info 
     //for(unsigned int iMc=0; iMc<fMcParticles->size(); iMc++){
     //    
     //}
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
      fMCParticleNHitsU.clear();
      fMCParticleNHitsV.clear();
      fMCParticleNHitsW.clear();
      fMCParticleTrackID.clear();
      fMCParticleMotherID.clear();
      fMCParticleMotherPosition.clear();
      fMCParticleMotherPdgCode.clear();
      fMCParticlePdgCode.clear();
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
      fPFPShowerCompleteness.clear();
      fPFPPurity.clear();
      fPFPShowerPurity.clear();
      fPFPShowerLength.clear();
      fPFPShowerDirectionX.clear();
      fPFPShowerDirectionY.clear();
      fPFPShowerDirectionZ.clear();
      fPFPShowerCollectionPlaneEnergy.clear();
      fPFPShowerOpeningAngle.clear();
      fPFPShowerNHits.clear();
      fPFPShowerNHitsU.clear();
      fPFPShowerNHitsV.clear();
      fPFPShowerNHitsW.clear();
      fPFPShowerTrueParticleMatchedId.clear();
      fPFPShowerTrueParticleMatchedPosition.clear();
      fPFPShowerTrueParticleMatchedPdgCode.clear();
      fPFPShowerTrueParticleNHits.clear();
      fPFPShowerTrueParticleNHitsU.clear();
      fPFPShowerTrueParticleNHitsV.clear();
      fPFPShowerTrueParticleNHitsW.clear();
      fPFPShowerNSharedTrueParticleHits.clear();
      //Pi0 vectors

}

int test::Pi0Analysis::GetVectorPositionForTrackID(int trackID, std::vector<int> trackIDvector){
      int pos(999999);
      for(unsigned int ipos=0; ipos<trackIDvector.size(); ipos++) {
        if(trackIDvector.at(ipos)==trackID){pos=ipos;} 
      }
      return pos;
}

DEFINE_ART_MODULE(test::Pi0Analysis)
