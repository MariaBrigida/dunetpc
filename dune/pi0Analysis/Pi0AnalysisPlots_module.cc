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

double xTolerance = 20;
double yTolerance = 20;
double zTolerance = 20;
double lowBoundaryX;
double lowBoundaryY;
double lowBoundaryZ;
double highBoundaryX;
double highBoundaryY;
double highBoundaryZ;

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

  bool PassesTruePhotonQualityCuts(int nTrueHitsU, int nTrueHitsV, int nTrueHitsW, double trueEnergy, double trueVertexX, double trueVertexY, double trueVertexZ);
  bool PassesTrueFiducialVolumeQualityCuts(double trueVertexX, double trueVertexY, double trueVertexZ);
  bool PassesRecoShowerQualityCuts(/*int nSharedHits, */double purity, double completeness);

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
  std::cout << "Start analysis" << std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  // Implementation of optional member function here.
  fTestTree = tfs->make<TTree>("Test","Test Tree");
  TH1D *hTruePi0InvMass, *hRecoPi0InvMass, *hRecoPi0InvMass_trueEnergy, *hRecoPi0InvMass_trueDirection, *hEnergyResolution, *hTruePhotonEnergy, *hTruePhotonsOpeningAngle, *hRecoPhotonsOpeningAngle, *hPhotonsOpeningAngleResolution;
  TH2I *hDebugTrackIDVsPfoMatchedID;
  TH2D *hTrueVsRecoPhotonEnergy, *hTrueVsRecoPhotonMomentumX, *hTrueVsRecoPhotonMomentumY, *hTrueVsRecoPhotonMomentumZ, *hShowerCompletenessVsEnergyResolution, *hShowerPurityVsEnergyResolution, *hTrueVsRecoPhotonsOpeningAngle;
  TH2D *hShowerNHitsVsEnergyResolution, *hShowerNHitsVsCompleteness;
  TH2D *hTrueEnergyVsEnergyResolution, *hTrueEnergyVsCompleteness;
  TH1I *hTruePhotonsPerPi0_trueQualityCuts, *hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts;
  //TH1I *hTruePi0s;
  TH1D *hEfficiency; 
  TH2D *hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts;
  TH2D *hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts;
  TH1D *hPi0TruePhotonsOpeningAngleEffi;
  TH2D *hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts;
  TH2D *hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts;
  TH1D *hPi0TruePhotonsMinimumEnergyEffi;

  hDebugTrackIDVsPfoMatchedID = tfs->make<TH2I>("hDebugTrackIDVsPfoMatchedID","hDebugTrackIDVsPfoMatchedID",1000,0,1000,1000,0,1000);
  hTruePi0InvMass = tfs->make<TH1D>("hTruePi0InvMass","hTruePi0InvMass",100,-0.05,0.45);
  hRecoPi0InvMass = tfs->make<TH1D>("hRecoPi0InvMass","hRecoPi0InvMass",100,-0.05,0.45);
  hRecoPi0InvMass_trueEnergy = tfs->make<TH1D>("hRecoPi0InvMass_trueEnergy","hRecoPi0InvMass_trueEnergy",100,-0.05,0.45);
  hRecoPi0InvMass_trueDirection = tfs->make<TH1D>("hRecoPi0InvMass_trueDirection","hRecoPi0InvMass_trueDirection",100,-0.05,0.45);
  hEnergyResolution = tfs->make<TH1D>("hEnergyResolution","hEnergyResolution",200,-4,4);
  hTrueVsRecoPhotonEnergy = tfs->make<TH2D>("hTrueVsRecoPhotonEnergy","hTrueVsRecoPhotonEnergy",100,0,4,100,0,4);
  hShowerCompletenessVsEnergyResolution = tfs->make<TH2D>("hShowerCompletenessVsEnergyResolution","hShowerCompletenessVsEnergyResolution",100,0,1,160,-2,4);
  hShowerPurityVsEnergyResolution = tfs->make<TH2D>("hShowerPurityVsEnergyResolution","hShowerPurityVsEnergyResolution",100,0,1,160,-2,4);
  hTruePhotonEnergy = tfs->make<TH1D>("hTruePhotonEnergy", "hTruePhotonEnergy", 100,0,10);
  hTrueEnergyVsEnergyResolution = tfs->make<TH2D>("hTrueEnergyVsEnergyResolution","hTrueEnergyVsEnergyResolution",100,0,20,160,-2,4);
  hTrueEnergyVsCompleteness = tfs->make<TH2D>("hTrueEnergyVsCompleteness","hTrueEnergyVsCompleteness",100,0,20,100,0,1);
  hShowerNHitsVsEnergyResolution = tfs->make<TH2D>("hShowerNHitsVsEnergyResolution","hShowerNHitsVsEnergyResolution",200,0,1000,160,-2,4);
  hShowerNHitsVsCompleteness = tfs->make<TH2D>("hShowerNHitsVsCompleteness","hShowerNHitsVsCompleteness",200,0,1000,100,0,1);
  //hMaxShowerCompletenessVsEnergyResolution = tfs->make<TH2D>("hMaxShowerCompletenessVsEnergyResolution","hMaxShowerCompletenessVsEnergyResolution",100,0,1,160,-2,4);
  //hMinShowerCompletenessVsEnergyResolution = tfs->make<TH2D>("hMinShowerCompletenessVsEnergyResolution","hMinShowerCompletenessVsEnergyResolution",100,0,1,160,-2,4);
  //hMeanShowerCompletenessVsEnergyResolution = tfs->make<TH2D>("hMeanShowerCompletenessVsEnergyResolution","hMeanShowerCompletenessVsEnergyResolution",100,0,1,160,-2,4);
  hTruePhotonsOpeningAngle = tfs->make<TH1D>("hTruePhotonsOpeningAngle","hTruePhotonsOpeningAngle",100,0,180);
  hRecoPhotonsOpeningAngle = tfs->make<TH1D>("hRecoPhotonsOpeningAngle","hRecoPhotonsOpeningAngle",100,0,180);
  hPhotonsOpeningAngleResolution = tfs->make<TH1D>("hPhotonsOpeningAngleResolution","hPhotonsOpeningAngleResolution",150,-1.5,6);
  hTrueVsRecoPhotonsOpeningAngle = tfs->make<TH2D>("hTrueVsRecoPhotonsOpeningAngle","hTrueVsRecoPhotonsOpeningAngle",100,0,180,100,0,180);
  hTrueVsRecoPhotonMomentumX = tfs->make<TH2D>("hTrueVsRecoPhotonMomentumX","hTrueVsRecoPhotonMomentumX",100,0,1,100,0,1);
  hTrueVsRecoPhotonMomentumY = tfs->make<TH2D>("hTrueVsRecoPhotonMomentumY","hTrueVsRecoPhotonMomentumY",100,0,1,100,0,1);
  hTrueVsRecoPhotonMomentumZ = tfs->make<TH2D>("hTrueVsRecoPhotonMomentumZ","hTrueVsRecoPhotonMomentumZ",100,0,1,100,0,1);
  hEfficiency = tfs->make<TH1D>("hEfficiency","hEfficiency",1,0.,1.);
  //hTruePi0s = tfs->make<TH1I>("hTruePi0s","hTruePi0s",1,0,1);
  hTruePhotonsPerPi0_trueQualityCuts = tfs->make<TH1I>("hTruePhotonsPerPi0_trueQualityCuts","hTruePhotonsPerPi0_trueQualityCuts",3,0,3);
  hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts = tfs->make<TH1I>("hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts","hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts",3,0,3);
  hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts = tfs->make<TH2D>("hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts","hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts",45,0,180,3,0,3);
  hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts = tfs->make<TH2D>("hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts","hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts",45,0,180,3,0,3);
  hPi0TruePhotonsOpeningAngleEffi = tfs->make<TH1D>("hPi0TruePhotonsOpeningAngleEffi","hPi0TruePhotonsOpeningAngleEffi",45,0,180);
  hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts = tfs->make<TH2D>("hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts","hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts",40,0,1,3,0,3);
  hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts = tfs->make<TH2D>("hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts","hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts",40,0,1,3,0,3);
  hPi0TruePhotonsMinimumEnergyEffi = tfs->make<TH1D>("hPi0TruePhotonsMinimumEnergyEffi","hPi0TruePhotonsMinimumEnergyEffi",40,0,1);

  TTree *eventTree = (TTree*)fTestFile->Get("ana/Event");
  TTree *globalTree = (TTree*)fTestFile->Get("ana/Global");
  //tree->SetBranchAddress("mcParticleVertexX",&mcParticleVertexX, &b_mcParticleVertexX);
 
  unsigned int globalEventID = 0; 
  unsigned int eventID = 0; 
  unsigned int runID = 0; 
  unsigned int subrunID = 0; 
  std::vector<int> *mcParticleTrackID = 0;
  std::vector<int> *mcParticlePdgCode = 0;
  std::vector<int> *mcParticleMotherPdgCode = 0;
  std::vector<int> *mcParticleMotherPosition = 0;
  std::vector<double> *mcParticleStartMomentumX = 0;
  std::vector<double> *mcParticleStartMomentumY = 0;
  std::vector<double> *mcParticleStartMomentumZ = 0;
  std::vector<double> *mcParticleEnergy = 0;
  std::vector<double> *mcParticleTrueVertexX = 0;
  std::vector<double> *mcParticleTrueVertexY = 0;
  std::vector<double> *mcParticleTrueVertexZ = 0;

  std::vector<double> *pfpShowerCollectionPlaneEnergy = 0;
  std::vector<double> *pfpShowerDirectionX = 0;
  std::vector<double> *pfpShowerDirectionY = 0;
  std::vector<double> *pfpShowerDirectionZ = 0;
  std::vector<double> *pfpShowerCompleteness = 0;
  std::vector<double> *pfpShowerPurity = 0;
  std::vector<int> *pfpShowerNHits = 0;
  std::vector<int> *pfpShowerNHitsU = 0;
  std::vector<int> *pfpShowerNHitsV = 0;
  std::vector<int> *pfpShowerNHitsW = 0;
  std::vector<int> *pfpShowerTrueParticleNHits = 0;
  std::vector<int> *pfpShowerTrueParticleNHitsU = 0;
  std::vector<int> *pfpShowerTrueParticleNHitsV = 0;
  std::vector<int> *pfpShowerTrueParticleNHitsW = 0;
  std::vector<int> *pfpShowerTrueParticleMatchedPosition = 0;
  std::vector<int> *pfpShowerTrueParticleMatchedId = 0;
 
  double lowBoundaryXCenter = 0;
  double highBoundaryXCenter = 0;
  double lowBoundaryYCenter = 0;
  double highBoundaryYCenter = 0;
  double lowBoundaryZCenter = 0;
  double highBoundaryZCenter = 0;
  double halfWidth = 0;
  double halfHeight = 0;
  double halfLength = 0;

  
  //Define Event tree branches
  TBranch *b_globalEventID, *b_eventID, *b_runID, *b_subrunID;
  TBranch *b_mcParticleTrackID, *b_mcParticlePdgCode, *b_mcParticleMotherPdgCode, *b_mcParticleMotherPosition;
  TBranch *b_mcParticleEnergy, *b_mcParticleTrueVertexX, *b_mcParticleTrueVertexY, *b_mcParticleTrueVertexZ;
  TBranch *b_mcParticleStartMomentumX, *b_mcParticleStartMomentumY, *b_mcParticleStartMomentumZ;
  TBranch *b_pfpShowerCollectionPlaneEnergy, *b_pfpShowerDirectionX, *b_pfpShowerDirectionY, *b_pfpShowerDirectionZ, *b_pfpShowerTrueParticleMatchedPosition, *b_pfpShowerTrueParticleMatchedId, *b_pfpShowerCompleteness, *b_pfpShowerPurity, *b_pfpShowerNHits, *b_pfpShowerNHitsU, *b_pfpShowerNHitsV, *b_pfpShowerNHitsW;
  TBranch *b_pfpShowerTrueParticleNHits, *b_pfpShowerTrueParticleNHitsU, *b_pfpShowerTrueParticleNHitsV, *b_pfpShowerTrueParticleNHitsW;
  TBranch *b_halfWidth, *b_halfHeight, *b_halfLength, *b_lowBoundaryXCenter, *b_lowBoundaryYCenter, *b_lowBoundaryZCenter, *b_highBoundaryXCenter, *b_highBoundaryYCenter, *b_highBoundaryZCenter;

  std::cout << "Read branches" << std::endl;
  //Set Event tree branch addresses
  eventTree->SetBranchAddress("eventID",&eventID, &b_eventID);
  eventTree->SetBranchAddress("runID",&runID, &b_runID);
  eventTree->SetBranchAddress("subrunID",&subrunID, &b_subrunID);
  eventTree->SetBranchAddress("globalEventID",&globalEventID, &b_globalEventID);
  eventTree->SetBranchAddress("mcParticleTrackID",&mcParticleTrackID, &b_mcParticleTrackID);
  eventTree->SetBranchAddress("mcParticlePdgCode",&mcParticlePdgCode, &b_mcParticlePdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPdgCode",&mcParticleMotherPdgCode, &b_mcParticleMotherPdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPosition",&mcParticleMotherPosition, &b_mcParticleMotherPosition);
  eventTree->SetBranchAddress("mcParticleStartMomentumX",&mcParticleStartMomentumX, &b_mcParticleStartMomentumX);
  eventTree->SetBranchAddress("mcParticleStartMomentumY",&mcParticleStartMomentumY, &b_mcParticleStartMomentumY);
  eventTree->SetBranchAddress("mcParticleStartMomentumZ",&mcParticleStartMomentumZ, &b_mcParticleStartMomentumZ);
  eventTree->SetBranchAddress("mcParticleEnergy",&mcParticleEnergy, &b_mcParticleEnergy);
  eventTree->SetBranchAddress("mcParticleVertexX",&mcParticleTrueVertexX, &b_mcParticleTrueVertexX);
  eventTree->SetBranchAddress("mcParticleVertexY",&mcParticleTrueVertexY, &b_mcParticleTrueVertexY);
  eventTree->SetBranchAddress("mcParticleVertexZ",&mcParticleTrueVertexZ, &b_mcParticleTrueVertexZ);

  eventTree->SetBranchAddress("pfpShowerCollectionPlaneEnergy",&pfpShowerCollectionPlaneEnergy, &b_pfpShowerCollectionPlaneEnergy);
  eventTree->SetBranchAddress("pfpShowerDirectionX",&pfpShowerDirectionX, &b_pfpShowerDirectionX);
  eventTree->SetBranchAddress("pfpShowerDirectionY",&pfpShowerDirectionY, &b_pfpShowerDirectionY);
  eventTree->SetBranchAddress("pfpShowerDirectionZ",&pfpShowerDirectionZ, &b_pfpShowerDirectionZ);
  eventTree->SetBranchAddress("pfpShowerTrueParticleMatchedPosition",&pfpShowerTrueParticleMatchedPosition, &b_pfpShowerTrueParticleMatchedPosition);
  eventTree->SetBranchAddress("pfpShowerTrueParticleMatchedId",&pfpShowerTrueParticleMatchedId, &b_pfpShowerTrueParticleMatchedId);
  eventTree->SetBranchAddress("pfpShowerCompleteness",&pfpShowerCompleteness, &b_pfpShowerCompleteness);
  eventTree->SetBranchAddress("pfpShowerPurity",&pfpShowerPurity, &b_pfpShowerPurity);
  eventTree->SetBranchAddress("pfpShowerNHits",&pfpShowerNHits, &b_pfpShowerNHits);
  eventTree->SetBranchAddress("pfpShowerNHitsU",&pfpShowerNHitsU, &b_pfpShowerNHitsU);
  eventTree->SetBranchAddress("pfpShowerNHitsV",&pfpShowerNHitsV, &b_pfpShowerNHitsV);
  eventTree->SetBranchAddress("pfpShowerNHitsW",&pfpShowerNHitsW, &b_pfpShowerNHitsW);
  eventTree->SetBranchAddress("pfpShowerTrueParticleNHits",&pfpShowerTrueParticleNHits, &b_pfpShowerTrueParticleNHits);
  eventTree->SetBranchAddress("pfpShowerTrueParticleNHitsU",&pfpShowerTrueParticleNHitsU, &b_pfpShowerTrueParticleNHitsU);
  eventTree->SetBranchAddress("pfpShowerTrueParticleNHitsV",&pfpShowerTrueParticleNHitsV, &b_pfpShowerTrueParticleNHitsV);
  eventTree->SetBranchAddress("pfpShowerTrueParticleNHitsW",&pfpShowerTrueParticleNHitsW, &b_pfpShowerTrueParticleNHitsW);

  //Geometry
  globalTree->SetBranchAddress("maxXCenter",&highBoundaryXCenter, &b_highBoundaryXCenter);
  globalTree->SetBranchAddress("minXCenter",&lowBoundaryXCenter, &b_lowBoundaryXCenter);
  globalTree->SetBranchAddress("maxYCenter",&highBoundaryYCenter, &b_highBoundaryYCenter);
  globalTree->SetBranchAddress("minYCenter",&lowBoundaryYCenter, &b_lowBoundaryYCenter);
  globalTree->SetBranchAddress("maxZCenter",&highBoundaryZCenter, &b_highBoundaryZCenter);
  globalTree->SetBranchAddress("minZCenter",&lowBoundaryZCenter, &b_lowBoundaryZCenter);
  globalTree->SetBranchAddress("halfWidth", &halfWidth, &b_halfWidth);
  globalTree->SetBranchAddress("halfHeight", &halfHeight, &b_halfHeight);
  globalTree->SetBranchAddress("halfLength", &halfLength, &b_halfLength);

  globalTree->GetEntry(0);
  lowBoundaryX = lowBoundaryXCenter - halfWidth + xTolerance;
  highBoundaryX = highBoundaryXCenter + halfWidth - xTolerance;
  lowBoundaryY = lowBoundaryYCenter - halfHeight + yTolerance;
  highBoundaryY = highBoundaryYCenter + halfHeight - yTolerance;
  lowBoundaryZ = lowBoundaryZCenter - halfLength + zTolerance;
  highBoundaryZ = highBoundaryZCenter + halfLength - zTolerance;

  /*std::cout << " lowBoundaryXCenter = " << lowBoundaryXCenter << " highBoundaryXCenter = " << highBoundaryXCenter << std::endl;
  std::cout << " lowBoundaryYCenter = " << lowBoundaryYCenter << " highBoundaryYCenter = " << highBoundaryYCenter << std::endl;
  std::cout << " lowBoundaryZCenter = " << lowBoundaryZCenter << " highBoundaryZCenter = " << highBoundaryZCenter << std::endl;

  std::cout << " halfWidth = " << halfWidth << " halfHeight = " << halfHeight << " halfLength = " << halfLength << std::endl;

  std::cout << " lowBoundaryX = " << lowBoundaryX << " highBoundaryX = " << highBoundaryX << std::endl;
  std::cout << " lowBoundaryY = " << lowBoundaryY << " highBoundaryY = " << highBoundaryY << std::endl;
  std::cout << " lowBoundaryZ = " << lowBoundaryZ << " highBoundaryZ = " << highBoundaryZ << std::endl;*/

  std::vector<int> pi0PhotonGlobalEventID, pi0PhotonMotherPosInMCParticleVect, pi0PhotonMatchedPfpShowerPos, pi0PhotonTrackId, pi0PhotonMatchedPfpShowerMatchedTrackID,pi0PhotonMotherPosInPi0Vect;
  std::vector<int> pi0PhotonShowerNHits, pi0PhotonShowerNHitsU, pi0PhotonShowerNHitsV, pi0PhotonShowerNHitsW;
  std::vector<int> pi0PhotonTrueParticleNHits, pi0PhotonTrueParticleNHitsU, pi0PhotonTrueParticleNHitsV, pi0PhotonTrueParticleNHitsW;
  std::vector<double> pi0PhotonTrueEnergy, pi0PhotonTrueDirectionX, pi0PhotonTrueDirectionY, pi0PhotonTrueDirectionZ, pi0PhotonTrueVertexX, pi0PhotonTrueVertexY, pi0PhotonTrueVertexZ;
  std::vector<double> pi0PhotonMatchedPfpShowerCollectionPlaneEnergy, pi0PhotonMatchedPfpShowerDirectionX, pi0PhotonMatchedPfpShowerDirectionY, pi0PhotonMatchedPfpShowerDirectionZ, pi0PhotonShowerCompleteness, pi0PhotonShowerPurity;


  //Loop over eventTree entries
  int nTruePi0s(0);
  int nTruePi0sInThisEvent(0);
  //std::vector<int> pi0TrackID;//vectors that save how many photons passing quality cuts are associated to each pi0
  std::vector<int> pi0GlobalEventNumber;
  std::vector<int> pi0PosInMCParticleVector;
  std::vector<int> pi0NPhotons_trueQualityCuts;;
  std::vector<int> pi0NPhotons_trueAndRecoQualityCuts;

  std::cout << "N Tree Entries = " << eventTree->GetEntries() << std::endl;
  //int nPhotons(0);
  for(int iEntry=0; iEntry<eventTree->GetEntries(); iEntry++){
    nTruePi0sInThisEvent=0;
    //pi0PosInMCParticleVector.clear();
    //pi0NPhotons_trueQualityCuts.clear();
    //pi0NPhotons_trueAndRecoQualityCuts.clear();
    eventTree->GetEntry(iEntry);
    //if(globalEventID!=314) continue;
    std::cout << "----------------------------------------> global ev number = " << globalEventID << "<----------------------------------------" << std::endl;
    //std::cout << "global event number = " << globalEventID << " event number = " <<eventID << " run number = " << runID << " subrun ID = " << subrunID << std::endl;
    //std::cout << " pfpShowerTrueParticleMatchedPosition->size() = " << pfpShowerTrueParticleMatchedPosition->size() << " mcParticleTrackID->size() = " << mcParticleTrackID->size() << std::endl;
    
    for(std::vector<int>::size_type iMCPart = 0; iMCPart<mcParticleTrackID->size(); iMCPart++){
      if(mcParticlePdgCode->at(iMCPart)==111) {
        //std::cout << "iMCPart = " << iMCPart << " globalEventID = " << globalEventID << std::endl;
        //hNTruePi0s->Fill(1);
        //pi0TrackID.push_back(mcParticleTrackID->at(iMCPart));
        pi0PosInMCParticleVector.push_back(iMCPart);
        pi0GlobalEventNumber.push_back(globalEventID);
        pi0NPhotons_trueQualityCuts.push_back(0);
        pi0NPhotons_trueAndRecoQualityCuts.push_back(0);
        nTruePi0s++;
	nTruePi0sInThisEvent++;
        //std::cout << "true pi0 found at event " << globalEventID << std::endl;
        //std::cout << " nr. pi0 = " << pi0NPhotons_trueQualityCuts.size() << std::endl;
      }
    }

    for(std::vector<int>::size_type iMCPart = 0; iMCPart<mcParticleTrackID->size(); iMCPart++){
      if(mcParticlePdgCode->at(iMCPart)==22 && mcParticleMotherPdgCode->at(iMCPart)==111){
          pi0PhotonMotherPosInMCParticleVect.push_back(mcParticleMotherPosition->at(iMCPart));
          pi0PhotonGlobalEventID.push_back(globalEventID);
          pi0PhotonTrackId.push_back(mcParticleTrackID->at(iMCPart));
          //std::cout << "debug iMCPart = " << iMCPart << " mcParticleTrackID->size = " << mcParticleTrackID->size() << std::endl;
          //std::cout << "debug photon from pi0. iMCPart = " << iMCPart << std::endl;
          pi0PhotonTrueEnergy.push_back(mcParticleEnergy->at(iMCPart));
          pi0PhotonTrueVertexX.push_back(mcParticleTrueVertexX->at(iMCPart));
          pi0PhotonTrueVertexY.push_back(mcParticleTrueVertexY->at(iMCPart));
          pi0PhotonTrueVertexZ.push_back(mcParticleTrueVertexZ->at(iMCPart));
          std::cout << "photon iMCPart = " << iMCPart << " true energy = " << mcParticleEnergy->at(iMCPart) << " track ID = " << mcParticleTrackID->at(iMCPart) << std::endl;
           //std::cout << "deb0" << std::endl;
          pi0PhotonTrueDirectionX.push_back(mcParticleStartMomentumX->at(iMCPart)); 
          // std::cout << "deb1" << std::endl;
          pi0PhotonTrueDirectionY.push_back(mcParticleStartMomentumY->at(iMCPart)); 
           //std::cout << "deb2" << std::endl;
          pi0PhotonTrueDirectionZ.push_back(mcParticleStartMomentumZ->at(iMCPart));
          // std::cout << "deb3" << std::endl;
          //std::cout << "pfpShowerTrueParticleMatchedPosition->size() = " << pfpShowerTrueParticleMatchedPosition->size() << std::endl;
          //std::cout << "pfpShowerTrueParticleMatchedId->size() = " << pfpShowerTrueParticleMatchedId->size() << std::endl;
          int nMatchedRecoShowers(0);
          for(std::vector<int>::size_type iPfp=0; iPfp<pfpShowerTrueParticleMatchedId->size(); iPfp++){
            if(pfpShowerCollectionPlaneEnergy->at(iPfp)<0) continue; //Filter out Pfps where energy value is not good!
            //std::cout << " nMatchedRecoShowers = " << nMatchedRecoShowers << " pfpShowerTrueParticleMatchedId->at(iPfp) = " << pfpShowerTrueParticleMatchedId->at(iPfp) << " mcParticleTrackID->at(iMCPart) = " << mcParticleTrackID->at(iMCPart) << std::endl;
            if(nMatchedRecoShowers==0 && (pfpShowerTrueParticleMatchedId->at(iPfp)==mcParticleTrackID->at(iMCPart))) {
              //std::cout << "And it matched to pfp at position = " << iPfp << std::endl;
              nMatchedRecoShowers++;
              pi0PhotonMatchedPfpShowerPos.push_back(iPfp);
              pi0PhotonMatchedPfpShowerMatchedTrackID.push_back(pfpShowerTrueParticleMatchedId->at(iPfp)); 
              pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.push_back(pfpShowerCollectionPlaneEnergy->at(iPfp));
              pi0PhotonMatchedPfpShowerDirectionX.push_back(pfpShowerDirectionX->at(iPfp));
              pi0PhotonMatchedPfpShowerDirectionY.push_back(pfpShowerDirectionY->at(iPfp));
              pi0PhotonMatchedPfpShowerDirectionZ.push_back(pfpShowerDirectionZ->at(iPfp));
              pi0PhotonShowerCompleteness.push_back(pfpShowerCompleteness->at(iPfp));
              pi0PhotonShowerPurity.push_back(pfpShowerPurity->at(iPfp));
              pi0PhotonShowerNHits.push_back(pfpShowerNHits->at(iPfp));
              pi0PhotonShowerNHitsU.push_back(pfpShowerNHitsU->at(iPfp));
              pi0PhotonShowerNHitsV.push_back(pfpShowerNHitsV->at(iPfp));
              pi0PhotonShowerNHitsW.push_back(pfpShowerNHitsW->at(iPfp));
              pi0PhotonTrueParticleNHits.push_back(pfpShowerTrueParticleNHits->at(iPfp));
              pi0PhotonTrueParticleNHitsU.push_back(pfpShowerTrueParticleNHitsU->at(iPfp));
              pi0PhotonTrueParticleNHitsV.push_back(pfpShowerTrueParticleNHitsV->at(iPfp));
              pi0PhotonTrueParticleNHitsW.push_back(pfpShowerTrueParticleNHitsW->at(iPfp));
              //std::cout <<"iPfp = " << iPfp << " completeness = " << pfpShowerCompleteness->at(iPfp) << " reco energy = " << pfpShowerCollectionPlaneEnergy->at(iPfp) << " matched trackID = " << pfpShowerTrueParticleMatchedId->at(iPfp) << std::endl;
            }
            //If I already found a pfp matching to this MC particle, only replace if purity*completeness is higher
            else if(nMatchedRecoShowers==1 && (pfpShowerTrueParticleMatchedId->at(iPfp)==mcParticleTrackID->at(iMCPart))){
              //std::cout << "found a better match. replacing..." << std::endl;
              double purity = pfpShowerPurity->at(iPfp);
              double completeness = pfpShowerCompleteness->at(iPfp);
              if(completeness*purity>pi0PhotonShowerCompleteness.back()*pi0PhotonShowerPurity.back()){
                pi0PhotonMatchedPfpShowerPos.back() = iPfp;
                pi0PhotonMatchedPfpShowerMatchedTrackID.back() = pfpShowerTrueParticleMatchedId->at(iPfp); 
                pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.back() = pfpShowerCollectionPlaneEnergy->at(iPfp);
                pi0PhotonMatchedPfpShowerDirectionX.back() = pfpShowerDirectionX->at(iPfp);
                pi0PhotonMatchedPfpShowerDirectionY.back() = pfpShowerDirectionY->at(iPfp);
                pi0PhotonMatchedPfpShowerDirectionZ.back() = pfpShowerDirectionZ->at(iPfp);
                pi0PhotonShowerCompleteness.back() = pfpShowerCompleteness->at(iPfp);
                pi0PhotonShowerPurity.back() = pfpShowerPurity->at(iPfp);
                pi0PhotonShowerNHits.back()= pfpShowerNHits->at(iPfp);
                pi0PhotonShowerNHitsU.back()= pfpShowerNHitsU->at(iPfp);
                pi0PhotonShowerNHitsW.back()= pfpShowerNHitsV->at(iPfp);
                pi0PhotonShowerNHitsV.back()= pfpShowerNHitsW->at(iPfp);
                pi0PhotonTrueParticleNHits.back() = pfpShowerTrueParticleNHits->at(iPfp);
                pi0PhotonTrueParticleNHitsU.back() = pfpShowerTrueParticleNHitsU->at(iPfp);
                pi0PhotonTrueParticleNHitsV.back() = pfpShowerTrueParticleNHitsV->at(iPfp);
                pi0PhotonTrueParticleNHitsW.back() = pfpShowerTrueParticleNHitsW->at(iPfp);
                //std::cout <<"REPLACING iPfp = " << iPfp << " completeness = " << pfpShowerCompleteness->at(iPfp) << " reco energy = " << pfpShowerCollectionPlaneEnergy->at(iPfp) << " matched trackID = " << pfpShowerTrueParticleMatchedId->at(iPfp) << std::endl;
              }
            }
          }
          if(!nMatchedRecoShowers){
              //std::cout << "And it did not match to any pfp" << std::endl;
              pi0PhotonMatchedPfpShowerPos.push_back(-999);
              pi0PhotonMatchedPfpShowerMatchedTrackID.push_back(-999); 
              pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.push_back(-999);
              pi0PhotonMatchedPfpShowerDirectionX.push_back(-999);
              pi0PhotonMatchedPfpShowerDirectionY.push_back(-999);
              pi0PhotonMatchedPfpShowerDirectionZ.push_back(-999);
              pi0PhotonShowerCompleteness.push_back(-999);
              pi0PhotonShowerPurity.push_back(-999);
              pi0PhotonShowerNHits.push_back(-999);
              pi0PhotonShowerNHitsU.push_back(-999);
              pi0PhotonShowerNHitsV.push_back(-999);
              pi0PhotonShowerNHitsW.push_back(-999);
              pi0PhotonTrueParticleNHits.push_back(-999);
              pi0PhotonTrueParticleNHitsU.push_back(-999);
              pi0PhotonTrueParticleNHitsV.push_back(-999);
              pi0PhotonTrueParticleNHitsW.push_back(-999);
          }

          //Find pi0 mother position
          int pi0ID(-999);
          for(std::vector<int>::size_type iMotherPos=0; iMotherPos<pi0PosInMCParticleVector.size(); iMotherPos++){
            if(pi0PosInMCParticleVector.at(iMotherPos)==pi0PhotonMotherPosInMCParticleVect.back() && pi0GlobalEventNumber.at(iMotherPos)==pi0PhotonGlobalEventID.back()) {pi0ID=iMotherPos;std::cout << "globalEventID = " << globalEventID << " photonID  = "<<  pi0PhotonTrueEnergy.size() << " energu = " << pi0PhotonTrueEnergy.back() << " mother pos in pi0 vect= " << pi0ID << std::endl;} 
          }
	  pi0PhotonMotherPosInPi0Vect.push_back(pi0ID);
	  //std::cout << "photon of energy = "<< pi0PhotonTrueEnergy.back() << "  at event = " << pi0PhotonGlobalEventID.back() << "passes MC cuts? " << PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.back(),pi0PhotonTrueParticleNHitsV.back(),pi0PhotonTrueParticleNHitsW.back(), pi0PhotonTrueEnergy.back()) << std::endl;
          if(pi0ID>=0 && PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.back(),pi0PhotonTrueParticleNHitsV.back(),pi0PhotonTrueParticleNHitsW.back(), pi0PhotonTrueEnergy.back(),pi0PhotonTrueVertexX.back(),pi0PhotonTrueVertexY.back(),pi0PhotonTrueVertexZ.back())){pi0NPhotons_trueQualityCuts.at(pi0ID)++;std::cout << "globalEventID = " << globalEventID << " photonID passing true cuts = "<<  pi0PhotonTrueEnergy.size() << " energu = " << pi0PhotonTrueEnergy.back() << " mother pos in pi0 vect= " << pi0ID << std::endl;}
          if(pi0ID>=0 && PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.back(),pi0PhotonTrueParticleNHitsV.back(),pi0PhotonTrueParticleNHitsW.back(), pi0PhotonTrueEnergy.back(),pi0PhotonTrueVertexX.back(),pi0PhotonTrueVertexY.back(),pi0PhotonTrueVertexZ.back()) && PassesRecoShowerQualityCuts(pi0PhotonShowerPurity.back(),pi0PhotonShowerCompleteness.back())) {pi0NPhotons_trueAndRecoQualityCuts.at(pi0ID)++; std::cout << "globalEventID = " << globalEventID << " photonID passing all cuts = "<<  pi0PhotonTrueEnergy.size() << " energu = " << pi0PhotonTrueEnergy.back() << " mother pos in pi0 vect= " << pi0ID << std::endl;}
        }
      }

      //IDENTIFY EVENTS FOR WHICH a specific pi0 has 0,1, or 2 reconstructed photons, for debugging purposes (making ev displays)
      int nPhotonsPerPi0(0),nPhotonsPerPi0PassingTrueQualityCuts(0),nPhotonsPerPi0PassingTrueAndRecoQualityCuts(0), photon1IDPerPi0(-999), photon2IDPerPi0(-999);
      std::cout << "nTruePi0sInThisEvent = " << nTruePi0sInThisEvent << std::endl;
      for(long unsigned int iPi0=(pi0PosInMCParticleVector.size()-nTruePi0sInThisEvent); iPi0<pi0PosInMCParticleVector.size(); iPi0++){
        nPhotonsPerPi0=0;nPhotonsPerPi0PassingTrueQualityCuts=0;nPhotonsPerPi0PassingTrueAndRecoQualityCuts=0;
        //Let's find the daughter photons associated to this pi0
        for(long unsigned int iPi0Phot=0; iPi0Phot<pi0PhotonMotherPosInPi0Vect.size(); iPi0Phot++){
          if(pi0GlobalEventNumber.at(iPi0)==pi0PhotonGlobalEventID.at(iPi0Phot))std::cout << " photonID = " << iPi0Phot << " energy = " << pi0PhotonTrueEnergy.at(iPi0Phot) << " mother pos in pi0 vect = " << pi0PhotonMotherPosInPi0Vect.at(iPi0Phot) << " pi0PosInMCParticleVector in pi0 vect = " << pi0PosInMCParticleVector.at(iPi0) << " pi0GlobalEventNumber.at(iPi0) = " << pi0GlobalEventNumber.at(iPi0) << " pi0PhotonGlobalEventID.at(iPi0Phot) = " << pi0PhotonGlobalEventID.at(iPi0Phot) << std::endl;
	  if(pi0PhotonMotherPosInPi0Vect.at(iPi0Phot)<0)continue;
          if((int)iPi0==pi0PhotonMotherPosInPi0Vect.at(iPi0Phot) && pi0GlobalEventNumber.at(iPi0)==pi0PhotonGlobalEventID.at(iPi0Phot)){
            nPhotonsPerPi0++;
            if(photon1IDPerPi0<0)photon1IDPerPi0=iPi0Phot;
	    else photon2IDPerPi0=iPi0Phot;
            if(PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(iPi0Phot),pi0PhotonTrueParticleNHitsV.at(iPi0Phot),pi0PhotonTrueParticleNHitsW.at(iPi0Phot), pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonTrueVertexX.at(iPi0Phot),pi0PhotonTrueVertexY.at(iPi0Phot),pi0PhotonTrueVertexZ.at(iPi0Phot)))nPhotonsPerPi0PassingTrueQualityCuts++;
            if(PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(iPi0Phot),pi0PhotonTrueParticleNHitsV.at(iPi0Phot),pi0PhotonTrueParticleNHitsW.at(iPi0Phot), pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonTrueVertexX.at(iPi0Phot),pi0PhotonTrueVertexY.at(iPi0Phot),pi0PhotonTrueVertexZ.at(iPi0Phot)) && PassesRecoShowerQualityCuts(pi0PhotonShowerPurity.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot)))nPhotonsPerPi0PassingTrueAndRecoQualityCuts++;
            //std::cout << "iPi0 = " << iPi0 << " iPi0Phit = " << iPi0Phot << " true photon energy = " << pi0PhotonTrueEnergy.at(iPi0Phot) << std::endl;
          }
        }
        /*if(nPhotonsPerPi0==0)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " NOT RECONSTRUCTED" << std::endl; 
        else if(nPhotonsPerPi0==1)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " ONLY 1 RECONSTRUCTED PHOTON" << std::endl;
        else if(nPhotonsPerPi0==2)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " BOTH PHOTONS RECONSTRUCTED!!!" << std::endl;
        if(nPhotonsPerPi0PassingTrueQualityCuts==0)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " HAS NO PHOTON PASSING TRUE QUALITY CUTS" << std::endl; 
        if(nPhotonsPerPi0PassingTrueQualityCuts==1)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " HAS ONLY 1 PHOTON PASSING TRUE QUALITY CUTS" << std::endl; 
        if(nPhotonsPerPi0PassingTrueQualityCuts==2)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " BOTH PHOTONS PASS TRUE QUALITY CUTS!!!" << std::endl; 
        if(nPhotonsPerPi0PassingTrueAndRecoQualityCuts==0)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " HAS NO PHOTON PASSING TRUE AND RECO QUALITY CUTS" << std::endl; 
        if(nPhotonsPerPi0PassingTrueAndRecoQualityCuts==1)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " HAS ONLY 1 PHOTON PASSING TRUE AND RECO QUALITY CUTS" << std::endl; 
        if(nPhotonsPerPi0PassingTrueAndRecoQualityCuts==2)std::cout << "DIAGNOSTICS: pi0 at pos = " << iPi0 << " in entry = " << iEntry << " BOTH PHOTONS PASS TRUE AND RECO QUALITY CUTS!!!" << std::endl;
        std::cout << "DEBUG photon1IDPerPi0 = " << photon1IDPerPi0 << " photon2IDPerPi0 = " << photon2IDPerPi0 << std::endl;*/
	if(photon1IDPerPi0>0 && photon2IDPerPi0>0){ //I think the only case in which this won't work is the Dalitz case
          TVector3 photon1trueDir(pi0PhotonTrueDirectionX.at(photon1IDPerPi0),pi0PhotonTrueDirectionY.at(photon1IDPerPi0),pi0PhotonTrueDirectionZ.at(photon1IDPerPi0));
	  TVector3 photon2trueDir(pi0PhotonTrueDirectionX.at(photon2IDPerPi0),pi0PhotonTrueDirectionY.at(photon2IDPerPi0),pi0PhotonTrueDirectionZ.at(photon2IDPerPi0));
          double pi0PhotonsTrueAngle=photon1trueDir.Angle(photon2trueDir)*180/TMath::Pi();
          double minimumPhotonEnergy=TMath::Min(pi0PhotonTrueEnergy.at(photon1IDPerPi0),pi0PhotonTrueEnergy.at(photon2IDPerPi0));
          if(minimumPhotonEnergy>0.6) std::cout << "high energy photons!! event number = " << pi0GlobalEventNumber.at(iPi0) << std::endl; 
          hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->Fill(pi0PhotonsTrueAngle,nPhotonsPerPi0PassingTrueQualityCuts);   
          hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts->Fill(pi0PhotonsTrueAngle,nPhotonsPerPi0PassingTrueAndRecoQualityCuts);  
	  hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts->Fill(minimumPhotonEnergy,nPhotonsPerPi0PassingTrueQualityCuts);
	  hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts->Fill(minimumPhotonEnergy,nPhotonsPerPi0PassingTrueAndRecoQualityCuts);
        }
      }
    }

  //MAYBE THE BLOCK BELOW IS NOT REALLY NEEDED, AND CAN DO ALL EFFICIENCY IN THE ABOVE PI0 LOOPS! ACTUALLY MAYBE EVERYTHING THAT FOLLOWS CAN BE DONE ABOVE!!! AND MAYBE IT CAN ALL BE DONE IN THE VERY FIRST LOOP?
  //std::cout << "pi0PhotonMotherPos size = " << pi0PhotonMotherPosInMCParticleVect.size() << std::endl;
  int nPi0_noPhot_trueQualityCuts= std::count_if(pi0NPhotons_trueQualityCuts.begin(), pi0NPhotons_trueQualityCuts.end(), [](int i){return i == 0;}); 
  int nPi0_onePhot_trueQualityCuts= std::count_if(pi0NPhotons_trueQualityCuts.begin(), pi0NPhotons_trueQualityCuts.end(), [](int i){return i == 1;}); 
  int nPi0_twoPhot_trueQualityCuts= std::count_if(pi0NPhotons_trueQualityCuts.begin(), pi0NPhotons_trueQualityCuts.end(), [](int i){return i == 2;}); 
  int nPi0_noPhot_trueAndRecoQualityCuts = std::count_if(pi0NPhotons_trueAndRecoQualityCuts.begin(), pi0NPhotons_trueAndRecoQualityCuts.end(), [](int i){return i == 0;}); 
  int nPi0_onePhot_trueAndRecoQualityCuts = std::count_if(pi0NPhotons_trueAndRecoQualityCuts.begin(), pi0NPhotons_trueAndRecoQualityCuts.end(), [](int i){return i == 1;}); 
  int nPi0_twoPhot_trueAndRecoQualityCuts = std::count_if(pi0NPhotons_trueAndRecoQualityCuts.begin(), pi0NPhotons_trueAndRecoQualityCuts.end(), [](int i){return i == 2;}); 

  std::cout << "N true pi0s = " << pi0NPhotons_trueQualityCuts.size() << " nr pi0s with zero/one/two photons passing true quality cuts = " << nPi0_noPhot_trueQualityCuts<< " / " << nPi0_onePhot_trueQualityCuts<< " / " << nPi0_twoPhot_trueQualityCuts<< std::endl;
  std::cout << "N true pi0s = " << pi0NPhotons_trueQualityCuts.size() << " nr pi0s with zero/one/two photons passing true and reco quality cuts = " << nPi0_noPhot_trueAndRecoQualityCuts << " / " << nPi0_onePhot_trueAndRecoQualityCuts << " / " << nPi0_twoPhot_trueAndRecoQualityCuts << std::endl;
  //hTruePi0s->SetBinContent(0,pi0NPhotons_trueQualityCuts.size());
  hTruePhotonsPerPi0_trueQualityCuts->SetBinContent(1,nPi0_noPhot_trueQualityCuts);
  hTruePhotonsPerPi0_trueQualityCuts->SetBinContent(2,nPi0_onePhot_trueQualityCuts);
  hTruePhotonsPerPi0_trueQualityCuts->SetBinContent(3,nPi0_twoPhot_trueQualityCuts);
  hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts->SetBinContent(1,nPi0_noPhot_trueAndRecoQualityCuts);
  hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts->SetBinContent(2,nPi0_onePhot_trueAndRecoQualityCuts);
  hReconstructedPhotonsPerPi0_trueAndRecoQualityCuts->SetBinContent(3,nPi0_twoPhot_trueAndRecoQualityCuts);
  //double effi = (double)nPi0_twoPhot_trueAndRecoQualityCuts/pi0NPhotons_trueQualityCuts.size();
  double effi = (double)nPi0_twoPhot_trueAndRecoQualityCuts/(double)nPi0_twoPhot_trueQualityCuts;
  double effiError = TMath::Sqrt(effi*(1-effi)/(double)nPi0_twoPhot_trueQualityCuts);
  std::cout << "effi = " << effi << " +- " << effiError << std::endl;
  //std::cout << " pi0NPhotons_trueQualityCuts.size() = " << pi0NPhotons_trueQualityCuts.size() << " nPi0_noPhot_trueQualityCuts= " << nPi0_noPhot_trueQualityCuts<< " nPi0_onePhot_trueQualityCuts= " << nPi0_onePhot_trueQualityCuts<< " nPi0_twoPhot_trueQualityCuts= " << nPi0_twoPhot_trueQualityCuts<< " efficiency = " << effi << std::endl;
  hEfficiency->SetBinContent(1,effi);
  //Make pi0 mass plots

  //Pi0 vectors
  std::vector<double> pi0Daughter1TrueEnergy, pi0Daughter1TrueDirectionX, pi0Daughter1TrueDirectionY, pi0Daughter1TrueDirectionZ;
  std::vector<double> pi0Daughter2TrueEnergy, pi0Daughter2TrueDirectionX, pi0Daughter2TrueDirectionY, pi0Daughter2TrueDirectionZ;
  
  int nPhotonPairs(0); //for debug
   //std::cout << " globalEventID = " << globalEventID << " DEBUG pi0PhotonMotherPosInMCParticleVect.size() = " << pi0PhotonMotherPosInMCParticleVect.size() << std::endl;

  for(long unsigned int iPi0Phot=0; iPi0Phot<pi0PhotonMotherPosInMCParticleVect.size(); iPi0Phot++){
     //if(pi0PhotonTrueParticleNHitsU.at(iPi0Phot)<0) continue;

     //std::cout << "globalEventNumber = " << pi0PhotonGlobalEventID.at(iPi0Phot) << std::endl;
    //std::cout << "photon's mother position = " << pi0PhotonMotherPosInMCParticleVect.at(iPi0Phot) << std::endl;
    for(long unsigned int jPi0Phot=iPi0Phot+1; jPi0Phot<pi0PhotonMotherPosInMCParticleVect.size(); jPi0Phot++){

     //if(pi0PhotonTrueParticleNHitsU.at(jPi0Phot)<0) continue;
     if(pi0PhotonMotherPosInMCParticleVect.at(iPi0Phot)==pi0PhotonMotherPosInMCParticleVect.at(jPi0Phot) && pi0PhotonGlobalEventID.at(iPi0Phot)==pi0PhotonGlobalEventID.at(jPi0Phot) && pi0PhotonMatchedPfpShowerPos.at(iPi0Phot)>=0 && pi0PhotonMatchedPfpShowerPos.at(jPi0Phot)>=0){
     

     double recoPhot1Energy = pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.at(iPi0Phot)/1000;
     double energyResValuePhot1 = (recoPhot1Energy-pi0PhotonTrueEnergy.at(iPi0Phot))/pi0PhotonTrueEnergy.at(iPi0Phot);
     hTruePhotonEnergy->Fill(pi0PhotonTrueEnergy.at(iPi0Phot));
     hTrueEnergyVsEnergyResolution->Fill(pi0PhotonTrueEnergy.at(iPi0Phot),energyResValuePhot1);
     hTrueEnergyVsCompleteness->Fill(pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot));
     hShowerNHitsVsEnergyResolution->Fill(pi0PhotonShowerNHits.at(iPi0Phot), energyResValuePhot1);
     hShowerNHitsVsCompleteness->Fill(pi0PhotonShowerNHits.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot));
     hShowerCompletenessVsEnergyResolution->Fill(pi0PhotonShowerCompleteness.at(iPi0Phot),energyResValuePhot1);
     hShowerPurityVsEnergyResolution->Fill(pi0PhotonShowerPurity.at(iPi0Phot),energyResValuePhot1);
     TVector3 phot1mom(pi0PhotonTrueDirectionX.at(iPi0Phot),pi0PhotonTrueDirectionY.at(iPi0Phot),pi0PhotonTrueDirectionZ.at(iPi0Phot));
     TVector3 recoPhot1Mom(pi0PhotonMatchedPfpShowerDirectionX.at(iPi0Phot),pi0PhotonMatchedPfpShowerDirectionY.at(iPi0Phot),pi0PhotonMatchedPfpShowerDirectionZ.at(iPi0Phot));

     double recoPhot2Energy = pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.at(jPi0Phot)/1000;
     double energyResValuePhot2 = (recoPhot2Energy-pi0PhotonTrueEnergy.at(jPi0Phot))/pi0PhotonTrueEnergy.at(jPi0Phot);
     hTruePhotonEnergy->Fill(pi0PhotonTrueEnergy.at(jPi0Phot));
     hTrueEnergyVsEnergyResolution->Fill(pi0PhotonTrueEnergy.at(jPi0Phot),energyResValuePhot2);
     hTrueEnergyVsCompleteness->Fill(pi0PhotonTrueEnergy.at(jPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot));
     hShowerNHitsVsEnergyResolution->Fill(pi0PhotonShowerNHits.at(jPi0Phot), energyResValuePhot2);
     hShowerNHitsVsCompleteness->Fill(pi0PhotonShowerNHits.at(jPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot));
     hShowerCompletenessVsEnergyResolution->Fill(pi0PhotonShowerCompleteness.at(jPi0Phot),energyResValuePhot2);
     hShowerPurityVsEnergyResolution->Fill(pi0PhotonShowerPurity.at(jPi0Phot),energyResValuePhot2);
     hEnergyResolution->Fill(energyResValuePhot1);
     hEnergyResolution->Fill(energyResValuePhot2);
     TVector3 phot2mom(pi0PhotonTrueDirectionX.at(jPi0Phot),pi0PhotonTrueDirectionY.at(jPi0Phot),pi0PhotonTrueDirectionZ.at(jPi0Phot));
     TVector3 recoPhot2Mom(pi0PhotonMatchedPfpShowerDirectionX.at(jPi0Phot),pi0PhotonMatchedPfpShowerDirectionY.at(jPi0Phot),pi0PhotonMatchedPfpShowerDirectionZ.at(jPi0Phot));

     if(!PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(iPi0Phot),pi0PhotonTrueParticleNHitsV.at(iPi0Phot),pi0PhotonTrueParticleNHitsW.at(iPi0Phot), pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonTrueVertexX.at(iPi0Phot),pi0PhotonTrueVertexY.at(iPi0Phot),pi0PhotonTrueVertexZ.at(iPi0Phot)) || !PassesRecoShowerQualityCuts(pi0PhotonShowerPurity.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot))) continue;
     if(!PassesTruePhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(jPi0Phot),pi0PhotonTrueParticleNHitsV.at(jPi0Phot),pi0PhotonTrueParticleNHitsW.at(jPi0Phot),pi0PhotonTrueEnergy.at(jPi0Phot),pi0PhotonTrueVertexX.at(jPi0Phot),pi0PhotonTrueVertexY.at(jPi0Phot),pi0PhotonTrueVertexZ.at(jPi0Phot)) || !PassesRecoShowerQualityCuts(pi0PhotonShowerPurity.at(jPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot))) continue;
        nPhotonPairs++;

        double angle = phot1mom.Angle(phot2mom);
        double angleDeg = angle/TMath::Pi()*180;
        double invMass = TMath::Sqrt(2*pi0PhotonTrueEnergy.at(iPi0Phot)*pi0PhotonTrueEnergy.at(jPi0Phot)*(1-TMath::Cos(angle))); 
           //std::cout << "debC" << std::endl;
           //std::cout << "debD" << std::endl;
           //std::cout << "debE" << std::endl;
        if(pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.at(iPi0Phot)<0 || pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.at(jPi0Phot)<0) continue;
           //std::cout << "debF" << std::endl;
        double recoAngle = recoPhot1Mom.Angle(recoPhot2Mom);
        double recoAngleDeg = recoAngle/TMath::Pi()*180;
        double recoInvMass = TMath::Sqrt(2*recoPhot1Energy*recoPhot2Energy*(1-TMath::Cos(recoAngle)));
        double recoInvMassTrueEnergy = TMath::Sqrt(2*pi0PhotonTrueEnergy.at(iPi0Phot)*pi0PhotonTrueEnergy.at(jPi0Phot)*(1-TMath::Cos(recoAngle)));
        double recoInvMassTrueDirection = TMath::Sqrt(2*recoPhot1Energy*recoPhot2Energy*(1-TMath::Cos(angle))); 
        //double minCompleteness = TMath::Min(pi0PhotonShowerCompleteness.at(iPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot));
        //double maxCompleteness = TMath::Max(pi0PhotonShowerCompleteness.at(iPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot));
        //double meanCompleteness = (pi0PhotonShowerCompleteness.at(iPi0Phot)+pi0PhotonShowerCompleteness.at(jPi0Phot))/2;
        //std::cout << "true phot 1 trackID E mom: " << pi0PhotonTrackId.at(iPi0Phot) << " " << pi0PhotonTrueEnergy.at(iPi0Phot) << " " << pi0PhotonTrueDirectionX.at(iPi0Phot) << " " << pi0PhotonTrueDirectionY.at(iPi0Phot) << " " << pi0PhotonTrueDirectionZ.at(iPi0Phot) << std::endl;
        //std::cout << "true phot 2 trackID E mom: " << pi0PhotonTrackId.at(jPi0Phot) << " " << pi0PhotonTrueEnergy.at(jPi0Phot) << " " << pi0PhotonTrueDirectionX.at(jPi0Phot) << " " << pi0PhotonTrueDirectionY.at(jPi0Phot) << " " << pi0PhotonTrueDirectionZ.at(jPi0Phot) << std::endl;
        //std::cout << "true angle = " << angle << " " << " true inv mass = " << invMass <<std::endl;
        //std::cout << "reco phot 1 matchedID E mom: " << pi0PhotonMatchedPfpShowerMatchedTrackID.at(iPi0Phot) << " " << recoPhot1Energy << " " << pi0PhotonMatchedPfpShowerDirectionX.at(iPi0Phot) << " " << pi0PhotonMatchedPfpShowerDirectionY.at(iPi0Phot) << " " << pi0PhotonMatchedPfpShowerDirectionZ.at(iPi0Phot) << " completeness = " << pi0PhotonShowerCompleteness.at(iPi0Phot) << " purity = " << pi0PhotonShowerPurity.at(iPi0Phot) <<std::endl;
        //std::cout << "reco phot 2 matchedID E mom: " << pi0PhotonMatchedPfpShowerMatchedTrackID.at(jPi0Phot) << " " << recoPhot2Energy << " " << pi0PhotonMatchedPfpShowerDirectionX.at(jPi0Phot) << " " << pi0PhotonMatchedPfpShowerDirectionY.at(jPi0Phot) << " " << pi0PhotonMatchedPfpShowerDirectionZ.at(jPi0Phot) <<" completeness = " << pi0PhotonShowerCompleteness.at(jPi0Phot) << " purity = " << pi0PhotonShowerPurity.at(jPi0Phot) << std::endl;
        //std::cout << "global ev ID = " << pi0PhotonGlobalEventID.at(iPi0Phot) << " " << pi0PhotonGlobalEventID.at(jPi0Phot) << " reco angle = " << recoAngle << " reco inv mass = " << recoInvMass << std::endl;
        //Fill plots
        //if(recoInvMass>0.2) std::cout << "inv Mass = " << recoInvMass << " globalEventID = " << pi0PhotonGlobalEventID.at(iPi0Phot) << std::endl;
        //std::cout << "inv Mass = " << recoInvMass << " globalEventID = " << pi0PhotonGlobalEventID.at(iPi0Phot) << std::endl;

        hTruePi0InvMass->Fill(invMass);
        hRecoPi0InvMass->Fill(recoInvMass);
        hRecoPi0InvMass_trueEnergy->Fill(recoInvMassTrueEnergy);
        hRecoPi0InvMass_trueDirection->Fill(recoInvMassTrueDirection);
	hDebugTrackIDVsPfoMatchedID->Fill(pi0PhotonTrackId.at(iPi0Phot),pi0PhotonMatchedPfpShowerMatchedTrackID.at(iPi0Phot));
	hDebugTrackIDVsPfoMatchedID->Fill(pi0PhotonTrackId.at(jPi0Phot),pi0PhotonMatchedPfpShowerMatchedTrackID.at(jPi0Phot));
        hTrueVsRecoPhotonEnergy->Fill(pi0PhotonTrueEnergy.at(iPi0Phot),recoPhot1Energy);
        hTrueVsRecoPhotonEnergy->Fill(pi0PhotonTrueEnergy.at(jPi0Phot),recoPhot2Energy);
        hTrueVsRecoPhotonMomentumX->Fill(phot1mom.X()/phot1mom.Mag(),recoPhot1Mom.X()/recoPhot1Mom.Mag());
        hTrueVsRecoPhotonMomentumX->Fill(phot2mom.X()/phot2mom.Mag(),recoPhot2Mom.X()/recoPhot2Mom.Mag());
        hTrueVsRecoPhotonMomentumY->Fill(phot1mom.Y()/phot1mom.Mag(),recoPhot1Mom.Y()/recoPhot1Mom.Mag());
        hTrueVsRecoPhotonMomentumY->Fill(phot2mom.Y()/phot2mom.Mag(),recoPhot2Mom.Y()/recoPhot2Mom.Mag());
        hTrueVsRecoPhotonMomentumZ->Fill(phot1mom.Z()/phot1mom.Mag(),recoPhot1Mom.Z()/recoPhot1Mom.Mag());
        hTrueVsRecoPhotonMomentumZ->Fill(phot2mom.Z()/phot2mom.Mag(),recoPhot2Mom.Z()/recoPhot2Mom.Mag());
        hTruePhotonsOpeningAngle->Fill(angleDeg);
        hRecoPhotonsOpeningAngle->Fill(recoAngleDeg);
        hPhotonsOpeningAngleResolution->Fill((recoAngleDeg-angleDeg)/angleDeg);
        hTrueVsRecoPhotonsOpeningAngle->Fill(angleDeg,recoAngleDeg);

        //Differentiate between photon1 and photon2 and other reco energy planes, add quality cuts
      }
    }
  }
  //hNTruePi0s_twoRecoPhotons->Fill(1);
  //std::cout << "nPhotonPairs = " << nPhotonPairs << std::endl;
  //
  //
  //MAKE PI0 EFFICIENCY PLOTS
  int nBinsAngle=hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->GetNbinsX();
  int nBinsEnergy=hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts->GetNbinsX();
  double angleEffi(0), energyEffi(0), angleEffiError(0), energyEffiError(0);
  for(int iBin=0; iBin<nBinsAngle; iBin++){
    if(hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsAngle,3,3)){
      angleEffi=hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts->Integral(iBin,nBinsAngle,3,3)/hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsAngle,3,3);
      angleEffiError=TMath::Sqrt(angleEffi*(1-angleEffi)/hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsAngle,3,3));
    }
    else {angleEffi=0;angleEffiError=0;}
    //std::cout << " angle effi for bin = " << iBin << " = " << angleEffi << " num = " << hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueAndRecoQualityCuts->Integral(iBin,nBinsAngle,3,3) << " den = " << hPi0TruePhotonsOpeningAngleVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsAngle,3,3) << std::endl;
    hPi0TruePhotonsOpeningAngleEffi->SetBinContent(iBin,angleEffi);
    hPi0TruePhotonsOpeningAngleEffi->SetBinError(iBin,angleEffiError);
  }
  for(int iBin=0; iBin<nBinsEnergy; iBin++){
    if(hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsEnergy,3,3)){
      energyEffi=hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueAndRecoQualityCuts->Integral(iBin,nBinsEnergy,3,3)/hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsEnergy,3,3);
      energyEffiError=TMath::Sqrt(energyEffi*(1-energyEffi)/hPi0TruePhotonsMinimumEnergyVsNTruePhotonsPassingTrueQualityCuts->Integral(iBin,nBinsEnergy,3,3));
    }
    else {energyEffi=0;energyEffiError=0;}
    //std::cout << " energy effi for bin = " << iBin << " = " << energyEffi << std::endl;
    hPi0TruePhotonsMinimumEnergyEffi->SetBinContent(iBin,energyEffi);
    hPi0TruePhotonsMinimumEnergyEffi->SetBinError(iBin,energyEffiError);
  } 

}

void test::Pi0AnalysisPlots::endJob()
{
  // Implementation of optional member function here.
}

bool test::Pi0AnalysisPlots::PassesTruePhotonQualityCuts(int nTrueHitsU, int nTrueHitsV, int nTrueHitsW, double trueEnergy, double trueVertexX, double trueVertexY, double trueVertexZ)
{
  if(nTrueHitsU+nTrueHitsV+nTrueHitsW<15) return false;
  if((nTrueHitsU<5 && nTrueHitsV<5) || (nTrueHitsU<5 && nTrueHitsW<5) || (nTrueHitsV<5 && nTrueHitsW<5)) return false;
  if(!PassesTrueFiducialVolumeQualityCuts(trueVertexX,trueVertexY,trueVertexZ)) return false;
  //if(trueEnergy<0.1) return false; //100 MeV cut
  return true;
}

bool test::Pi0AnalysisPlots::PassesRecoShowerQualityCuts(double purity, double completeness)
{
  if(purity<0.5) return false;
  if(completeness<0.5) return false;
  return true;
}

bool test::Pi0AnalysisPlots::PassesTrueFiducialVolumeQualityCuts(double trueVertexX, double trueVertexY, double trueVertexZ)
{
  if(trueVertexX>highBoundaryX-xTolerance || trueVertexX<lowBoundaryX+xTolerance) return false;
  if(trueVertexY>highBoundaryY-yTolerance || trueVertexY<lowBoundaryY+yTolerance) return false;
  if(trueVertexZ>highBoundaryZ-zTolerance || trueVertexZ<lowBoundaryZ+zTolerance) return false;
  return true;
}

DEFINE_ART_MODULE(test::Pi0AnalysisPlots)
