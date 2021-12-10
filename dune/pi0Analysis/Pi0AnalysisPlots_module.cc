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

  bool PassesPhotonQualityCuts(int nHitsU, int nHitsV, int nHitsW, double energy, double purity, double completeness);

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
  TH1D *hTruePi0InvMass, *hRecoPi0InvMass, *hRecoPi0InvMass_trueEnergy, *hRecoPi0InvMass_trueDirection, *hEnergyResolution, *hTruePhotonEnergy, *hTruePhotonsOpeningAngle, *hRecoPhotonsOpeningAngle, *hPhotonsOpeningAngleResolution;
  TH2I *hDebugTrackIDVsPfoMatchedID;
  TH2D *hTrueVsRecoPhotonEnergy, *hTrueVsRecoPhotonMomentumX, *hTrueVsRecoPhotonMomentumY, *hTrueVsRecoPhotonMomentumZ, *hShowerCompletenessVsEnergyResolution, *hShowerPurityVsEnergyResolution, *hTrueVsRecoPhotonsOpeningAngle;
  TH2D *hShowerNHitsVsEnergyResolution, *hShowerNHitsVsCompleteness;
  TH2D *hTrueEnergyVsEnergyResolution, *hTrueEnergyVsCompleteness;
  TH1I *hReconstructedPhotonsPerPi0, *hReconstructedPhotonsPerPi0_qualityCuts;
  TH1D *hEfficiency; 

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
  hReconstructedPhotonsPerPi0 = tfs->make<TH1I>("hReconstructedPhotonsPerPi0","hReconstructedPhotonsPerPi0",3,0,3);
  hReconstructedPhotonsPerPi0_qualityCuts = tfs->make<TH1I>("hReconstructedPhotonsPerPi0_qualityCuts","hReconstructedPhotonsPerPi0_qualityCuts",3,0,3);

  TTree *eventTree = (TTree*)fTestFile->Get("ana/Event");
  //eventTree->Print();
  //std::vector<double> *mcParticleVertexX = 0;
  //TBranch *b_mcParticleVertexX;
  //tree->SetBranchAddress("mcParticleVertexX",&mcParticleVertexX, &b_mcParticleVertexX);
 
  unsigned int globalEventID = 0; 
  std::vector<int> *mcParticleTrackID = 0;
  std::vector<int> *mcParticlePdgCode = 0;
  std::vector<int> *mcParticleMotherPdgCode = 0;
  std::vector<int> *mcParticleMotherPosition = 0;
  std::vector<double> *mcParticleStartMomentumX = 0;
  std::vector<double> *mcParticleStartMomentumY = 0;
  std::vector<double> *mcParticleStartMomentumZ = 0;
  std::vector<double> *mcParticleEnergy = 0;

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

  
  //Define Event tree branches
  TBranch *b_globalEventID;
  TBranch *b_mcParticleTrackID, *b_mcParticlePdgCode, *b_mcParticleMotherPdgCode, *b_mcParticleMotherPosition;
  TBranch *b_mcParticleEnergy;
  TBranch *b_mcParticleStartMomentumX, *b_mcParticleStartMomentumY, *b_mcParticleStartMomentumZ;
  TBranch *b_pfpShowerCollectionPlaneEnergy, *b_pfpShowerDirectionX, *b_pfpShowerDirectionY, *b_pfpShowerDirectionZ, *b_pfpShowerTrueParticleMatchedPosition, *b_pfpShowerTrueParticleMatchedId, *b_pfpShowerCompleteness, *b_pfpShowerPurity, *b_pfpShowerNHits, *b_pfpShowerNHitsU, *b_pfpShowerNHitsV, *b_pfpShowerNHitsW;
  TBranch *b_pfpShowerTrueParticleNHits, *b_pfpShowerTrueParticleNHitsU, *b_pfpShowerTrueParticleNHitsV, *b_pfpShowerTrueParticleNHitsW;

  //Set Event tree branch addresses
  eventTree->SetBranchAddress("globalEventID",&globalEventID, &b_globalEventID);
  eventTree->SetBranchAddress("mcParticleTrackID",&mcParticleTrackID, &b_mcParticleTrackID);
  eventTree->SetBranchAddress("mcParticlePdgCode",&mcParticlePdgCode, &b_mcParticlePdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPdgCode",&mcParticleMotherPdgCode, &b_mcParticleMotherPdgCode);
  eventTree->SetBranchAddress("mcParticleMotherPosition",&mcParticleMotherPosition, &b_mcParticleMotherPosition);
  eventTree->SetBranchAddress("mcParticleStartMomentumX",&mcParticleStartMomentumX, &b_mcParticleStartMomentumX);
  eventTree->SetBranchAddress("mcParticleStartMomentumY",&mcParticleStartMomentumY, &b_mcParticleStartMomentumY);
  eventTree->SetBranchAddress("mcParticleStartMomentumZ",&mcParticleStartMomentumZ, &b_mcParticleStartMomentumZ);
  eventTree->SetBranchAddress("mcParticleEnergy",&mcParticleEnergy, &b_mcParticleEnergy);

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

  std::vector<int> pi0PhotonGlobalEventID, pi0PhotonMotherPos, pi0PhotonMatchedPfpShowerPos, pi0PhotonTrackId, pi0PhotonMatchedPfpShowerMatchedTrackID;
  std::vector<int> pi0PhotonShowerNHits, pi0PhotonShowerNHitsU, pi0PhotonShowerNHitsV, pi0PhotonShowerNHitsW;
  std::vector<int> pi0PhotonTrueParticleNHits, pi0PhotonTrueParticleNHitsU, pi0PhotonTrueParticleNHitsV, pi0PhotonTrueParticleNHitsW;
  std::vector<double> pi0PhotonTrueEnergy, pi0PhotonTrueDirectionX, pi0PhotonTrueDirectionY, pi0PhotonTrueDirectionZ;
  std::vector<double> pi0PhotonMatchedPfpShowerCollectionPlaneEnergy, pi0PhotonMatchedPfpShowerDirectionX, pi0PhotonMatchedPfpShowerDirectionY, pi0PhotonMatchedPfpShowerDirectionZ, pi0PhotonShowerCompleteness, pi0PhotonShowerPurity;


  //Loop over eventTree entries
  int nTruePi0s(0);
  //std::vector<int> pi0TrackID;//vectors that save how many photons passing quality cuts are associated to each pi0
  std::vector<int> pi0GlobalEventNumber;
  std::vector<int> pi0Pos;
  std::vector<int> pi0NPhotons;
  std::vector<int> pi0NPhotons_qualityCuts;

  std::cout << "N Tree Entries = " << eventTree->GetEntries() << std::endl;
  for(int iEntry=0; iEntry<eventTree->GetEntries(); iEntry++){
    //pi0Pos.clear();
    //pi0NPhotons.clear();
    //pi0NPhotons_qualityCuts.clear();
    eventTree->GetEntry(iEntry);
    //if(globalEventID>20000) break;
    //std::cout << "global ev number = " << globalEventID << std::endl;
    
    for(std::vector<int>::size_type iMCPart = 0; iMCPart<mcParticleTrackID->size(); iMCPart++){
      if(mcParticlePdgCode->at(iMCPart)==111) {
        //std::cout << "iMCPart = " << iMCPart << " globalEventID = " << globalEventID << std::endl;
        //hNTruePi0s->Fill(1);
        //pi0TrackID.push_back(mcParticleTrackID->at(iMCPart));
        pi0Pos.push_back(iMCPart);
        pi0GlobalEventNumber.push_back(globalEventID);
        pi0NPhotons.push_back(0);
        pi0NPhotons_qualityCuts.push_back(0);
        nTruePi0s++;
        //std::cout << " nr. pi0 = " << pi0NPhotons.size() << std::endl;
      }
    }

    for(std::vector<int>::size_type iMCPart = 0; iMCPart<mcParticleTrackID->size(); iMCPart++){
      if(mcParticlePdgCode->at(iMCPart)==22 && mcParticleMotherPdgCode->at(iMCPart)==111){
          pi0PhotonMotherPos.push_back(mcParticleMotherPosition->at(iMCPart));
          pi0PhotonGlobalEventID.push_back(globalEventID);
          pi0PhotonTrackId.push_back(mcParticleTrackID->at(iMCPart));
          //std::cout << "debug iMCPart = " << iMCPart << " mcParticleTrackID->size = " << mcParticleTrackID->size() << std::endl;
          pi0PhotonTrueEnergy.push_back(mcParticleEnergy->at(iMCPart));
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
          for(std::vector<int>::size_type iPfp=0; iPfp<pfpShowerTrueParticleMatchedPosition->size(); iPfp++){
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
            }
            //If I already found a pfp matching to this MC particle, only replace if purity*completeness is higher
            else if(nMatchedRecoShowers==1 && (pfpShowerTrueParticleMatchedId->at(iPfp)==mcParticleTrackID->at(iMCPart))){
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

              }
            else continue;
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
              pi0PhotonTrueParticleNHits.push_back(-999);
              pi0PhotonTrueParticleNHitsU.push_back(-999);
              pi0PhotonTrueParticleNHitsV.push_back(-999);
              pi0PhotonTrueParticleNHitsW.push_back(-999);
          }

          //std::cout << " photon index = " << pi0PhotonMotherPos.size() -1 << " pi0PhotonGlobalEventID = " << (int)pi0PhotonGlobalEventID.back() << " pi0PhotonTrueParticleNHitsU = " << pi0PhotonTrueParticleNHitsU.back() << " pi0PhotonMotherPos.size() = " << pi0PhotonMotherPos.size() << " pi0PhotonGlobalEventID.size() = " << (int)pi0PhotonGlobalEventID.size() << " pi0PhotonTrueParticleNHits.size() = " << pi0PhotonTrueParticleNHits.size() << std::endl;
          //Fill pi0 Efficiency calculation vectors
          /*if(((pi0PhotonMotherPos.size()-1)==118) || ((pi0PhotonMotherPos.size()-1)==119)){
		std::cout << " I am looking for this mother position: " << pi0PhotonMotherPos.back() << " in the pi0 vectors. I will loop over them." << std::endl;
          	for(std::vector<int>::size_type iPi0Vect=0; iPi0Vect<pi0Pos.size(); iPi0Vect++){
            	std::cout << "iPi0Vect = " << iPi0Vect << " pi0Pos.at(iPi0Vect) = " << pi0Pos.at(iPi0Vect) << " pi0GlobalEventNumber.at(iPi0Vect) = " << pi0GlobalEventNumber.at(iPi0Vect) << std::endl;
          	}
          }*/

          //Find pi0 mother position
          int pi0ID(-999);
          for(std::vector<int>::size_type iMotherPos=0; iMotherPos<pi0Pos.size(); iMotherPos++){
            if(pi0Pos.at(iMotherPos)==pi0PhotonMotherPos.back() && pi0GlobalEventNumber.at(iMotherPos)==pi0PhotonGlobalEventID.back()) pi0ID=iMotherPos; 
          }

          //std::vector<int>::iterator it = std::find(pi0Pos.begin(), pi0Pos.end(), pi0PhotonMotherPos.back());
          //int pi0ID = std::distance(pi0Pos.begin(), it);

          //if(((pi0PhotonMotherPos.size()-1)==118) || ((pi0PhotonMotherPos.size()-1)==119))std::cout << " photon index = "  << pi0PhotonMotherPos.size()-1 << " pi0PhotonTrueParticleNHitsU " << pi0PhotonTrueParticleNHitsU.at(pi0PhotonMotherPos.size()-1) << " pi0PhotonTrueParticleNHitsV " << pi0PhotonTrueParticleNHitsV.at(pi0PhotonMotherPos.size()-1) << " pi0PhotonTrueParticleNHitsW " << pi0PhotonTrueParticleNHitsW.at(pi0PhotonMotherPos.size()-1) << " pi0PhotonTrueEnergy " << pi0PhotonTrueEnergy.at(pi0PhotonMotherPos.size()-1) << " pi0PhotonShowerPurity " << pi0PhotonShowerPurity.at(pi0PhotonMotherPos.size()-1) << " pi0PhotonShowerCompleteness " << pi0PhotonShowerCompleteness.at(pi0PhotonMotherPos.size()-1) << " Passes element cuts = " <<PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.back(),pi0PhotonTrueParticleNHitsV.back(),pi0PhotonTrueParticleNHitsW.back(), pi0PhotonTrueEnergy.back(),pi0PhotonShowerPurity.back(),pi0PhotonShowerCompleteness.back())<< " (int)pi0GlobalEventNumber.at(pi0ID) = " << (int)pi0GlobalEventNumber.at(pi0ID) << " (int)globalEventID = " << (int)globalEventID << std::endl;
          //std::cout << "Pi0 daughter photon! iMCPart = " << iMCPart << " iPhoton = " << pi0PhotonTrueParticleNHitsU.size() - 1 << " pi0ID = " << pi0ID << " pi0GlobalEventNumber = " << (int)pi0GlobalEventNumber.at(pi0ID) << std::endl;
          if(pi0ID>=0)pi0NPhotons.at(pi0ID)++;
          if(pi0ID>=0 && PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.back(),pi0PhotonTrueParticleNHitsV.back(),pi0PhotonTrueParticleNHitsW.back(), pi0PhotonTrueEnergy.back(),pi0PhotonShowerPurity.back(),pi0PhotonShowerCompleteness.back())) {pi0NPhotons_qualityCuts.at(pi0ID)++;}
            //std::cout << "iMCPart = "  << iMCPart << " iPhoton = " << pi0PhotonTrueParticleNHitsU.size() - 1 << " pi0PhotonTrueParticleNHitsU " << pi0PhotonTrueParticleNHitsU.back() << " pi0PhotonTrueParticleNHitsV " << pi0PhotonTrueParticleNHitsV.back() << " pi0PhotonTrueParticleNHitsW " << pi0PhotonTrueParticleNHitsW.back() << " pi0PhotonTrueEnergy " << pi0PhotonTrueEnergy.back() << " pi0PhotonShowerPurity " << pi0PhotonShowerPurity.back() << " pi0PhotonShowerCompleteness " << pi0PhotonShowerCompleteness.back() << std::endl;
          //std::cout << "iMc = " << iMCPart << " pi0PhotonMotherPos= " << pi0PhotonMotherPos.back() << std::endl;
          //std::cout << "debug: pi0PhotonMatchedPfpShowerMatchedTrackID = " << pi0PhotonMatchedPfpShowerMatchedTrackID.back() << " pi0PhotonTrackId = " << pi0PhotonTrackId.back() << std::endl;

        }
      }
      //int nPhot(0); for(long unsigned int iPi0=0; iPi0<pi0NPhotons.size(); iPi0++){nPhot+=pi0NPhotons.at(iPi0);}
      //for(long unsigned int iPi0=0; iPi0<pi0NPhotons.size(); iPi0++){
        //std::cout << "debug: iPi0 = " << iPi0 << " pi0NPhotons = " << pi0NPhotons.at(iPi0) << std::endl;
      //}
      //std::cout << " nr photons = " << nPhot << " size of photon array = " << pi0PhotonMatchedPfpShowerPos.size() << std::endl;
      
    }
  //std::cout << "pi0PhotonMotherPos size = " << pi0PhotonMotherPos.size() << std::endl;
  int nPi0_noPhot = std::count_if(pi0NPhotons.begin(), pi0NPhotons.end(), [](int i){return i == 0;}); 
  int nPi0_onePhot = std::count_if(pi0NPhotons.begin(), pi0NPhotons.end(), [](int i){return i == 1;}); 
  int nPi0_twoPhot = std::count_if(pi0NPhotons.begin(), pi0NPhotons.end(), [](int i){return i == 2;}); 
  int nPi0_noPhot_qualityCuts = std::count_if(pi0NPhotons_qualityCuts.begin(), pi0NPhotons_qualityCuts.end(), [](int i){return i == 0;}); 
  int nPi0_onePhot_qualityCuts = std::count_if(pi0NPhotons_qualityCuts.begin(), pi0NPhotons_qualityCuts.end(), [](int i){return i == 1;}); 
  int nPi0_twoPhot_qualityCuts = std::count_if(pi0NPhotons_qualityCuts.begin(), pi0NPhotons_qualityCuts.end(), [](int i){return i == 2;}); 

  std::cout << "N true pi0s = " << pi0NPhotons.size() << " nr pi0s with zero/one/two reco photon = " << nPi0_noPhot << " / " << nPi0_onePhot << " / " << nPi0_twoPhot << std::endl;
  std::cout << "N true pi0s = " << pi0NPhotons.size() << " nr pi0s with zero/one/two reco photons with quality cuts = " << nPi0_noPhot_qualityCuts << " / " << nPi0_onePhot_qualityCuts << " / " << nPi0_twoPhot_qualityCuts << std::endl;
  hReconstructedPhotonsPerPi0->SetBinContent(1,nPi0_noPhot);
  hReconstructedPhotonsPerPi0->SetBinContent(2,nPi0_onePhot);
  hReconstructedPhotonsPerPi0->SetBinContent(3,nPi0_twoPhot);
  hReconstructedPhotonsPerPi0_qualityCuts->SetBinContent(1,nPi0_noPhot_qualityCuts);
  hReconstructedPhotonsPerPi0_qualityCuts->SetBinContent(2,nPi0_onePhot_qualityCuts);
  hReconstructedPhotonsPerPi0_qualityCuts->SetBinContent(3,nPi0_twoPhot_qualityCuts);
  double effi = (double)nPi0_twoPhot_qualityCuts/pi0NPhotons.size();
  std::cout << " pi0NPhotons.size() = " << pi0NPhotons.size() << " nPi0_noPhot = " << nPi0_noPhot << " nPi0_onePhot = " << nPi0_onePhot << " nPi0_twoPhot = " << nPi0_twoPhot << " efficiency = " << effi << std::endl;
  hEfficiency->SetBinContent(1,effi);
  //Make pi0 mass plots

  //Pi0 vectors
  std::vector<double> pi0Daughter1TrueEnergy, pi0Daughter1TrueDirectionX, pi0Daughter1TrueDirectionY, pi0Daughter1TrueDirectionZ;
  std::vector<double> pi0Daughter2TrueEnergy, pi0Daughter2TrueDirectionX, pi0Daughter2TrueDirectionY, pi0Daughter2TrueDirectionZ;
  
  int nPhotonPairs(0); //for debug
  for(long unsigned int iPi0Phot=0; iPi0Phot<pi0PhotonMotherPos.size(); iPi0Phot++){
     if(pi0PhotonTrueParticleNHitsU.at(iPi0Phot)<0) continue;

     //std::cout << "globalEventNumber = " << pi0PhotonGlobalEventID.at(iPi0Phot) << std::endl;
    //std::cout << "photon's mother position = " << pi0PhotonMotherPos.at(iPi0Phot) << std::endl;
    for(long unsigned int jPi0Phot=iPi0Phot+1; jPi0Phot<pi0PhotonMotherPos.size(); jPi0Phot++){
     if(pi0PhotonTrueParticleNHitsU.at(jPi0Phot)<0) continue;
     if(pi0PhotonMotherPos.at(iPi0Phot)==pi0PhotonMotherPos.at(jPi0Phot) && pi0PhotonGlobalEventID.at(iPi0Phot)==pi0PhotonGlobalEventID.at(jPi0Phot) && pi0PhotonMatchedPfpShowerPos.at(iPi0Phot)>=0 && pi0PhotonMatchedPfpShowerPos.at(jPi0Phot)>=0){

     double recoPhot1Energy = pi0PhotonMatchedPfpShowerCollectionPlaneEnergy.at(iPi0Phot)/1000;
     double energyResValuePhot1 = (recoPhot1Energy-pi0PhotonTrueEnergy.at(iPi0Phot))/pi0PhotonTrueEnergy.at(iPi0Phot);
     hTruePhotonEnergy->Fill(pi0PhotonTrueEnergy.at(iPi0Phot));
     hTrueEnergyVsEnergyResolution->Fill(pi0PhotonTrueEnergy.at(iPi0Phot),energyResValuePhot1);
     hTrueEnergyVsCompleteness->Fill(pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot));
     hShowerNHitsVsEnergyResolution->Fill(pi0PhotonShowerNHits.at(iPi0Phot), energyResValuePhot1);
     hShowerNHitsVsCompleteness->Fill(pi0PhotonShowerNHits.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot));
     hShowerCompletenessVsEnergyResolution->Fill(pi0PhotonShowerCompleteness.at(iPi0Phot),energyResValuePhot1);
     hShowerPurityVsEnergyResolution->Fill(pi0PhotonShowerPurity.at(iPi0Phot),energyResValuePhot1);


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

     if(!PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(iPi0Phot),pi0PhotonTrueParticleNHitsV.at(iPi0Phot),pi0PhotonTrueParticleNHitsW.at(iPi0Phot), pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonShowerPurity.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot))) continue;
      if(!PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(jPi0Phot),pi0PhotonTrueParticleNHitsV.at(jPi0Phot),pi0PhotonTrueParticleNHitsW.at(jPi0Phot),pi0PhotonTrueEnergy.at(jPi0Phot),pi0PhotonShowerPurity.at(jPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot))) continue;
     //if(iPi0Phot==118)std::cout << "iPi0Phot = "  << iPi0Phot << " pi0PhotonTrueParticleNHitsU " << pi0PhotonTrueParticleNHitsU.at(iPi0Phot) << " pi0PhotonTrueParticleNHitsV " << pi0PhotonTrueParticleNHitsV.at(iPi0Phot) << " pi0PhotonTrueParticleNHitsW " << pi0PhotonTrueParticleNHitsW.at(iPi0Phot) << " pi0PhotonTrueEnergy " << pi0PhotonTrueEnergy.at(iPi0Phot) << " pi0PhotonShowerPurity " << pi0PhotonShowerPurity.at(iPi0Phot) << " pi0PhotonShowerCompleteness " << pi0PhotonShowerCompleteness.at(iPi0Phot) << " passes element cuts = " << PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(iPi0Phot),pi0PhotonTrueParticleNHitsV.at(iPi0Phot),pi0PhotonTrueParticleNHitsW.at(iPi0Phot), pi0PhotonTrueEnergy.at(iPi0Phot),pi0PhotonShowerPurity.at(iPi0Phot),pi0PhotonShowerCompleteness.at(iPi0Phot)) << std::endl;
    // if(jPi0Phot==119)std::cout << "jPi0Phot = "  << jPi0Phot << " pi0PhotonTrueParticleNHitsU " << pi0PhotonTrueParticleNHitsU.at(jPi0Phot) << " pi0PhotonTrueParticleNHitsV " << pi0PhotonTrueParticleNHitsV.at(jPi0Phot) << " pi0PhotonTrueParticleNHitsW " << pi0PhotonTrueParticleNHitsW.at(jPi0Phot) << " pi0PhotonTrueEnergy " << pi0PhotonTrueEnergy.at(jPi0Phot) << " pi0PhotonShowerPurity " << pi0PhotonShowerPurity.at(jPi0Phot) << " pi0PhotonShowerCompleteness " << pi0PhotonShowerCompleteness.at(jPi0Phot) << " passes element cuts = " << PassesPhotonQualityCuts(pi0PhotonTrueParticleNHitsU.at(jPi0Phot),pi0PhotonTrueParticleNHitsV.at(jPi0Phot),pi0PhotonTrueParticleNHitsW.at(jPi0Phot),pi0PhotonTrueEnergy.at(jPi0Phot),pi0PhotonShowerPurity.at(jPi0Phot),pi0PhotonShowerCompleteness.at(jPi0Phot)) << std::endl;
        nPhotonPairs++;
        //hNTruePi0s_twoRecoPhotons->Fill(1);
        //std::cout << "PASSED CUTS! i= " << iPi0Phot << " j= " << jPi0Phot << " pi0PhotonGlobalEventID = " << pi0PhotonGlobalEventID.at(jPi0Phot) << " pi0PhotonTrueParticleNHitsU.at(iPi0Phot) = " << pi0PhotonTrueParticleNHitsU.at(iPi0Phot) << " pi0PhotonTrueParticleNHitsU.at(jPi0Phot) = " << pi0PhotonTrueParticleNHitsU.at(jPi0Phot) << std::endl;
        //std::cout << "i= " << iPi0Phot << " j= " << jPi0Phot << " pi position = " << pi0PhotonMotherPos.at(iPi0Phot) << std::endl;
        TVector3 phot1mom(pi0PhotonTrueDirectionX.at(iPi0Phot),pi0PhotonTrueDirectionY.at(iPi0Phot),pi0PhotonTrueDirectionZ.at(iPi0Phot));
           //std::cout << "debA" << std::endl;
        TVector3 phot2mom(pi0PhotonTrueDirectionX.at(jPi0Phot),pi0PhotonTrueDirectionY.at(jPi0Phot),pi0PhotonTrueDirectionZ.at(jPi0Phot));
           //std::cout << "debB" << std::endl;
        double angle = phot1mom.Angle(phot2mom);
        double angleDeg = angle/TMath::Pi()*180;
        double invMass = TMath::Sqrt(2*pi0PhotonTrueEnergy.at(iPi0Phot)*pi0PhotonTrueEnergy.at(jPi0Phot)*(1-TMath::Cos(angle))); 
           //std::cout << "debC" << std::endl;
        TVector3 recoPhot1Mom(pi0PhotonMatchedPfpShowerDirectionX.at(iPi0Phot),pi0PhotonMatchedPfpShowerDirectionY.at(iPi0Phot),pi0PhotonMatchedPfpShowerDirectionZ.at(iPi0Phot));
           //std::cout << "debD" << std::endl;
        TVector3 recoPhot2Mom(pi0PhotonMatchedPfpShowerDirectionX.at(jPi0Phot),pi0PhotonMatchedPfpShowerDirectionY.at(jPi0Phot),pi0PhotonMatchedPfpShowerDirectionZ.at(jPi0Phot));
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
        //std::cout << "reco angle = " << recoAngle << " reco inv mass = " << recoInvMass << std::endl;
        //Fill plots
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
  std::cout << "nPhotonPairs = " << nPhotonPairs << std::endl;
}

void test::Pi0AnalysisPlots::endJob()
{
  // Implementation of optional member function here.
}

bool test::Pi0AnalysisPlots::PassesPhotonQualityCuts(int nHitsU, int nHitsV, int nHitsW, double energy, double purity, double completeness)
{
  if(nHitsU+nHitsV+nHitsW<15) return false;
  if((nHitsU<5 && nHitsV<5) || (nHitsU<5 && nHitsW<5) || (nHitsV<5 && nHitsW<5)) return false;
  if(energy<0.1) return false; //100 MeV cut
  if(purity<0.5) return false;
  if(completeness<0.5) return false;
  return true;
}

DEFINE_ART_MODULE(test::Pi0AnalysisPlots)
