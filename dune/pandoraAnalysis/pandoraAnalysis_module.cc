////////////////////////////////////////////////////////////////////////
// Class:       pandoraAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        pandoraAnalysis_module.cc
//
// Generated at Fri Aug  7 15:01:35 2020 by Maria Brigida Brunetti using cetskelgen
// from cetlib version v3_10_00.
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

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"


#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
  class pandoraAnalysis;
}


class test::pandoraAnalysis : public art::EDAnalyzer {
public:
  explicit pandoraAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pandoraAnalysis(pandoraAnalysis const&) = delete;
  pandoraAnalysis(pandoraAnalysis&&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis const&) = delete;
  pandoraAnalysis& operator=(pandoraAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree *fTree;

  unsigned int fEventID;
  unsigned int fNMCParticles;
  unsigned int fNPFParticles;

  std::vector<bool> *fMCIsPrimary;
  std::vector<int> *fMCParticlePdgCode;
  std::vector<double> *fMCParticleTrueEnergy; 
  std::vector<int> *fMCParticleTrackID; 
  std::vector<int> *fMCParticleMotherTrackID; 
  std::vector<std::string> *fMCParticleStartProcess; 
  std::vector<std::string> *fMCParticleEndProcess; 
  std::vector<int> *fMCParticleNTrajectoryPoint; 
  std::vector<TLorentzVector> *fMCParticleStartPosition;
  std::vector<TLorentzVector> *fMCParticleStartMomentum;
  std::vector<TLorentzVector> *fMCParticleEndPosition;
  std::vector<TLorentzVector> *fMCParticleEndMomentum;


  std::vector<bool> *fPFPIsPrimary;
  std::vector<int> *fPFPParentID;
  std::vector<int> *fPFPPdgCode;
  std::vector<int> *fPFPNDaughters;
  std::vector<int> *fPFPNClusters;
  std::vector<int> *fPFPNTracks;
  std::vector<std::vector<int>> *fPFPCluPlane;
  std::vector<std::vector<int>> *fPFPCluView;
  std::vector<std::vector<int>> *fPFPCluNHits;
  std::vector<std::vector<double>> *fPFPCluIntegral;


  std::vector<int> *fPFPTrackID;
  std::vector<double> *fPFPTrackLength;
  std::vector<double> *fPFPTrackStartX;
  std::vector<double> *fPFPTrackStartY;
  std::vector<double> *fPFPTrackStartZ;
  std::vector<double> *fPFPTrackVertexX;
  std::vector<double> *fPFPTrackVertexY;
  std::vector<double> *fPFPTrackVertexZ;
  std::vector<double> *fPFPTrackEndX;
  std::vector<double> *fPFPTrackEndY;
  std::vector<double> *fPFPTrackEndZ;
  std::vector<double> *fPFPTrackTheta;
  std::vector<double> *fPFPTrackPhi;
  std::vector<double> *fPFPTrackZenithAngle;
  std::vector<double> *fPFPTrackAzimuthAngle;
  std::vector<double> *fPFPTrackStartDirectionX;
  std::vector<double> *fPFPTrackStartDirectionY;
  std::vector<double> *fPFPTrackStartDirectionZ;
  std::vector<double> *fPFPTrackVertexDirectionX;
  std::vector<double> *fPFPTrackVertexDirectionY;
  std::vector<double> *fPFPTrackVertexDirectionZ;
  std::vector<double> *fPFPTrackEndDirectionX;
  std::vector<double> *fPFPTrackEndDirectionY;
  std::vector<double> *fPFPTrackEndDirectionZ;
  std::vector<float> *fPFPTrackChi2;
  std::vector<int> *fPFPTrackNdof;

  std::vector<int>    *fPFPShowerID;
  std::vector<int>    *fPFPShowerBestPlane;
  std::vector<double> *fPFPShowerDirectionX;
  std::vector<double> *fPFPShowerDirectionY;
  std::vector<double> *fPFPShowerDirectionZ;
  std::vector<double> *fPFPShowerDirectionErrX;
  std::vector<double> *fPFPShowerDirectionErrY;
  std::vector<double> *fPFPShowerDirectionErrZ;
  std::vector<double> *fPFPShowerStartX;
  std::vector<double> *fPFPShowerStartY;
  std::vector<double> *fPFPShowerStartZ;
  std::vector<double> *fPFPShowerStartErrX;
  std::vector<double> *fPFPShowerStartErrY;
  std::vector<double> *fPFPShowerStartErrZ;
  std::vector<std::vector<double>> *fPFPShowerEnergy;		//A value per shower per plane!
  std::vector<std::vector<double>> *fPFPShowerEnergyErr;
  std::vector<std::vector<double>> *fPFPShowerMIPEnergy;
  std::vector<std::vector<double>> *fPFPShowerMIPEnergyErr;
  std::vector<double> *fPFPShowerLength;
  std::vector<double> *fPFPShowerOpenAngle;
  std::vector<std::vector<double>> *fPFPShowerdEdx;
  std::vector<std::vector<double>> *fPFPShowerdEdxErr;

  std::string fTruthLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;

  const geo::Geometry* fGeom;
};


test::pandoraAnalysis::pandoraAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fMCIsPrimary(nullptr),
    fMCParticlePdgCode(nullptr),
    fMCParticleTrueEnergy(nullptr),
    fMCParticleTrackID(nullptr),
    fMCParticleMotherTrackID(nullptr),
    fMCParticleStartProcess(nullptr),
    fMCParticleEndProcess(nullptr),
    fMCParticleNTrajectoryPoint(nullptr),
    fPFPIsPrimary(nullptr),
    fPFPParentID(nullptr),
    fPFPPdgCode(nullptr),
    fPFPNDaughters(nullptr),
    fPFPNClusters(nullptr),
    fPFPNTracks(nullptr),
    fPFPCluPlane(nullptr),
    fPFPCluView(nullptr),
    fPFPCluNHits(nullptr),
    fPFPCluIntegral(nullptr),
    fPFPTrackID(nullptr),
    fPFPTrackLength(nullptr),
    fPFPTrackStartX(nullptr),
    fPFPTrackStartY(nullptr),
    fPFPTrackStartZ(nullptr),
    fPFPTrackVertexX(nullptr),
    fPFPTrackVertexY(nullptr),
    fPFPTrackVertexZ(nullptr),
    fPFPTrackEndX(nullptr),
    fPFPTrackEndY(nullptr),
    fPFPTrackEndZ(nullptr),
    fPFPTrackTheta(nullptr),
    fPFPTrackPhi(nullptr),
    fPFPTrackZenithAngle(nullptr),
    fPFPTrackAzimuthAngle(nullptr),
    fPFPTrackStartDirectionX(nullptr),
    fPFPTrackStartDirectionY(nullptr),
    fPFPTrackStartDirectionZ(nullptr),
    fPFPTrackVertexDirectionX(nullptr),
    fPFPTrackVertexDirectionY(nullptr),
    fPFPTrackVertexDirectionZ(nullptr),
    fPFPTrackEndDirectionX(nullptr),
    fPFPTrackEndDirectionY(nullptr),
    fPFPTrackEndDirectionZ(nullptr),
    fPFPTrackChi2(nullptr),
    fPFPTrackNdof(nullptr),
    fPFPShowerID(nullptr),
    fPFPShowerDirectionX(nullptr),
    fPFPShowerDirectionY(nullptr),
    fPFPShowerDirectionZ(nullptr),
    fPFPShowerDirectionErrX(nullptr),
    fPFPShowerDirectionErrY(nullptr),
    fPFPShowerDirectionErrZ(nullptr),
    fPFPShowerStartX(nullptr),
    fPFPShowerStartY(nullptr),
    fPFPShowerStartZ(nullptr),
    fPFPShowerStartErrX(nullptr),
    fPFPShowerStartErrY(nullptr),
    fPFPShowerStartErrZ(nullptr),
    fPFPShowerEnergy(nullptr),      //An energy value per shower per plane
    fPFPShowerEnergyErr(nullptr),
    fPFPShowerMIPEnergy(nullptr),
    fPFPShowerMIPEnergyErr(nullptr),
    fPFPShowerLength(nullptr),
    fPFPShowerOpenAngle(nullptr),
    fPFPShowerdEdx(nullptr),
    fPFPShowerdEdxErr(nullptr)

    //fTrackdEdx(nullptr),
    //fTrackResidualRange(nullptr)
  // More initializers here.
{
  fTruthLabel = p.get<std::string>("TruthLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fShowerLabel = p.get<std::string>("ShowerLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}

void test::pandoraAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  fNMCParticles = 0;
  fNPFParticles = 0;
  //fNPrimaries = 0;
  //fNPrimaryDaughters = 0;
  //fDaughterTrackLengths->clear();
  fPFPIsPrimary->clear();
  fMCParticlePdgCode->clear();
  fPFPPdgCode->clear();
  fPFPNDaughters->clear();
  fPFPParentID->clear();
  fMCParticleTrueEnergy->clear();
  fMCParticleTrackID->clear();
  fMCParticleMotherTrackID->clear();
  fMCParticleStartProcess->clear();
  fMCParticleEndProcess->clear();
  fMCParticleNTrajectoryPoint->clear();
  fPFPNClusters->clear();
  fPFPNTracks->clear();
  fPFPCluPlane->clear();
  fPFPCluView->clear();
  fPFPCluNHits->clear();
  fPFPCluIntegral->clear();
  fPFPTrackLength->clear();
  fPFPTrackStartX->clear();
  fPFPTrackStartY->clear();
  fPFPTrackStartZ->clear();
  fPFPTrackEndX->clear();
  fPFPTrackEndY->clear();
  fPFPTrackEndZ->clear();
  fPFPTrackID->clear();
  
  //Access the truth information
  if(!e.isRealData()){
    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
      fNMCParticles=mcParticles->size();
      bool isMCPrimary(false);
      for(unsigned int iMc=0; iMc< mcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        fMCParticleTrueEnergy->push_back(trueParticle.E());
        fMCParticlePdgCode->push_back(trueParticle.PdgCode());
        fMCParticleTrackID->push_back(trueParticle.TrackId());
        fMCParticleMotherTrackID->push_back(trueParticle.Mother());
        fMCParticleStartProcess->push_back(trueParticle.Process());
        fMCParticleEndProcess->push_back(trueParticle.EndProcess());
        fMCParticleNTrajectoryPoint->push_back(trueParticle.NumberTrajectoryPoints());
        trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
        fMCIsPrimary->push_back(isMCPrimary);
      }
    }
  }

  //Access the PFParticles from Pandora
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e,fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();
  if(!fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    return;
  }

  //Access the Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVect;
  if (e.getByLabel(fPFParticleLabel,clusterHandle))
    art::fill_ptr_vector(clusterVect,clusterHandle);
  art::FindManyP<recob::Cluster> clusterAssoc(pfparticleVect, e, fPFParticleLabel);

  //Access the hits from pandora
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVect;
  if (e.getByLabel(fPFParticleLabel,hitHandle)) 
    art::fill_ptr_vector(hitVect,hitHandle); 
  //Hits from tracks:
  //art::FindManyP<recob::Hit> trackHitAssoc(trackVect, e, fTrackLabel);
  

  std::vector<int> pfpCluPlane, pfpCluView, pfpCluNHits;
  std::vector<double> pfpCluIntegral;
  std::vector<int> pfpTrackID;

  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
    iPfp++;
    fPFPIsPrimary->push_back(pfp->IsPrimary());
    fPFPPdgCode->push_back(pfp->PdgCode());
    fPFPNDaughters->push_back(pfp->NumDaughters());
    (pfp->IsPrimary())? fPFPParentID->push_back(-1) : fPFPParentID->push_back(pfp->Parent());

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterAssoc.at(pfp.key());    
    fPFPNClusters->push_back(pfpClusters.size());
    if(!pfpClusters.empty()){
      int iClu(0);
      for(const art::Ptr<recob::Cluster> &clu:pfpClusters){
	pfpCluView.push_back(clu->View());
        pfpCluNHits.push_back(clu->NHits());
	pfpCluIntegral.push_back(clu->Integral());
        iClu++;
      }
    }

    art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp,e,fPFParticleLabel, fTrackLabel); 
    if(track){
      fPFPTrackID->push_back(track->ID());
      fPFPTrackLength->push_back(track->Length());
      fPFPTrackStartX->push_back(track->Start().X());
      fPFPTrackStartY->push_back(track->Start().Y());
      fPFPTrackStartZ->push_back(track->Start().Z());
      fPFPTrackVertexX->push_back(track->Vertex().X());
      fPFPTrackVertexY->push_back(track->Vertex().Y());
      fPFPTrackVertexZ->push_back(track->Vertex().Z());
      fPFPTrackEndX->push_back(track->End().X());
      fPFPTrackEndY->push_back(track->End().Y());
      fPFPTrackEndZ->push_back(track->End().Z());
      fPFPTrackTheta->push_back(track->Theta());
      fPFPTrackPhi->push_back(track->Phi());
      fPFPTrackZenithAngle->push_back(track->ZenithAngle());
      fPFPTrackAzimuthAngle->push_back(track->AzimuthAngle());
      fPFPTrackStartDirectionX->push_back(track->StartDirection().X());
      fPFPTrackStartDirectionY->push_back(track->StartDirection().Y());
      fPFPTrackStartDirectionZ->push_back(track->StartDirection().Z());
      fPFPTrackVertexDirectionX->push_back(track->VertexDirection().X());
      fPFPTrackVertexDirectionY->push_back(track->VertexDirection().Y());
      fPFPTrackVertexDirectionZ->push_back(track->VertexDirection().Z());
      fPFPTrackEndDirectionX->push_back(track->EndDirection().X());
      fPFPTrackEndDirectionY->push_back(track->EndDirection().Y());
      fPFPTrackEndDirectionZ->push_back(track->EndDirection().Z());
      fPFPTrackChi2->push_back(track->Chi2());
      fPFPTrackNdof->push_back(track->Ndof());
    }

    art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp,e,fPFParticleLabel, fShowerLabel); 
    if(shower){
      fPFPShowerID->push_back(shower->ID());
      fPFPShowerBestPlane->push_back(shower->best_plane());
      fPFPShowerDirectionX->push_back(shower->Direction().X());
      fPFPShowerDirectionY->push_back(shower->Direction().Y());
      fPFPShowerDirectionZ->push_back(shower->Direction().Z());
      fPFPShowerDirectionErrX->push_back(shower->DirectionErr().X());
      fPFPShowerDirectionErrY->push_back(shower->DirectionErr().Y());
      fPFPShowerDirectionErrZ->push_back(shower->DirectionErr().Z());
      fPFPShowerStartX->push_back(shower->ShowerStart().X());
      fPFPShowerStartY->push_back(shower->ShowerStart().Y());
      fPFPShowerStartZ->push_back(shower->ShowerStart().Z());
      fPFPShowerStartErrX->push_back(shower->ShowerStartErr().X());
      fPFPShowerStartErrY->push_back(shower->ShowerStartErr().Y());
      fPFPShowerStartErrZ->push_back(shower->ShowerStartErr().Z());
      fPFPShowerEnergy->push_back(shower->Energy());
      fPFPShowerEnergyErr->push_back(shower->EnergyErr());
      fPFPShowerMIPEnergy->push_back(shower->MIPEnergy());
      fPFPShowerMIPEnergyErr->push_back(shower->MIPEnergyErr());
      fPFPShowerLength->push_back(shower->Length());
      fPFPShowerOpenAngle->push_back(shower->OpenAngle());
      fPFPShowerdEdx->push_back(shower->dEdx());
      fPFPShowerdEdxErr->push_back(shower->dEdxErr());
    }

    fPFPCluPlane->push_back(pfpCluPlane);
    fPFPCluView->push_back(pfpCluView);
    fPFPCluNHits->push_back(pfpCluNHits);
    fPFPCluIntegral->push_back(pfpCluIntegral);
    pfpCluPlane.clear();
    pfpCluView.clear();
    pfpCluNHits.clear();
    pfpCluIntegral.clear();
  }

  fTree->Fill();
}

void test::pandoraAnalysis::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  //Event branches
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  //MC truth branches
  fTree->Branch("mcIsMCPrimary",&fMCIsPrimary);
  fTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode);
  fTree->Branch("mcParticleTrueEnergy",&fMCParticleTrueEnergy);
  fTree->Branch("mcParticleTrackID",&fMCParticleTrackID);
  fTree->Branch("mcParticleMotherTrackID",&fMCParticleMotherTrackID);
  fTree->Branch("mcParticleStartProcess",&fMCParticleStartProcess);
  fTree->Branch("mcParticleEndProcess",&fMCParticleEndProcess);
  fTree->Branch("mcParticleNTrajectoryPoints",&fMCParticleNTrajectoryPoint);
  fTree->Branch("mcParticleStartPosition",&fMCParticleStartPosition);
  fTree->Branch("mcParticleStartMomentum",&fMCParticleStartMomentum);
  fTree->Branch("mcParticleEndPosition",&fMCParticleEndPosition);
  fTree->Branch("mcParticleEndMomentum",&fMCParticleEndMomentum);

  //PFP branches
  fTree->Branch("pfpIsPrimary",&fPFPIsPrimary);
  fTree->Branch("pfpParentID",&fPFPParentID);
  fTree->Branch("pfpPdgCode",&fPFPPdgCode);
  fTree->Branch("pfpNDaughters",&fPFPNDaughters);
  fTree->Branch("pfpNClusters",&fPFPNClusters);
  fTree->Branch("pfpNTracks",&fPFPNTracks);
  fTree->Branch("pfpCluPlane",&fPFPCluPlane);
  fTree->Branch("pfpCluView",&fPFPCluView);
  fTree->Branch("pfpCluNHits",&fPFPCluNHits);
  fTree->Branch("pfpCluIntegral",&fPFPCluIntegral);
  fTree->Branch("pfpTrackLength",&fPFPTrackLength);
  fTree->Branch("pfpTrackStartX",&fPFPTrackStartX);
  fTree->Branch("pfpTrackStartY",&fPFPTrackStartY);
  fTree->Branch("pfpTrackStartZ",&fPFPTrackStartZ);
  fTree->Branch("pfpTrackEndX",&fPFPTrackEndX);
  fTree->Branch("pfpTrackEndY",&fPFPTrackEndY);
  fTree->Branch("pfpTrackEndZ",&fPFPTrackEndZ);
  fTree->Branch("pfpTrackID",&fPFPTrackID);
  //fTree->Branch("trackdEdx",&fTrackdEdx);
  //fTree->Branch("trackResidualRange",&fTrackResidualRange);

  //fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  //fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/I");


  //fPFPTrackLengthHist = tfs->make<TH1D>("trackLengthHist","reconstructed track lengths",110,-100,1000); 

}

void test::pandoraAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::pandoraAnalysis)
