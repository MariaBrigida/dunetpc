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
  std::vector<float> *fMCParticleTrueEnergy; 
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
  //int fNPrimaryDaughters;
  //std::vector<float> *fDaughterTrackLengths;
  std::vector<std::vector<float>> *fPFPTrackLength;
  std::vector<std::vector<float>> *fPFPTrackStartX;
  std::vector<std::vector<float>> *fPFPTrackStartY;
  std::vector<std::vector<float>> *fPFPTrackStartZ;
  std::vector<std::vector<float>> *fPFPTrackEndX;
  std::vector<std::vector<float>> *fPFPTrackEndY;
  std::vector<std::vector<float>> *fPFPTrackEndZ;
  std::vector<std::vector<int>> *fPFPTrackID;
  //std::vector<std::vector<float>> *fTrackdEdx;
  //std::vector<std::vector<float>> *fTrackResidualRange;
  std::string fTruthLabel;
  std::string fTrackLabel;
  //std::string fClusterLabel;
  //std::string fCalorimetryLabel;
  std::string fPFParticleLabel;

  const geo::Geometry* fGeom;
  //protoana::ProtoDUNETrackUtils trackUtil;
};


test::pandoraAnalysis::pandoraAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fNMCParticles(nullptr),
    fNPFParticles(nullptr),	
    fPFPIsPrimary(nullptr),
    fMCIsPrimary(nullptr),
    fMCParticlePdgCode(nullptr),
    fPFPPdgCode(nullptr),
    fPFPNDaughters(nullptr),
    fPFPParentID(nullptr),
    fMCParticleTrueEnergy(nullptr),
    fMCParticleTrackID(nullptr),
    fMCParticleMotherTrackID(nullptr),
    fMCParticleStartProcess(nullptr),
    fMCParticleEndProcess(nullptr),
    fMCParticleNTrajectoryPoint(nullptr),
    fPFPNClusters(nullptr),
    fPFPNTracks(nullptr),
    fPFPCluPlane(nullptr),
    fPFPCluView(nullptr),
    fPFPCluNHits(nullptr),
    fPFPCluIntegral(nullptr),
    fPFPTrackLength(nullptr),
    fPFPTrackStartX(nullptr),
    fPFPTrackStartY(nullptr),
    fPFPTrackStartZ(nullptr),
    fPFPTrackEndX(nullptr),
    fPFPTrackEndY(nullptr),
    fPFPTrackEndZ(nullptr),
    fPFPTrackID(nullptr)
    //fTrackdEdx(nullptr),
    //fTrackResidualRange(nullptr)
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fTruthLabel = p.get<std::string>("TruthLabel");
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  /*fCalorimetryLabel = p.get<std::string>("CalorimetryLabel");
  fClusterLabel = p.get<std::string>("ClusterLabel");*/
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  /*for(int iplane=0; iplane<3; iplane ++){
    fTrackHitsTPCs.at(iplane)->clear();
  }*/
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
  //std::cout << "e.isRealData = " << e.isRealData() << std::endl;
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
  art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfparticleVect;
  if (e.getByLabel(fPFParticleLabel,pfparticleHandle)) //make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect,pfparticleHandle); //fill the vector with art::Ptr PFParticles
  //if(!pfparticleVect.size()) return;
  fNPFParticles = pfparticleVect.size();

  //Access the Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVect;
  if (e.getByLabel(fPFParticleLabel,clusterHandle)) //make sure the handle is valid
    art::fill_ptr_vector(clusterVect,clusterHandle); //fill the vector
  art::FindManyP<recob::Cluster> clusterAssoc(pfparticleVect, e, fPFParticleLabel);

  //Access the Tracks from pandoraTrack
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> trackVect;
  //fTree->Branch("parentID",&fPFPParentID);
  if (e.getByLabel(fTrackLabel,trackHandle)) //make sure the handle is valid
    art::fill_ptr_vector(trackVect,trackHandle); //fill the vector with art::Ptr PFParticles
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect, e, fTrackLabel);

  //Access the hits from pandora
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVect;
  if (e.getByLabel(fPFParticleLabel,hitHandle)) //make sure the handle is valid
    art::fill_ptr_vector(hitVect,hitHandle); //fill the vector with art::Ptr PFParticles
  //Hits from tracks:
  art::FindManyP<recob::Hit> trackHitAssoc(trackVect, e, fTrackLabel);
  

  //art::FindManyP<anab::Calorimetry> calorimetryAssoc(trackVect, e, fCalorimetryLabel);
  std::vector<int> pfpCluPlane, pfpCluView, pfpCluNHits;
  std::vector<double> pfpCluIntegral;
  std::vector<float> pfpTrackLenght, pfpTrackStartX, pfpTrackStartY, pfpTrackStartZ, pfpTrackEndX, pfpTrackEndY, pfpTrackEndZ;
  std::vector<int> pfpTrackID;

//std::cout << "debug6" << std::endl;
  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
    //std::cout <<"----------------------------------->iPfp = " << iPfp << std::endl;
    iPfp++;
    fPFPIsPrimary->push_back(pfp->IsPrimary());
    fPFPPdgCode->push_back(pfp->PdgCode());
    fPFPNDaughters->push_back(pfp->NumDaughters());
    (pfp->IsPrimary())? fPFPParentID->push_back(-1) : fPFPParentID->push_back(pfp->Parent());
    //(pfp->IsPrimary())? fPFPParentID->push_back(-1) : fPFPParentID->push_back(-2);
    //fPFPParentID->push_back(-1);
    //std::cout << "self = " << pfp->Self() << std::endl;
    //std::cout << "isprimary? " << fPFPIsPrimary->back() << std::endl;
    //if(!(pfp->IsPrimary() && (std::abs(pfp->GetPdgCode())==14 || std::abs(pfp->GetPdgCode())==12))) continue;
    //if(!(pfp->IsPrimary())) continue;
    //fNPrimaries++;

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterAssoc.at(pfp.key());    
    fPFPNClusters->push_back(pfpClusters.size());
    if(!pfpClusters.empty()){
      int iClu(0);
      for(const art::Ptr<recob::Cluster> &clu:pfpClusters){
       //std::cout << "iClu = " << iClu <<" view = " << clu->View() << " plane = " << clu->Plane().asPlaneID().toString() << " tpc: " << clu->Plane().asTPCID().toString()<< " cryostat: " << clu->Plane().asCryostatID().toString() <<  " width = " << clu->Width() << std::endl; 
        //std::cout << "start charge: " << clu->StartCharge() << " angle: " << clu->StartAngle() << " opening angle: " << clu->StartOpeningAngle() << std::endl;
        //std::cout << "end charge: " << clu->EndCharge() << " angle: " << clu->EndAngle() << " opening angle: " << clu->EndOpeningAngle() << std::endl;
        //std::cout << "integral: " << clu->Integral() << " integral std dev: " << clu->IntegralStdDev() << " integral average: " << clu->IntegralAverage() << std::endl;
        //std::cout << "nhits: " << clu->NHits() << std::endl;
        //pfpCluPlane.push_back(clu->Plane().toString());
	pfpCluView.push_back(clu->View());
        pfpCluNHits.push_back(clu->NHits());
	pfpCluIntegral.push_back(clu->Integral());
        iClu++;

      }
    }


    std::vector<art::Ptr<recob::Track>> pfpTracks = trackAssoc.at(pfp.key());    
    fPFPNTracks->push_back(pfpTracks.size());
    std::cout << "tracks size = " << pfpTracks.size() << std::endl;
    std::vector<unsigned int> trkTrackHitsTPCs_Plane0,trkTrackHitsTPCs_Plane1,trkTrackHitsTPCs_Plane2;
    if(!pfpTracks.empty()){
      int iTrk(0);
      for(const art::Ptr<recob::Track> &trk:pfpTracks){
        //std::cout << "iTrk = " << iTrk << " fPFPTrackLength size = " << fPFPTrackLength->size() << std::endl;
        std::cout << "--------------new track------------------" << std::endl;
	iTrk++;
        pfpTrackLenght.push_back(trk->Length());
        pfpTrackStartX.push_back(trk->Start().X());
        pfpTrackStartY.push_back(trk->Start().Y());
        pfpTrackStartZ.push_back(trk->Start().Z());
        pfpTrackEndX.push_back(trk->End().X());
        pfpTrackEndY.push_back(trk->End().Y());
        pfpTrackEndZ.push_back(trk->End().Z());
        pfpTrackID.push_back(trk->ID());
	//fPFPTrackLengthHist->Fill(trk->Length());

        /*std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(trk.key());
        if(!trackHits.empty()){
          int iTrkHit(0);
          //std::cout << "Track n. " << iTrk << " has NHits = " << trackHits.size() << std::endl;
          for(const art::Ptr<recob::Hit> &trkHit:trackHits){
            
            //std::cout << "iTrkHit = " << iTrkHit << " integral = " << trkHit->Integral() << " view = " << trkHit->View()<< std::endl; 
            iTrkHit++;
      
          }
        } */

		

        //std::vector<art::Ptr<anab::Calorimetry>> trackCalo = calorimetryAssoc.at(trk.key());
        //for(const art::Ptr<anab::Calorimetry> &cal:trackCalo){
	  //if(!cal->PlaneID().isValid) continue;	
	  //int planenum = cal->PlaneID().Plane;
	  //if(planenum!=2)continue;
	  //fTrackdEdx->push_back(cal->dEdx());
	  //fTrackResidualRange->push_back(cal->ResidualRange());
        //}
      }
    }


    fPFPCluPlane->push_back(pfpCluPlane);
    fPFPCluView->push_back(pfpCluView);
    fPFPCluNHits->push_back(pfpCluNHits);
    fPFPCluIntegral->push_back(pfpCluIntegral);
    fPFPTrackLength->push_back(pfpTrackLenght);
    fPFPTrackStartX->push_back(pfpTrackStartX);
    fPFPTrackStartY->push_back(pfpTrackStartY);
    fPFPTrackStartZ->push_back(pfpTrackStartZ);
    fPFPTrackEndX->push_back(pfpTrackEndX);
    fPFPTrackEndY->push_back(pfpTrackEndY);
    fPFPTrackEndZ->push_back(pfpTrackEndZ);
    fPFPTrackID->push_back(pfpTrackID);
    pfpCluPlane.clear();
    pfpCluView.clear();
    pfpCluNHits.clear();
    pfpCluIntegral.clear();
    pfpTrackLenght.clear();
    pfpTrackStartX.clear();
    pfpTrackStartY.clear();
    pfpTrackStartZ.clear();
    pfpTrackEndX.clear();
    pfpTrackEndY.clear();
    pfpTrackEndZ.clear();
  }

  //fNPrimaryDaughters = pfp->NumDaughters();
  //std::cout << "ID = " << neutrinoID <<std::endl;  
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
  fTree->Branch("mcParticleEndPosition",&fMCParticleEndPosition)
  fTree->Branch("mcParticleEndMomentum",&fMCParticleEndMomentum)

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
