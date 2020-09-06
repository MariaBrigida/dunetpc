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
  //TH1D *fTrackLengthHist;

  unsigned int fEventID;
  unsigned int fNPFParticles;
  //unsigned int fNPrimaries;
  std::vector<bool> *fIsPrimary;
  std::vector<bool> *fIsMCPrimary;
  std::vector<int> *fPdgCode;
  std::vector<int> *fNumDaughters;
  std::vector<int> *fParentID;
  std::vector<float> *fTrueEnergy; 
  std::vector<int> *fNClusters;
  std::vector<int> *fNTracks;
  std::vector<std::vector<int>> *fCluPlane;
  std::vector<std::vector<int>> *fCluView;
  std::vector<std::vector<int>> *fCluNHits;
  std::vector<std::vector<double>> *fCluIntegral;
  //int fNPrimaryDaughters;
  //std::vector<float> *fDaughterTrackLengths;
  std::vector<std::vector<float>> *fTrackLength;
  std::vector<std::vector<float>> *fTrackStartX;
  std::vector<std::vector<float>> *fTrackStartY;
  std::vector<std::vector<float>> *fTrackStartZ;
  std::vector<std::vector<float>> *fTrackEndX;
  std::vector<std::vector<float>> *fTrackEndY;
  std::vector<std::vector<float>> *fTrackEndZ;
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
    fIsPrimary(nullptr),
    fIsMCPrimary(nullptr),
    fPdgCode(nullptr),
    fNumDaughters(nullptr),
    fParentID(nullptr),
    fTrueEnergy(nullptr),
    fNClusters(nullptr),
    fNTracks(nullptr),
    fCluPlane(nullptr),
    fCluView(nullptr),
    fCluNHits(nullptr),
    fCluIntegral(nullptr),
    fTrackLength(nullptr),
    fTrackStartX(nullptr),
    fTrackStartY(nullptr),
    fTrackStartZ(nullptr),
    fTrackEndX(nullptr),
    fTrackEndY(nullptr),
    fTrackEndZ(nullptr)
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
}

void test::pandoraAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();
  //std::cout << "debug" << std::endl;
  fNPFParticles = 0;
  //fNPrimaries = 0;
  //fNPrimaryDaughters = 0;
  //fDaughterTrackLengths->clear();
  fIsPrimary->clear();
  fPdgCode->clear();
  fNumDaughters->clear();
  fParentID->clear();
  //fTrueEnergy->clear();	//this creates a seg fault!
  fNClusters->clear();
  fNTracks->clear();
  fCluPlane->clear();
  fCluView->clear();
  fCluNHits->clear();
  fCluIntegral->clear();
  fTrackLength->clear();
  fTrackStartX->clear();
  fTrackStartY->clear();
  fTrackStartZ->clear();
  fTrackEndX->clear();
  fTrackEndY->clear();
  fTrackEndZ->clear();
  //fTrackdEdx->clear();
  //fTrackResidualRange->clear();

  //std::cout << "debug2" << std::endl;

  //Access the truth information
  if(!e.isRealData()){

    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()){
      bool isMCPrimary(false);
      for(unsigned int iMc=0; iMc< mcParticles->size(); iMc++){
        const simb::MCParticle trueParticle = mcParticles->at(iMc);
        fTrueEnergy->push_back(trueParticle.E());
        trueParticle.Process()=="primary"? isMCPrimary=true:isMCPrimary=false;
        fIsMCPrimary->push_back(isMCPrimary);
      }
    }
 }
//
//std::cout << "debug1" << std::endl;
  //Access the PFParticles from Pandora
  art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfparticleVect;
  if (e.getByLabel(fPFParticleLabel,pfparticleHandle)) //make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect,pfparticleHandle); //fill the vector with art::Ptr PFParticles
  //if(!pfparticleVect.size()) return;
  fNPFParticles = pfparticleVect.size();
  //std::cout << "debug3" << std::endl;

//std::cout << "debug2" << std::endl;
  //Access the Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVect;
  if (e.getByLabel(fPFParticleLabel,clusterHandle)) //make sure the handle is valid
    art::fill_ptr_vector(clusterVect,clusterHandle); //fill the vector
  art::FindManyP<recob::Cluster> clusterAssoc(pfparticleVect, e, fPFParticleLabel);

//std::cout << "debug3" << std::endl;
  //Access the Tracks from pandoraTrack
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> trackVect;
  //fTree->Branch("parentID",&fParentID);
  if (e.getByLabel(fTrackLabel,trackHandle)) //make sure the handle is valid
    art::fill_ptr_vector(trackVect,trackHandle); //fill the vector with art::Ptr PFParticles
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect, e, fTrackLabel);

//std::cout << "debug4" << std::endl;
  //Access the hits from pandora
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVect;
  if (e.getByLabel(fPFParticleLabel,hitHandle)) //make sure the handle is valid
    art::fill_ptr_vector(hitVect,hitHandle); //fill the vector with art::Ptr PFParticles
  //Hits from tracks:
  art::FindManyP<recob::Hit> trackHitAssoc(trackVect, e, fTrackLabel);
  

//std::cout << "debug5" << std::endl;
  //art::FindManyP<anab::Calorimetry> calorimetryAssoc(trackVect, e, fCalorimetryLabel);
  //std::cout << "i event = " << fEventID << " n pfparticles = " << fNPFParticles << std::endl;
  std::vector<int> pfpCluPlane, pfpCluView, pfpCluNHits;
  std::vector<double> pfpCluIntegral;
  std::vector<float> pfpTrackLenght, pfpTrackStartX, pfpTrackStartY, pfpTrackStartZ, pfpTrackEndX, pfpTrackEndY, pfpTrackEndZ;


//std::cout << "debug6" << std::endl;
  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
    //std::cout <<"----------------------------------->iPfp = " << iPfp << std::endl;
    iPfp++;
    fIsPrimary->push_back(pfp->IsPrimary());
    fPdgCode->push_back(pfp->PdgCode());
    fNumDaughters->push_back(pfp->NumDaughters());
    (pfp->IsPrimary())? fParentID->push_back(-1) : fParentID->push_back(pfp->Parent());
    //(pfp->IsPrimary())? fParentID->push_back(-1) : fParentID->push_back(-2);
    //fParentID->push_back(-1);
    //std::cout << "self = " << pfp->Self() << std::endl;
    //std::cout << "isprimary? " << fIsPrimary->back() << std::endl;
    //if(!(pfp->IsPrimary() && (std::abs(pfp->GetPdgCode())==14 || std::abs(pfp->GetPdgCode())==12))) continue;
    //if(!(pfp->IsPrimary())) continue;
    //fNPrimaries++;

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterAssoc.at(pfp.key());    
    fNClusters->push_back(pfpClusters.size());
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


//std::cout << "debug7" << std::endl;
    std::vector<art::Ptr<recob::Track>> pfpTracks = trackAssoc.at(pfp.key());    
    fNTracks->push_back(pfpTracks.size());
    std::cout << "tracks size = " << pfpTracks.size() << std::endl;
    if(!pfpTracks.empty()){
      int iTrk(0);
      for(const art::Ptr<recob::Track> &trk:pfpTracks){
        //std::cout << "iTrk = " << iTrk << " fTrackLength size = " << fTrackLength->size() << std::endl;
        std::cout << "--------------new track------------------" << std::endl;
	iTrk++;
        pfpTrackLenght.push_back(trk->Length());
        pfpTrackStartX.push_back(trk->Start().X());
        pfpTrackStartY.push_back(trk->Start().Y());
        pfpTrackStartZ.push_back(trk->Start().Z());
        pfpTrackEndX.push_back(trk->End().X());
        pfpTrackEndY.push_back(trk->End().Y());
        pfpTrackEndZ.push_back(trk->End().Z());
	//fTrackLengthHist->Fill(trk->Length());

        /*std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssoc.at(trk.key());
        if(!trackHits.empty()){
          int iTrkHit(0);
          //std::cout << "Track n. " << iTrk << " has NHits = " << trackHits.size() << std::endl;
          for(const art::Ptr<recob::Hit> &trkHit:trackHits){
            
            //std::cout << "iTrkHit = " << iTrkHit << " integral = " << trkHit->Integral() << " view = " << trkHit->View()<< std::endl; 
            iTrkHit++;
      
          }
        } */

//DEBUG SLAVIC CRPs
//
    std::vector<unsigned> hitsTpcId( fGeom->Nplanes() );
    // loop over the planes 
    for(size_t i_plane = 0; i_plane<fGeom->Nplanes(); i_plane++) {
       // get hits in this plane
       //auto hits = trackUtil.GetRecoTrackHitsFromPlane( trk, e, fTrackModuleLabel, i_plane );


      std::vector<const recob::Hit*> trackHits;
      //if( planeID > 2 ){
      //  std::cout << "Please input plane 0, 1, or 2" << std::endl;
      //  return trackHits;
      //}
   
      auto recoTracks = e.getValidHandle<std::vector<recob::Track> >(fTrackLabel);
      art::FindManyP<recob::Hit> findHits(recoTracks,e,fTrackLabel);
      std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(trk->ID());
   
      for(const art::Ptr<recob::Hit> hit : inputHits){
        unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
        if( thePlane != i_plane ) continue;
         
        trackHits.push_back(hit.get());
   
      }

      //if( fLogLevel >= 3 ){
         //std::cout<<"Hits in plane "<<i_plane<<" "<<trackHits.size()<<std::endl;
      // }
       std::vector<unsigned> plane_hits_tpcid( trackHits.size() );
         
       // l(op over hits
       for(size_t i_hit = 0; i_hit<trackHits.size(); i_hit++ ){
         plane_hits_tpcid[i_hit] = trackHits[i_hit]->WireID().TPC;
         std::cout << "plane " << i_plane << " WireID().TPC = " << trackHits[i_hit]->WireID().TPC << std::endl;
       }
         
       int this_tpcid = -1;
       auto start = plane_hits_tpcid.begin();
       auto end   = plane_hits_tpcid.end();
         
       if( std::equal(start + 1, end, start)) {
         this_tpcid = (int)(*start);
       } //
       else {
       //  if( fLogLevel>= 2){
           std::cout<< this_tpcid << " hits for the same plane have mixed TPC IDs. Skipping...\n";
       //  }
       }
      
       /*if( this_tpcid < 0 ){
         hitsTpcId.clear();
         break;
       }
       else{
         hitsTpcId[i_plane] = (unsigned)this_tpcid;
       }*/
// end plane loop	
   }
//
///////////


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


//std::cout << "debug8" << std::endl;
    fCluPlane->push_back(pfpCluPlane);
    fCluView->push_back(pfpCluView);
    fCluNHits->push_back(pfpCluNHits);
    fCluIntegral->push_back(pfpCluIntegral);
    fTrackLength->push_back(pfpTrackLenght);
    fTrackStartX->push_back(pfpTrackStartX);
    fTrackStartY->push_back(pfpTrackStartY);
    fTrackStartZ->push_back(pfpTrackStartZ);
    fTrackEndX->push_back(pfpTrackEndX);
    fTrackEndY->push_back(pfpTrackEndY);
    fTrackEndZ->push_back(pfpTrackEndZ);
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
   
//std::cout << "debug9" << std::endl;

    //std::cout << "nclusters = " << pfpClusters.size() << " " << fNClusters->back() << std::endl;
    //std::cout << "nClusters = " << pfpClusters.size() << std::endl;



    //neutrinoID = pfp->Self();
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

  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
  fTree->Branch("isPrimary",&fIsPrimary);
  fTree->Branch("pdgCode",&fPdgCode);
  fTree->Branch("numDaughters",&fNumDaughters);
  fTree->Branch("parentID",&fParentID);
  fTree->Branch("isMCPrimary",&fIsMCPrimary);
  fTree->Branch("trueEnergy",&fTrueEnergy);
  fTree->Branch("nClusters",&fNClusters);
  fTree->Branch("nTracks",&fNTracks);
  fTree->Branch("cluPlane",&fCluPlane);
  fTree->Branch("cluView",&fCluView);
  fTree->Branch("cluNHits",&fCluNHits);
  fTree->Branch("cluIntegral",&fCluIntegral);
  //fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
  //fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/I");
  fTree->Branch("trackLength",&fTrackLength);
  fTree->Branch("trackStartX",&fTrackStartX);
  fTree->Branch("trackStartY",&fTrackStartY);
  fTree->Branch("trackStartZ",&fTrackStartZ);
  fTree->Branch("trackEndX",&fTrackEndX);
  fTree->Branch("trackEndY",&fTrackEndY);
  fTree->Branch("trackEndZ",&fTrackEndZ);
  //fTree->Branch("trackdEdx",&fTrackdEdx);
  //fTree->Branch("trackResidualRange",&fTrackResidualRange);

  //fTrackLengthHist = tfs->make<TH1D>("trackLengthHist","reconstructed track lengths",110,-100,1000); 

}

void test::pandoraAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::pandoraAnalysis)
