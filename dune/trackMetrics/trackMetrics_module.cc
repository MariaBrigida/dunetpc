////////////////////////////////////////////////////////////////////////
// Class:       trackMetrics
// Plugin Type: analyzer (art v3_05_01)
// File:        trackMetrics_module.cc
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
  class trackMetrics;
}


class test::trackMetrics : public art::EDAnalyzer {
public:
  explicit trackMetrics(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  trackMetrics(trackMetrics const&) = delete;
  trackMetrics(trackMetrics&&) = delete;
  trackMetrics& operator=(trackMetrics const&) = delete;
  trackMetrics& operator=(trackMetrics&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree *fTree;
  TH1D *fTrackCRPsize_Plane0;
  TH1D *fTrackCRPsize_5hits_Plane0;
  TH1D *fTrackCRPsize_10hits_Plane0;
  TH1D *fTrackCRPsize_Plane1;
  TH1D *fTrackCRPsize_5hits_Plane1;
  TH1D *fTrackCRPsize_10hits_Plane1;
  TH1D *fTrackCRPsize_Plane2;
  TH1D *fTrackCRPsize_5hits_Plane2;
  TH1D *fTrackCRPsize_10hits_Plane2;

  unsigned int fEventID;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCs_Plane0;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCs_Plane1;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCs_Plane2;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCcounts_Plane0;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCcounts_Plane1;
  std::vector<std::vector<unsigned int>> *fTrackHitsTPCcounts_Plane2;
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  const geo::Geometry* fGeom;
};


test::trackMetrics::trackMetrics(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fTrackHitsTPCs_Plane0(nullptr),
    fTrackHitsTPCs_Plane1(nullptr),
    fTrackHitsTPCs_Plane2(nullptr),
    fTrackHitsTPCcounts_Plane0(nullptr),
    fTrackHitsTPCcounts_Plane1(nullptr),
    fTrackHitsTPCcounts_Plane2(nullptr)
{
  fPFParticleLabel = p.get<std::string>("PFParticleLabel");
  fTrackLabel = p.get<std::string>("TrackLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}

void test::trackMetrics::analyze(art::Event const& e)
{

  fEventID = e.id().event();
  fTrackHitsTPCs_Plane0->clear();
  fTrackHitsTPCs_Plane1->clear();
  fTrackHitsTPCs_Plane2->clear();
  fTrackHitsTPCcounts_Plane0->clear();
  fTrackHitsTPCcounts_Plane1->clear();
  fTrackHitsTPCcounts_Plane2->clear();
  
  std::cout << "== EVENT " << fEventID <<" ======================================================" << std::endl;

  //Access the PFParticles from Pandora
  art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfparticleVect;
  if (e.getByLabel(fPFParticleLabel,pfparticleHandle)) //make sure the handle is valid
    art::fill_ptr_vector(pfparticleVect,pfparticleHandle); //fill the vector with art::Ptr PFParticles

  //Access the Tracks from pandoraTrack
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> trackVect;
  if (e.getByLabel(fTrackLabel,trackHandle)) //make sure the handle is valid
    art::fill_ptr_vector(trackVect,trackHandle); //fill the vector with art::Ptr PFParticles
  art::FindManyP<recob::Track> trackAssoc(pfparticleVect, e, fTrackLabel);

  //Loop over the PFPs
  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect){
    std::cout <<"-- PFP " << iPfp << " --------------------------------------" <<std::endl;
    iPfp++;

    std::vector<art::Ptr<recob::Track>> pfpTracks = trackAssoc.at(pfp.key());    
    std::vector<unsigned int> trkHitsTPCs_Plane0,trkHitsTPCs_Plane1,trkHitsTPCs_Plane2;
    std::vector<unsigned int> trkHitsTPCcounts_Plane0,trkHitsTPCcounts_Plane1,trkHitsTPCcounts_Plane2;
    if(!pfpTracks.empty()){
      int iTrk(0);
      //Loop over the tracks associated with each PFP
      for(const art::Ptr<recob::Track> &trk:pfpTracks){
        std::cout << "-- TRK " << iTrk << " ----------" << std::endl;
	iTrk++;

        std::vector<unsigned> hitsTpcId( fGeom->Nplanes() );
        //Loop over the planes 
        for(size_t i_plane = 0; i_plane<fGeom->Nplanes(); i_plane++) {
          std::cout << "-- PLANE " << i_plane << " ------" << std::endl;
          std::vector<const recob::Hit*> trackHits;
          auto recoTracks = e.getValidHandle<std::vector<recob::Track> >(fTrackLabel);
          art::FindManyP<recob::Hit> findHits(recoTracks,e,fTrackLabel);
          std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(trk->ID());
  
	  //Exclude hits that are not in this plane (sanity check) 
          for(const art::Ptr<recob::Hit> hit : inputHits){
            unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
            if( thePlane != i_plane ) continue;
            trackHits.push_back(hit.get());
          }

	  //Define TPC ID vector which will have an entry for every hit
          std::vector<int> plane_hits_tpcid( trackHits.size() );
         
          //Loop over hits and save TPC ID for each of them
          for(unsigned int i_hit = 0; i_hit<trackHits.size(); i_hit++ ){
            plane_hits_tpcid[i_hit] = int(trackHits[i_hit]->WireID().TPC);
            //std::cout << "plane " << i_plane << " WireID().TPC = " << trackHits[i_hit]->WireID().TPC << std::endl;
          }

	  //Loop over TPC vector and fill two new vectors for this plane:
	  //one for the non-repeating TPC IDs for this plane, and one of equal size for the number of hits in each
	  std::sort(plane_hits_tpcid.begin(), plane_hits_tpcid.end()); 
	  int occurrences(0), tpcID(-1);
          for(auto it : plane_hits_tpcid){
            if(it!=tpcID){
	      tpcID=it;
	      occurrences = std::count(plane_hits_tpcid.begin(), plane_hits_tpcid.end(), it) -1;
              std::cout << "iterator = " << it << " occurrences = " << occurrences << std::endl; 
              if(i_plane==0)trkHitsTPCs_Plane0.push_back(it);
              else if(i_plane==1)trkHitsTPCs_Plane1.push_back(it);
              else if(i_plane==2)trkHitsTPCs_Plane2.push_back(it);
              if(i_plane==0)trkHitsTPCcounts_Plane0.push_back(occurrences);
              else if(i_plane==1)trkHitsTPCcounts_Plane1.push_back(occurrences);
              else if(i_plane==2)trkHitsTPCcounts_Plane2.push_back(occurrences);
	    }
	  }	
        }
	
	//Fill tree branches with the vectors filled above, for each plane
        fTrackHitsTPCs_Plane0->push_back(trkHitsTPCs_Plane0);
        fTrackHitsTPCs_Plane1->push_back(trkHitsTPCs_Plane1);
        fTrackHitsTPCs_Plane2->push_back(trkHitsTPCs_Plane2);
        fTrackHitsTPCcounts_Plane0->push_back(trkHitsTPCcounts_Plane0);
        fTrackHitsTPCcounts_Plane1->push_back(trkHitsTPCcounts_Plane1);
        fTrackHitsTPCcounts_Plane2->push_back(trkHitsTPCcounts_Plane2);
	
	//Draw histograms for each plane: number of TPCs per track, same with a 5 and 10 hits threshold for each TPC
	fTrackCRPsize_Plane0->Fill(trkHitsTPCs_Plane0.size());
	fTrackCRPsize_Plane1->Fill(trkHitsTPCs_Plane1.size());
	fTrackCRPsize_Plane2->Fill(trkHitsTPCs_Plane2.size());
	Bool_t plane0_5hits(0), plane0_10hits(0), plane1_5hits(0), plane1_10hits(0), plane2_5hits(0), plane2_10hits(0);
	for(auto i: trkHitsTPCcounts_Plane0){if(i<5) plane0_5hits=false; break;}
	for(auto i: trkHitsTPCcounts_Plane0){if(i<10) plane0_5hits=false; break;}
	for(auto i: trkHitsTPCcounts_Plane1){if(i<5) plane1_5hits=false; break;}
	for(auto i: trkHitsTPCcounts_Plane1){if(i<10) plane1_10hits=false; break;}
	for(auto i: trkHitsTPCcounts_Plane2){if(i<5) plane2_5hits=false; break;}
	for(auto i: trkHitsTPCcounts_Plane2){if(i<10) plane2_10hits=false; break;}
	if(plane0_5hits)fTrackCRPsize_5hits_Plane0->Fill(trkHitsTPCs_Plane0.size());
	if(plane0_10hits)fTrackCRPsize_10hits_Plane0->Fill(trkHitsTPCs_Plane0.size());
	if(plane1_5hits)fTrackCRPsize_5hits_Plane1->Fill(trkHitsTPCs_Plane1.size());
	if(plane1_10hits)fTrackCRPsize_10hits_Plane1->Fill(trkHitsTPCs_Plane1.size());
	if(plane2_5hits)fTrackCRPsize_5hits_Plane2->Fill(trkHitsTPCs_Plane2.size());
	if(plane2_10hits)fTrackCRPsize_10hits_Plane2->Fill(trkHitsTPCs_Plane2.size());
        
        fTrackCRPs->Fill();
	
        trkHitsTPCs_Plane0.clear();
        trkHitsTPCs_Plane1.clear();
        trkHitsTPCs_Plane2.clear();
	trkHitsTPCcounts_Plane0.clear();
	trkHitsTPCcounts_Plane1.clear();
	trkHitsTPCcounts_Plane2.clear();
      }
    }
  }

  fTree->Fill();
}

void test::trackMetrics::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("trackHitsTPCs_Plane0",&fTrackHitsTPCs_Plane0);
  fTree->Branch("trackHitsTPCs_Plane1",&fTrackHitsTPCs_Plane1);
  fTree->Branch("trackHitsTPCs_Plane2",&fTrackHitsTPCs_Plane2);
  fTree->Branch("trackHitsTPCcounts_Plane0",&fTrackHitsTPCcounts_Plane0);
  fTree->Branch("trackHitsTPCcounts_Plane1",&fTrackHitsTPCcounts_Plane1);
  fTree->Branch("trackHitsTPCcounts_Plane2",&fTrackHitsTPCcounts_Plane2);

  fTrackCRPsize_Plane0 = tfs->make<TH1D>("TrackTrackCRPsize_Plane0","Number of Plane0 TPCs for hits in this track",5,0,5);
  fTrackCRPsize_5hits_Plane0 = tfs->make<TH1D>("TrackTrackCRPsize_5hits_Plane0","Number of Plane0 TPCs for hits in this track (5 hits/TPC min)",5,0,5);
  fTrackCRPsize_5hits_Plane0 = tfs->make<TH1D>("TrackTrackCRPsize_10hits_Plane0","Number of Plane0 TPCs for hits in this track (10 hits/TPC min)",5,0,5);
  fTrackCRPsize_Plane1 = tfs->make<TH1D>("TrackTrackCRPsize_Plane1","Number of Plane1 TPCs for hits in this track",5,0,5);
  fTrackCRPsize_5hits_Plane1 = tfs->make<TH1D>("TrackTrackCRPsize_5hits_Plane1","Number of Plane1 TPCs for hits in this track (5 hits/TPC min)",5,0,5);
  fTrackCRPsize_5hits_Plane1 = tfs->make<TH1D>("TrackTrackCRPsize_10hits_Plane1","Number of Plane1 TPCs for hits in this track (10 hits/TPC min)",5,0,5);
  fTrackCRPsize_Plane2 = tfs->make<TH1D>("TrackTrackCRPsize_Plane2","Number of Plane2 TPCs for hits in this track",5,0,5);
  fTrackCRPsize_5hits_Plane2 = tfs->make<TH1D>("TrackTrackCRPsize_5hits_Plane2","Number of Plane2 TPCs for hits in this track (5 hits/TPC min)",5,0,5);
  fTrackCRPsize_5hits_Plane2 = tfs->make<TH1D>("TrackTrackCRPsize_10hits_Plane2","Number of Plane2 TPCs for hits in this track (10 hits/TPC min)",5,0,5);

}

void test::trackMetrics::endJob()
{
}

DEFINE_ART_MODULE(test::trackMetrics)
