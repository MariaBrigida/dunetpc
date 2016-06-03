////////////////////////////////////////////////////////////////////////
// Class:       CheckGeometry
// Module Type: analyzer
// File:        CheckGeometry_module.cc
//
// Generated at Tue Jan  6 22:27:12 2015 by Tingjun Yang using artmod
// from cetpkgsupport v1_07_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TH2F.h"
#include "TLine.h"

#include <iostream>

constexpr unsigned short kMaxAuxDets = 100;
constexpr unsigned short kMaxTkIDs = 100;

namespace dune {
  class CheckGeometry;
}

class dune::CheckGeometry : public art::EDAnalyzer {
public:
  explicit CheckGeometry(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckGeometry(CheckGeometry const &) = delete;
  CheckGeometry(CheckGeometry &&) = delete;
  CheckGeometry & operator = (CheckGeometry const &) = delete;
  CheckGeometry & operator = (CheckGeometry &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.

};


dune::CheckGeometry::CheckGeometry(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  reconfigure(p);
}

void dune::CheckGeometry::analyze(art::Event const & evt)
{

  TCanvas *can = new TCanvas("c1","c1");
  can->cd();
  std::vector<TBox*> TPCBox;
  std::vector<TLine*> Wires;

  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  
  int nwires = 0;
  int nwires_tpc[8];
  for (int i = 0; i<8; ++i) nwires_tpc[i] = 0;
  art::ServiceHandle<geo::Geometry> geo;
  for (size_t t = 0; t<geo->NTPC(); ++t){
    //if (t%2==0) continue;
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(t);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-tpc.ActiveHalfWidth())
      minx = world[0]-tpc.ActiveHalfWidth();
    if (maxx<world[0]+tpc.ActiveHalfWidth())
      maxx = world[0]+tpc.ActiveHalfWidth();
    if (miny>world[1]-tpc.ActiveHalfHeight())
      miny = world[1]-tpc.ActiveHalfHeight();
    if (maxy<world[1]+tpc.ActiveHalfHeight())
      maxy = world[1]+tpc.ActiveHalfHeight();
    if (minz>world[2]-tpc.ActiveLength()/2.)
      minz = world[2]-tpc.ActiveLength()/2.;
    if (maxz<world[2]+tpc.ActiveLength()/2.)
      maxz = world[2]+tpc.ActiveLength()/2.;
    std::cout<<t<<" "<<world[0]-tpc.ActiveHalfWidth()
	     <<" "<<world[0]+tpc.ActiveHalfWidth()
	     <<" "<<world[1]-tpc.ActiveHalfHeight()
	     <<" "<<world[1]+tpc.ActiveHalfHeight()
	     <<" "<<world[2]-tpc.ActiveLength()/2.
	     <<" "<<world[2]+tpc.ActiveLength()/2.<<std::endl;
 
    TPCBox.push_back(new TBox(world[2]-tpc.ActiveLength()/2.,
			      world[1]-tpc.ActiveHalfHeight(),
			      world[2]+tpc.ActiveLength()/2.,
			      world[1]+tpc.ActiveHalfHeight()));
    TPCBox.back()->SetFillStyle(0);
    TPCBox.back()->SetLineStyle(2);
    TPCBox.back()->SetLineWidth(2);
    TPCBox.back()->SetLineColor(16);
    
    for (size_t p = 0; p<geo->Nplanes(t);++p){
      for (size_t w = 0; w<geo->Nwires(p,t); ++w){
	++nwires;
	++nwires_tpc[t];
	double xyz0[3];
	double xyz1[3];
	unsigned int c = 0;
//	if ((t==7&&p==0&&w==192)||
//	    (t==7&&p==1&&w==112)||
//	    (t==7&&p==2&&w==0)){
	if (true){
	  geo->WireEndPoints(c,t,p,w,xyz0,xyz1);
	  Wires.push_back(new TLine(xyz0[2],xyz0[1],xyz1[2],xyz1[1]));
	}
	//std::cout<<t<<" "<<p<<" "<<w<<" "<<xyz0[0]<<" "<<xyz0[1]<<" "<<xyz0[2]<<std::endl;
      }
    }
  }

  TH2F *frame = new TH2F("frame",";z (cm);y (cm)",100,minz,maxz,100,miny,maxy);
  frame->SetStats(0);
  frame->Draw();
  for (auto box: TPCBox) box->Draw();
  for (auto wire: Wires) wire->Draw();
  can->Print("wires.pdf");
  std::cout<<"N wires = "<<nwires<<std::endl;
  for (int i = 0; i<8; ++i){
    std::cout<<"TPC "<<i<<" has "<<nwires_tpc[i]<<" wires"<<std::endl;
  }
}

void dune::CheckGeometry::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(dune::CheckGeometry)
