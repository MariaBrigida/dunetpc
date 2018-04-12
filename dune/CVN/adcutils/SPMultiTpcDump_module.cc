//////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SPMultiTpcDump
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan (Dorota.Stefan@cern.ch), May 2016
//
// Produces "global image" of events in multi-tpc detector, could be FD workspace or ProtoDUNE (SP).
// Projections of ADC are saved to 1-byte / pixel (TH1C) format, int-size (TH2I) projection of PDG
// truth and float-size (TH2F) projection of true charge deposits can be saved as well. We use this
// to dump deconvoluted ADC images for preparation of various classifiers and CNN based pattern reco.
//
// Moved from larreco/RecoAlg/ImagePatternAlgs since it became specific to DUNE geometries.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SPMultiTpcDump_Module
#define SPMultiTpcDump_Module

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

// C++ Includes
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH2C.h" // ADC map
#include "TH2I.h" // PDG+vertex info map
#include "TH2F.h" // deposit map

namespace nnet
{

  struct NUVTX
  {
  	int interaction;
  	int nupdg;
  	int cryo;
  	int tpc;
  	TVector2 position[3];
  	TVector3 vtx3d;
  };

  class EventImageData // full image for one plane
  {
  public:
    EventImageData(size_t w, size_t d, bool saveDep) :
        fVtxX(-9999), fVtxY(-9999),
        fProjX(-9999), fProjY(-9999),
        fSaveDep(saveDep)
    {
        fAdc.resize(w, std::vector<float>(d, 0));
        if (saveDep) { fDeposit.resize(w, std::vector<float>(d, 0)); }
        fPdg.resize(w, std::vector<int>(d, 0));
    }

    void addTpc(const TrainingDataAlg & dataAlg, size_t gw, bool flipw, size_t gd, bool flipd);
    bool findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const;

    const std::vector< std::vector<float> > & adcData(void) const { return fAdc; }
    const std::vector<float> & wireAdc(size_t widx) const { return fAdc[widx]; }

    const std::vector< std::vector<float> > & depData(void) const { return fDeposit; }
    const std::vector<float> & wireDep(size_t widx) const { return fDeposit[widx]; }

    const std::vector< std::vector<int> > & pdgData(void) const { return fPdg; }
    const std::vector<int> & wirePdg(size_t widx) const { return fPdg[widx]; }

    void setProjXY(const TrainingDataAlg & dataAlg, float x, float y, size_t gw, bool flipw, size_t gd, bool flipd);
    float getProjX(void) const { return fProjX; }
    float getProjY(void) const { return fProjY; }

    int getVtxX(void) const { return fVtxX; }
    int getVtxY(void) const { return fVtxY; }

  private:
    std::vector< std::vector<float> > fAdc, fDeposit;
    std::vector< std::vector<int> > fPdg;
    int fVtxX, fVtxY;
    float fProjX, fProjY;
    bool fSaveDep;
  };

  void EventImageData::setProjXY(const TrainingDataAlg & dataAlg, float x, float y, size_t gw, bool flipw, size_t gd, bool flipd)
  {
    if (flipw) { fProjX = gw + dataAlg.NWires() - x - 1; }
    else { fProjX = gw + x; }

    if (flipd) { fProjY = gd + dataAlg.NScaledDrifts() - y/dataAlg.DriftWindow() - 1; }
    else { fProjY = gd + y/dataAlg.DriftWindow(); }
  }

  void EventImageData::addTpc(const TrainingDataAlg & dataAlg, size_t gw, bool flipw, size_t gd, bool flipd)
  {
    float zero = dataAlg.ZeroLevel();
    for (size_t w = 0; w < dataAlg.NWires(); ++w)
    {
        size_t drift_size = dataAlg.NScaledDrifts();
        std::vector<float> & dstAdc = fAdc[gw + w];
        std::vector<float> & dstDep = fDeposit[gw + w];
        std::vector<int> & dstPdg = fPdg[gw + w];
        const float* srcAdc = 0;
        const float* srcDep = 0;
        const int* srcPdg = 0;
        if (flipw)
        {
            srcAdc = dataAlg.wireData(dataAlg.NWires() - w - 1).data();
            srcDep = dataAlg.wireEdep(dataAlg.NWires() - w - 1).data();
            srcPdg = dataAlg.wirePdg(dataAlg.NWires() - w - 1).data();
        }
        else
        {
            srcAdc = dataAlg.wireData(w).data();
            srcDep = dataAlg.wireEdep(w).data();
            srcPdg = dataAlg.wirePdg(w).data();
        }
        if (flipd)
        {
            for (size_t d = 0; d < drift_size; ++d) { dstAdc[gd + d] += srcAdc[drift_size - d - 1] - zero; }
            if (fSaveDep) { for (size_t d = 0; d < drift_size; ++d) { dstDep[gd + d] += srcDep[drift_size - d - 1]; } }
            for (size_t d = 0; d < drift_size; ++d)
            {
                int code = srcPdg[drift_size - d - 1];
                int best_pdg = code & nnet::TrainingDataAlg::kPdgMask;
                int vtx_flags = (dstPdg[gd + d] | code) & nnet::TrainingDataAlg::kVtxMask;
                dstPdg[gd + d] = vtx_flags | best_pdg; // now just overwrite pdg and keep all vtx flags
                
                if (code & nnet::TrainingDataAlg::kNuPri) { fVtxX = gw + w; fVtxY = gd + d; } // dstAdc[gd + d] = 127; }
            }
        }
        else
        {
            for (size_t d = 0; d < drift_size; ++d) { dstAdc[gd + d] += srcAdc[d] - zero; }
            if (fSaveDep) { for (size_t d = 0; d < drift_size; ++d) { dstDep[gd + d] += srcDep[d]; } }
            for (size_t d = 0; d < drift_size; ++d)
            {
                int code = srcPdg[d];
                int best_pdg = code & nnet::TrainingDataAlg::kPdgMask;
                int vtx_flags = (dstPdg[gd + d] | code) & nnet::TrainingDataAlg::kVtxMask;
                dstPdg[gd + d] = vtx_flags | best_pdg; // now just overwrite pdg and keep all vtx flags

                if (code & nnet::TrainingDataAlg::kNuPri) { fVtxX = gw + w; fVtxY = gd + d; } // dstAdc[gd + d] = 127; }
            }
        }
    }
  }

  bool EventImageData::findCrop(size_t max_area_cut, unsigned int & w0, unsigned int & w1, unsigned int & d0, unsigned int & d1) const
  {
    size_t max_cut = max_area_cut / 4;
    float adcThr = 10;

    w0 = 0;
    size_t cut = 0;
    while (w0 < fAdc.size())
    {
        for (auto const d : fAdc[w0]) { if (d > adcThr) cut++; }
        if (cut < max_cut) w0++;
        else break;
    }
    w1 = fAdc.size() - 1;
    cut = 0;
    while (w1 > w0)
    {
        for (auto const d : fAdc[w1]) { if (d > adcThr) cut++; }
        if (cut < max_cut) w1--;
        else break;
    }
    w1++;

    d0 = 0;
    cut = 0;
    while (d0 < fAdc.front().size())
    {
        for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d0] > adcThr) cut++; }
        if (cut < max_cut) d0++;
        else break;
    }
    d1 = fAdc.front().size() - 1;
    cut = 0;
    while (d1 > d0)
    {
        for (size_t i = w0; i < w1; ++i) { if (fAdc[i][d1] > adcThr) cut++; }
        if (cut < max_cut) d1--;
        else break;
    }
    d1++;

    unsigned int margin = 32;
    if ((w1 - w0 > 8) && (d1 - d0 > 8))
    {
        if (w0 < margin) w0 = 0;
        else w0 -= margin;

        if (w1 > fAdc.size() - margin) w1 = fAdc.size();
        else w1 += margin;
        
        if (d0 < margin) d0 = 0;
        else d0 -= margin;
        
        if (d1 > fAdc.front().size() - margin) d1 = fAdc.front().size();
        else d1 += margin;
        
        return true;
    }
    else return false;
  }


  class SPMultiTpcDump : public art::EDAnalyzer
  {
  public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::TrainingDataAlg::Config> TrainingDataAlg { Name("TrainingDataAlg") };
		fhicl::Atom<art::InputTag> GenModuleLabel { Name("GenModuleLabel"), Comment("Neutrino generator label.") };
		fhicl::Atom<double> FidVolCut { Name("FidVolCut"), Comment("Take events with vertex inside this cut on volume.") };
		fhicl::Atom<bool> SaveDepositMap { Name("SaveDepositMap"), Comment("Save projections of the true energy depositions.") };
		fhicl::Atom<bool> SavePdgMap { Name("SavePdgMap"), Comment("Save vertex info and PDG codes map.") };
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit SPMultiTpcDump(Parameters const& config);
    
    void beginJob() override;

    void analyze(const art::Event& event) override;

  private:
  	void ResetVars();
  
  	bool prepareEv(const art::Event& event);
  	bool InsideFidVol(TVector3 const & pvtx) const;
  	void CorrOffset(TVector3& vec, const simb::MCParticle& particle);


	TVector2 GetProjVtx(TVector3 const & vtx3d, const size_t cryo, const size_t tpc, const size_t plane) const;

    nnet::TrainingDataAlg fTrainingDataAlg;
    art::InputTag fGenieGenLabel;

	TTree *fTree;
	TTree *fTree2D;
	int fEvent;     ///< number of the event being processed
	int fRun;       ///< number of the run being processed
	int fSubRun;    ///< number of the sub-run being processed
	
	bool fSaveDepositMap, fSavePdgMap;
		
	double fFidVolCut;
	
	NUVTX fPointid;
	int fCryo, fTpc, fPlane;
	int fPdg, fInteraction;
	int fPixX, fPixY;
	float fPosX, fPosY;

	geo::GeometryCore const* fGeometry;
	detinfo::DetectorProperties const* fDetProp;
  };

  //-----------------------------------------------------------------------
  SPMultiTpcDump::SPMultiTpcDump(SPMultiTpcDump::Parameters const& config) : art::EDAnalyzer(config),
	fTrainingDataAlg(config().TrainingDataAlg()),
	fGenieGenLabel(config().GenModuleLabel()),
	fSaveDepositMap(config().SaveDepositMap()),
	fSavePdgMap(config().SavePdgMap()),
	fFidVolCut(config().FidVolCut())
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  
  //-----------------------------------------------------------------------
  void SPMultiTpcDump::beginJob()
  {
		art::ServiceHandle<art::TFileService> tfs;

		fTree = tfs->make<TTree>("event","Event level info");
		fTree->Branch("fRun", &fRun, "fRun/I");
		fTree->Branch("fEvent", &fEvent, "fEvent/I");
		fTree->Branch("fCryo", &fCryo, "fCryo/I");
		fTree->Branch("fTpc", &fTpc, "fTpc/I");
		fTree->Branch("fPdg", &fPdg, "fPdg/I");
		fTree->Branch("fInteraction", &fInteraction, "fInteraction/I");
		
		fTree2D = tfs->make<TTree>("plane","Vertex 2D info");
		fTree2D->Branch("fRun", &fRun, "fRun/I");
		fTree2D->Branch("fEvent", &fEvent, "fEvent/I");
		fTree2D->Branch("fPlane", &fPlane, "fPlane/I");
		fTree2D->Branch("fPixX", &fPixX, "fPixX/I");
		fTree2D->Branch("fPixY", &fPixY, "fPixY/I");
		fTree2D->Branch("fPosX", &fPosX, "fPosX/F");
		fTree2D->Branch("fPosY", &fPosY, "fPosY/F");
  }
  
  //-----------------------------------------------------------------------
  bool SPMultiTpcDump::prepareEv(const art::Event& event)
  {
	art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
	if (!event.getByLabel(fGenieGenLabel, mctruthHandle)) { return false; }

	for (auto const & mc : (*mctruthHandle))
	{
		if (mc.Origin() == simb::kBeamNeutrino)
		{
			const TLorentzVector& pvtx = mc.GetNeutrino().Nu().Position();
			TVector3 vtx(pvtx.X(), pvtx.Y(), pvtx.Z());

			CorrOffset(vtx, mc.GetNeutrino().Nu());

			fPointid.nupdg = mc.GetNeutrino().Nu().PdgCode();
			fPointid.interaction = mc.GetNeutrino().CCNC();
			fPointid.vtx3d = vtx;

			if (InsideFidVol(vtx)) 
			{	
				double v[3] = {vtx.X(), vtx.Y(), vtx.Z()}; 
				unsigned int cryo = fGeometry->FindCryostatAtPosition(v);
				unsigned int tpc = fGeometry->FindTPCAtPosition(v).TPC;
				
				for (size_t i = 0; i < 3; ++i)
				{
					if (fGeometry->TPC(tpc, cryo).HasPlane(i))
					{ fPointid.position[i] = GetProjVtx(vtx, cryo, tpc, i); }
					else
					{ fPointid.position[i] = TVector2(0, 0); }
				}

				fPointid.cryo = cryo;
				fPointid.tpc = tpc;
				
				fCryo = cryo;
				fTpc = tpc;
				
				return true;	
			}
			else
			{
				fPointid.cryo = -1;
				fPointid.tpc = -1;
				fCryo = -1; fTpc = -1;
			    return false;
			}
		}
	}
	return false;
  }  

  //-----------------------------------------------------------------------
  void SPMultiTpcDump::analyze(const art::Event& event) 
  {
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();
    ResetVars();
    
    prepareEv(event);

	std::ostringstream os;
	os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;
	std::cout << "analyze " << os.str() << std::endl;

    size_t cryo = 0;
    size_t plane_align0[3] = { 1, 0, 2 }; // offsets to align signals in drift direction between planes
    size_t plane_align1[3] = { 2, 1, 0 };
    size_t max_align = 0;
    for (size_t i = 0; i < 3; ++i)
    {
        if (plane_align0[i] > max_align) max_align = plane_align0[i];
        if (plane_align1[i] > max_align) max_align = plane_align1[i];
    }

    size_t weff[3] = { 404, 404, 485 }; // effective offset in wires between consecutive tpc's, incl. dead spaces, per plane
    size_t tpcs[3] = { 2, 2, 6 };       // FD workspace dimension in #tpc's
    //size_t tpcs[3] = { 4, 1, 3 };     // the same for ProtoDUNE
    size_t apa_gap = 29;                // APA thickness (dead space in drift direction) in pixels

    bool goodEvent = false;
    unsigned int gd0 = 0, gd1 = 0;
	for (int p = 2; p >= 0; --p) // loop over planes and make global images
	{
	    size_t drift_size = fDetProp->NumberTimeSamples() / fTrainingDataAlg.DriftWindow(); // tpc size in drift direction
	    size_t wire_size = fGeometry->Nwires(p, 0, cryo);                                   // tpc size in wire direction
	    size_t max_offset = (p < 2) ? 752 : 0;                                              // top-bottom projection max offset

        size_t ntotw = (tpcs[2] - 1) * weff[p] + wire_size + tpcs[1] * max_offset;          // global image size in wire direction
        size_t ntotd = tpcs[0] * drift_size + (apa_gap + max_align) * (tpcs[0] / 2);        // global image size in drift direction

	    std::cout << "Plane: " << p << ", wires: " << ntotw << " drift: " << ntotd << std::endl;
	    EventImageData fullimg(ntotw, ntotd, fSaveDepositMap);

        size_t a, maxArea = 0;
	    for (size_t t = 0; t < fGeometry->NTPC(cryo); ++t) // loop over tpc's
	    {
            size_t tpc_z = t / (tpcs[0] * tpcs[1]);
            size_t tpc_y = (t / tpcs[0]) % tpcs[1];
            size_t tpc_x = t % tpcs[0];

            int dir = fGeometry->TPC(t, cryo).DetectDriftDirection();
            bool flip_d = (dir > 0); // need the flip in drift direction?

	        bool flip_w = false;     // need the flip in wire direction?
	        size_t eff_p = p;        // global plane
	        int offset = 0;
            int p_align = 0;

            if (p < 2) // no wire flip nor top-bottom alignments in Z (collection) plane
            {
	            if ((t % 4 == 0) || (t % 4 == 2))
	            {
	                eff_p = (p == 1) ? 0 : 1;     // swap U and V
	            }

                if ((t % 4 == 0) || (t % 4 == 1)) // bottom tpc's
                {
                    flip_w = (dir > 0) ? (eff_p == 1) : (eff_p == 0);
                }
                else                              // top tpc's
                {
                    flip_w = (dir < 0) ? (eff_p == 1) : (eff_p == 0);
                }

                offset = (p == 0) ? 48 : 752;     // top-bottopm offsets
            }
            
            if ((t % 4 == 0) || (t % 4 == 2))
            {
                p_align = plane_align0[p];
            }
            else
            {
                p_align = plane_align1[p];
            }

            size_t gw = tpc_z * weff[p] + tpc_y * offset;                          // global wire
            size_t gd = tpc_x * drift_size + apa_gap * (1 + tpc_x) / 2 + p_align;  // global drift

		    fTrainingDataAlg.setEventData(event, eff_p, t, cryo); // prepare ADC and truth projections for 1 plane in 1 tpc
		    std::cout << "   TPC: " << t << " wires: " << fTrainingDataAlg.NWires() << std::endl;
		    a = fTrainingDataAlg.getAdcArea();
		    if (a > 150)
		    {
		        fullimg.addTpc(fTrainingDataAlg, gw, flip_w, gd, flip_d);
		        if (a > maxArea)
		        {
                    TVector2 vtxProj = GetProjVtx(fPointid.vtx3d, cryo, t, p);
                    fullimg.setProjXY(fTrainingDataAlg, vtxProj.X(), vtxProj.Y(), gw, flip_w, gd, flip_d);
		            a = maxArea;
		        }
		    }
		}

        std::cout << "   find crop..." << std::endl;
        unsigned int w0, w1, d0, d1;
		if (fullimg.findCrop(80, w0, w1, d0, d1))
   		{
   		    if (goodEvent)
   		    {
   		        d0 = gd0; d1 = gd1;
   		    }
   		    else
   		    {
   		        goodEvent = true;
   		        gd0 = d0; gd1 = d1;
   		    }
   			std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
   		}
	    else { std::cout << "   skip empty event" << std::endl; break; }

		std::ostringstream ss1;
   		ss1 << os.str() << "_plane_" << p; // TH2's name

		art::ServiceHandle<art::TFileService> tfs;
   		TH2C* rawHist = tfs->make<TH2C>((ss1.str() + "_raw").c_str(), "ADC",
                    (int)(w1 - w0), (double)w0, (double)w1, (int)(d1 - d0), (double)d0, (double)d1);
   		
   		float zero = fTrainingDataAlg.ZeroLevel();
        for (size_t w = w0; w < w1; ++w)
        {
            auto const & raw = fullimg.wireAdc(w);
            for (size_t d = d0; d < d1; ++d)
            {
                rawHist->Fill(w, d, (char)(raw[d] + zero));
            }
   		}

   		if (fSaveDepositMap)
   		{
       		TH2F* depHist = tfs->make<TH2F>((ss1.str() + "_deposit").c_str(), "Deposit",
                    (int)(w1 - w0), (double)w0, (double)w1, (int)(d1 - d0), (double)d0, (double)d1);
            for (size_t w = w0; w < w1; ++w)
            {
                auto const & edep = fullimg.wireDep(w);
                for (size_t d = d0; d < d1; ++d) { depHist->Fill(w, d, edep[d]); }
       		}
       	}

       	if (fSavePdgMap)
       	{
   		TH2I* pdgHist = tfs->make<TH2I>((ss1.str() + "_pdg").c_str(), "PDG",
                    (int)(w1 - w0), (double)w0, (double)w1, (int)(d1 - d0), (double)d0, (double)d1);
            for (size_t w = w0; w < w1; ++w)
            {
                auto const & pdg = fullimg.wirePdg(w);
                for (size_t d = d0; d < d1; ++d) { pdgHist->Fill(w, d, pdg[d]); }
       		}
   		}

        if (goodEvent)
        {
    		fPlane = p;

  	    	if (fullimg.getVtxX() > -9999) fPixX = fullimg.getVtxX() - w0;
	    	if (fullimg.getVtxY() > -9999) fPixY = fullimg.getVtxY() - d0;

      		if (fullimg.getProjX() > -9999) fPosX = fullimg.getProjX() - w0;
    		if (fullimg.getProjY() > -9999) fPosY = fullimg.getProjY() - d0;

            std::cout << " *** plane:" << p << std::endl;
            std::cout << " ***   w0:" << w0 << ", w1:" << w1 << std::endl;
            std::cout << " ***   d0:" << d0 << ", d1:" << d1 << std::endl;
            std::cout << " ***  pix *** x:" << fPixX << ", y:" << fPixY << std::endl;
            std::cout << " *** proj *** x:" << fPosX << ", y:" << fPosY << std::endl;
            std::cout << " *** zero:" << zero << std::endl;

    		fTree2D->Fill();
    	}
	}
	if (goodEvent)
	{
   		fPdg = fPointid.nupdg;
		fInteraction = fPointid.interaction;
	    fTree->Fill();
	}
  }
  
  //-----------------------------------------------------------------------
  void SPMultiTpcDump::CorrOffset(TVector3& vec, const simb::MCParticle& particle)
  {
  	double vtx[3] = {vec.X(), vec.Y(), vec.Z()};
  	geo::TPCID tpcid = fGeometry->FindTPCAtPosition(vtx);

	if (tpcid.isValid)
	  {
	  	float corrt0x = particle.T() * 1.e-3 * fDetProp->DriftVelocity();
	  	if (fGeometry->TPC(tpcid).DetectDriftDirection() == 1) { corrt0x = corrt0x*(-1); }
	  	
	  	vtx[0] = vec.X() + corrt0x;
	  }

	vec.SetX(vtx[0]);
  }

  //-----------------------------------------------------------------------
  bool SPMultiTpcDump::InsideFidVol(TVector3 const & pvtx) const
  {
		double vtx[3] = {pvtx.X(), pvtx.Y(), pvtx.Z()};
		bool inside = false;

		if (!fGeometry->FindTPCAtPosition(vtx).isValid) return false;

		geo::TPCID idtpc = fGeometry->FindTPCAtPosition(vtx);

		if (fGeometry->HasTPC(idtpc))
		{		
			const geo::TPCGeo& tpcgeo = fGeometry->GetElement(idtpc);
			double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
			double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
			double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

			//x
			double dista = fabs(minx - pvtx.X());
			double distb = fabs(pvtx.X() - maxx); 

			if ((pvtx.X() > minx) && (pvtx.X() < maxx) &&
		 	(dista > fFidVolCut) && (distb > fFidVolCut))
			{ 
				inside = true;
			}
			else { inside = false; }

			//y
			dista = fabs(maxy - pvtx.Y());
			distb = fabs(pvtx.Y() - miny);
			if (inside && (pvtx.Y() > miny) && (pvtx.Y() < maxy) &&
		 		(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
			else inside = false;

			//z
			dista = fabs(maxz - pvtx.Z());
			distb = fabs(pvtx.Z() - minz);
			if (inside && (pvtx.Z() > minz) && (pvtx.Z() < maxz) &&
		 		(dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
			else inside = false;
		}
		
		return inside;
  }
  
  //-----------------------------------------------------------------------
  TVector2 SPMultiTpcDump::GetProjVtx(TVector3 const & vtx3d, const size_t cryo, const size_t tpc, const size_t plane) const
  {

		TVector2 vtx2d = pma::GetProjectionToPlane(vtx3d, plane, tpc, cryo);
		TVector2 vtxwd = pma::CmToWireDrift(vtx2d.X(), vtx2d.Y(), plane, tpc, cryo);
	
		return vtxwd;
  }

	void SPMultiTpcDump::ResetVars()
	{
		fCryo = 0;
		fTpc = 0;
		fPlane = 0;
		fPdg = 0;
		fInteraction = 0;
		fPosX = 0.0;
		fPosY = 0.0;
                fPixX = 0;
                fPixY = 0;
	}

DEFINE_ART_MODULE(SPMultiTpcDump)

}

#endif // SPMultiTpcDump_Module

