////////////////////////////////////////////////////////////////////////
/// \file    TFRegNetHandler.cxx
/// \brief   TFRegNetHandler for RegCVN modified from TFNetHandler.cxx
/// \author  Ilsoo Seong - iseong@uci.edu 
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <string>
#include "cetlib/getenv.h"

#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dune/RegCVN/art/TFRegNetHandler.h"
#include "dune/RegCVN/func/RegCVNImageUtils.h"

#include "TH2D.h"
#include "TCanvas.h"

namespace cvn
{

  TFRegNetHandler::TFRegNetHandler(const fhicl::ParameterSet& pset):
    fLibPath(cet::getenv(pset.get<std::string>("LibPath", ""))),
    //fLibPath((pset.get<std::string>("LibPath", ""))),
    fTFProtoBuf  (fLibPath+"/"+pset.get<std::string>("TFProtoBuf")),
    fInputs(pset.get<unsigned int>("NInputs")),
    fReverseViews(pset.get<std::vector<bool> >("ReverseViews"))
  {

    // Construct the TF Graph object. The empty vector {} is used since the protobuf
    // file gives the names of the output layer nodes
    mf::LogInfo("TFRegNetHandler") << "Loading network: " << fTFProtoBuf << std::endl;
    fTFGraph = tf::RegCVNGraph::create(fTFProtoBuf.c_str(),fInputs,{});
    if(!fTFGraph){
      art::Exception(art::errors::Unknown) << "Tensorflow model not found or incorrect";
    }

  }
  
  std::vector<float> TFRegNetHandler::Predict(const RegPixelMap& pm)
  {
   
    RegCVNImageUtils imageUtils;

    // Configure the image utility  
    imageUtils.SetViewReversal(fReverseViews);

    //std::cout << "Reverse views? [" << fReverseViews[0] << "," << fReverseViews[1] << "," << fReverseViews[2] << "]" << std::endl;

    ImageVectorF thisImage;
    imageUtils.ConvertPixelMapToImageVectorF(pm,thisImage);
    std::vector<ImageVectorF> vecForTF;
/*
    // Does this image look sensible?
    TCanvas *can = new TCanvas("can","",0,0,800,600);
    unsigned int fImageWires = pm.fNWire;
    unsigned int fImageTDCs = pm.fNTdc;
    TH2D* hView0 = new TH2D("hView0","",fImageWires,0,fImageWires,fImageTDCs,0,fImageTDCs);
    TH2D* hView1 = new TH2D("hView1","",fImageWires,0,fImageWires,fImageTDCs,0,fImageTDCs);
    TH2D* hView2 = new TH2D("hView2","",fImageWires,0,fImageWires,fImageTDCs,0,fImageTDCs);
    for(unsigned int w = 0; w < fImageWires; ++w){
      for(unsigned int t = 0; t < fImageTDCs; ++t){
        hView0->SetBinContent(w+1,t+1,thisImage[w][t][0]);
        hView1->SetBinContent(w+1,t+1,thisImage[w][t][1]);
        hView2->SetBinContent(w+1,t+1,thisImage[w][t][2]);
      }
    }
    hView0->Draw("colz");
    can->Print("view0a.png");
    hView1->Draw("colz");
    can->Print("view1a.png");
    hView2->Draw("colz");
    can->Print("view2a.png");
*/
    vecForTF.push_back(thisImage);

    auto cvnResults = fTFGraph->run(vecForTF, fInputs);

    //std::cout << "Number of CVN result vectors " << cvnResults.size() << " with " << cvnResults[0].size() << " categories" << std::endl;

    //std::cout << "summary: ";
    //for(auto const v : cvnResults[0]){
    //  std::cout << v << ", ";
    //}
    //std::cout << std::endl;

    return cvnResults[0];
  }
 
  std::vector<float> TFRegNetHandler::PredictNuEEnergy(const RegPixelMap& pm){
    std::vector<float> fullResults = this->Predict(pm);
    std::vector<float> nue_energy;
    nue_energy.push_back(fullResults[0]);
    return nue_energy;
  }

}

