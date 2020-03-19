////////////////////////////////////////////////////////////////////////
/// \file    CTPResult.cxx
/// \brief   Class storing the result from the convolutional track PID
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>

#include "dune/TrackPID/CTPResult.h"

namespace ctp
{


  CTPResult::CTPResult(){
    fMuonScore = 0.;
    fPionScore = 0.;
    fProtonScore = 0.;
  }

  CTPResult::CTPResult(const std::vector<float> &vals){
    fMuonScore = 0.;
    fPionScore = 0.;
    fProtonScore = 0.;
    if(vals.size() != 3){
      std::cout << "CTPResult Error: there should be three input values" << std::endl;
    }
    else{
      fMuonScore = vals.at(0);
      fPionScore = vals.at(1);
      fProtonScore = vals.at(2);
    }
  }

  CTPResult::~CTPResult(){

  }

  bool CTPResult::IsValid() const{
    return (fMuonScore > 1e-5) && (fPionScore > 1e-5) && (fProtonScore > 1e-5);
  }

}

