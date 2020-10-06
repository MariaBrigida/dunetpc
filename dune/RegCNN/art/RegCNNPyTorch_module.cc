// C/C++ includes
#include <iostream>
#include <sstream>

#include <torch/script.h>
#include <torch/torch.h>

#include "cetlib/getenv.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

#include "dune/RegCNN/func/RegCNNResult.h"
#include "dune/RegCNN/func/RegPixelMap3D.h"

namespace cnn {

    class RegCNNPyTorch : public art::EDProducer {

        public:
            explicit RegCNNPyTorch(fhicl::ParameterSet const& pset);
            ~RegCNNPyTorch();

            void produce(art::Event& evt);
            void beginJob();
            void endJob();

        private:

            std::string fLibPath;
            std::string fNetwork;
            std::string fPixelMapInput;
            std::string fResultLabel;
        
            torch::jit::script::Module module;
    }; // class RegCNNPyTorch

    RegCNNPyTorch::RegCNNPyTorch(fhicl::ParameterSet const& pset):
        EDProducer(pset),
        fLibPath       (cet::getenv(pset.get<std::string>      ("LibPath", ""))),
        fNetwork       (fLibPath + "/" + pset.get<std::string> ("Network")),
        fPixelMapInput (pset.get<std::string>                  ("PixelMapInput")),
        fResultLabel   (pset.get<std::string>                  ("ResultLabel"))
    {
        produces<std::vector<cnn::RegCNNResult> >(fResultLabel);
    }

    RegCNNPyTorch::~RegCNNPyTorch() {

    }

    void RegCNNPyTorch::beginJob() {
        std::cout<<"regcnn_torch job begins ...... "<<std::endl;
        try {
            // Deserialize the ScriptModule from a file using torch::jit::load().
            module = torch::jit::load(fNetwork);
        }
        catch (const c10::Error& e) {
            std::cerr<<"error loading the model\n";
            return;
        }
        std::cout<<"loaded model "<<fNetwork<<" ... ok\n"<<std::endl;
    }

    void RegCNNPyTorch::endJob() {
    }

    void RegCNNPyTorch::produce(art::Event& evt) {
        /// Define containers for the things we're going to produce
        std::unique_ptr< std::vector<RegCNNResult> >
                                      resultCol(new std::vector<RegCNNResult>);

        /// Load in 3D pixel map for direction reco.
        art::Handle< std::vector< cnn::RegPixelMap3D > > pixelmap3DListHandle;
        std::vector< art::Ptr< cnn::RegPixelMap3D > > pixelmap3Dlist;
        if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmap3DListHandle)) {
            art::fill_ptr_vector(pixelmap3Dlist, pixelmap3DListHandle);
        }

        if (pixelmap3Dlist.size() > 0) {
            std::cout<<"have pixelmap3Dlist"<<std::endl;
            RegPixelMap3D pm = *pixelmap3Dlist[0];

            at::Tensor t_pm = torch::from_blob(pm.fPE.data(), {1,1,32,32,32});

            std::vector<torch::jit::IValue> inputs_pm;
            inputs_pm.push_back(t_pm);
            at::Tensor torchOutput = module.forward(inputs_pm).toTensor();

            std::vector<float> networkOutput;
            for (unsigned int i= 0; i< 3; ++i) {
                networkOutput.push_back(torchOutput[0][i].item<float>());
                std::cout<<torchOutput[0][i].item<float>()<<std::endl;
            }

            //std::vector<float> networkOutput(3);
            resultCol->emplace_back(networkOutput);

            //std::cout<<output.slice(/*dim=*/1, /*start=*/0, /*end=*/10) << '\n';
            //std::cout<<output[0]<<std::endl;
        }

        evt.put(std::move(resultCol), fResultLabel);
    }

    DEFINE_ART_MODULE(cnn::RegCNNPyTorch)
} // end namespace cnn
