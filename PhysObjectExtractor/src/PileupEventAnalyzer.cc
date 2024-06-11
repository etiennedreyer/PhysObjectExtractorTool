// -*- C++ -*-
//
// Package:    PileupEventAnalyzer
// Class:      PileupEventAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// classes to extract PileupSimmaryInfor information
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Provenance/interface/EventID.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class PileupEventAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PileupEventAnalyzer(const edm::ParameterSet&);
      ~PileupEventAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);

      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

	TTree *tree;
    
    // PileUp Event information
    int numPU;
    std::vector<float> PU_zpositions;
    std::vector<ULong64_t> PU_eventIDs;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

PileupEventAnalyzer::PileupEventAnalyzer(const edm::ParameterSet& iConfig)
{
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("Events", "Events");

    // Event information
   tree->Branch("numPU", &numPU);
   tree->GetBranch("numPU")->SetTitle("Number of Pileup Events");
   tree->Branch("PU_zpositions", &PU_zpositions);
   tree->GetBranch("PU_zpositions")->SetTitle("Z positions of Pileup Events");
   tree->Branch("PU_eventIDs", &PU_eventIDs);
   tree->GetBranch("PU_eventIDs")->SetTitle("Event IDs of Pileup Events");
    
}

PileupEventAnalyzer::~PileupEventAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PileupEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   // Event information
   numPU = 0;
   PU_zpositions.clear();
   PU_eventIDs.clear();

   Handle<std::vector<PileupSummaryInfo>> pileupInfo;
   iEvent.getByLabel("addPileupInfo", pileupInfo);

   for (std::vector<PileupSummaryInfo>::const_iterator itPUInfo = pileupInfo->begin(); itPUInfo != pileupInfo->end(); ++itPUInfo) {
      numPU++;
      std::vector<float> zpositions = itPUInfo->getPU_zpositions();
      std::vector<EventID> eventIDs = itPUInfo->getPU_EventID();
      // std::cout << "Number of Pileup Events: " << itPUInfo->getPU_NumInteractions() << std::endl;
      // std::cout << "Number of Z positions: " << zpositions.size() << std::endl;
      // std::cout << "Number of Event IDs: " << eventIDs.size() << std::endl;
      for (std::vector<float>::const_iterator itZ = zpositions.begin(); itZ != zpositions.end(); ++itZ) {
         PU_zpositions.push_back(*itZ);
      }
      for (std::vector<EventID>::const_iterator itID = eventIDs.begin(); itID != eventIDs.end(); ++itID) {
         PU_eventIDs.push_back(itID->event());
         // std::cout << "Event ID: " << itID->event() << std::endl;
      }
   }

   tree->Fill();
   return;
}



// ------------ method called once each job just before starting event loop  ------------
void
PileupEventAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
PileupEventAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
PileupEventAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
PileupEventAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------

void
PileupEventAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
PileupEventAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupEventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupEventAnalyzer);
