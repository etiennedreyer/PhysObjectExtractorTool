// -*- C++ -*-
//
// Package:    GenHadronAnalyzer
// Class:      GenHadronAnalyzer
//

// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// classes to extract GenParticle information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// classes to save data
#include "TTree.h"
#include "TFile.h"
#include <vector>

//
// class declaration
//

class GenHadronAnalyzer : public edm::EDAnalyzer
{
public:
   explicit GenHadronAnalyzer(const edm::ParameterSet &);
   ~GenHadronAnalyzer();

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
   virtual void beginJob();
   virtual void analyze(const edm::Event &, const edm::EventSetup &);
   virtual void endJob();
   virtual void beginRun(edm::Run const &, edm::EventSetup const &);
   virtual void endRun(edm::Run const &, edm::EventSetup const &);
   virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
   virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

   // ----------member data ---------------------------

   TTree *mtree;

   int numGenPart;
   std::vector<int> GenPart_status;
   std::vector<float> GenPart_pt;
   std::vector<float> GenPart_eta;
   std::vector<float> GenPart_mass;
   std::vector<int> GenPart_pdgId;
   std::vector<float> GenPart_phi;
   std::vector<float> GenPart_vx;
   std::vector<float> GenPart_vy;
   std::vector<float> GenPart_vz;
   std::vector<int> GenPart_label;
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

GenHadronAnalyzer::GenHadronAnalyzer(const edm::ParameterSet &iConfig) 
{
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

   mtree->Branch("numHeavyHadron", &numGenPart);
   mtree->GetBranch("numHeavyHadron")->SetTitle("number of generator particles");
   mtree->Branch("HeavyHadron_pt", &GenPart_pt);
   mtree->GetBranch("HeavyHadron_pt")->SetTitle("generator particle transverse momentum");
   mtree->Branch("HeavyHadron_eta", &GenPart_eta);
   mtree->GetBranch("HeavyHadron_eta")->SetTitle("generator particle pseudorapidity");
   mtree->Branch("HeavyHadron_mass", &GenPart_mass);
   mtree->GetBranch("HeavyHadron_mass")->SetTitle("generator particle mass");
   mtree->Branch("HeavyHadron_pdgId", &GenPart_pdgId);
   mtree->GetBranch("HeavyHadron_pdgId")->SetTitle("generator particle PDG id");
   mtree->Branch("HeavyHadron_phi", &GenPart_phi);
   mtree->GetBranch("HeavyHadron_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
   mtree->Branch("HeavyHadron_status", &GenPart_status);
   mtree->GetBranch("HeavyHadron_status")->SetTitle("Particle status. 1=stable");
   
   mtree->Branch("HeavyHadron_label", &GenPart_label);
   mtree->GetBranch("HeavyHadron_label")->SetTitle("Particle label. 4 for C-hadron, 5 for B-hadron");

   mtree->Branch("HeavyHadron_vx", &GenPart_vx);
   mtree->GetBranch("HeavyHadron_vx")->SetTitle("generator particle x coordinate its vertex");
   mtree->Branch("HeavyHadron_vy", &GenPart_vy);
   mtree->GetBranch("HeavyHadron_vy")->SetTitle("generator particle y coordinate its vertex");
   mtree->Branch("HeavyHadron_vz", &GenPart_vz);
   mtree->GetBranch("HeavyHadron_vz")->SetTitle("generator particle z coordinate its vertex");
}

GenHadronAnalyzer::~GenHadronAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void GenHadronAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;
   using namespace std;

   numGenPart = 0;
   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_phi.clear();

   GenPart_status.clear();
   GenPart_label.clear();
   GenPart_vx.clear();
   GenPart_vy.clear();
   GenPart_vz.clear();


   Handle<reco::GenParticleCollection> gens;
   iEvent.getByLabel("genParticles", gens);

   if (gens.isValid())
   {
      for (reco::GenParticleCollection::const_iterator itGenPart = gens->begin(); itGenPart != gens->end(); ++itGenPart)
      {
         const reco::GenParticle & p = *itGenPart;

         int id = std::abs(p.pdgId());
         int status = p.status();

         // Quick filters first for performance
         // 1. Must be Status 2 (Decayed) for Pythia 6 heavy hadrons
         if (status != 2) continue;

         // 2. Filter only B or C Hadrons (check PDG ID structure)
         bool hasB = ((id/100)%10 == 5) || ((id/1000)%10 == 5);
         bool hasC = ((id/100)%10 == 4) || ((id/1000)%10 == 4);

         if (!hasB && !hasC) continue;

         // 3. "Last Copy" Check
         // We look at the daughters to see if this hadron "decays" into itself
         bool isLastCopy = true;
         size_t nDaughters = p.numberOfDaughters();
         for (size_t d = 0; d < nDaughters; ++d) {
            const reco::Candidate* daughter = p.daughter(d);
            int daughterId = std::abs(daughter->pdgId());
            
            // Exact match of ID implies this is just an intermediate step
            // (e.g. B0 -> B0 + gamma)
            if (daughterId == id) {
               isLastCopy = false;
               break;
            }
         }
         
         if (isLastCopy) {
            int label = hasB ? 5 : 4; // 5 for B-hadron, 4 for C-hadron
            GenPart_pt.push_back(p.pt());
            GenPart_eta.push_back(p.eta());
            GenPart_mass.push_back(p.mass());
            GenPart_pdgId.push_back(p.pdgId());
            GenPart_phi.push_back(p.phi());
            GenPart_status.push_back(p.status());
            GenPart_label.push_back(label);

            GenPart_vx.push_back(p.vx());
            GenPart_vy.push_back(p.vy());
            GenPart_vz.push_back(p.vz());
            ++numGenPart;
         }
         
      }
   }

   mtree->Fill();
   return;
}

// ------------ method called once each job just before starting event loop  ------------
void GenHadronAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenHadronAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void GenHadronAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a run  ------------
void GenHadronAnalyzer::endRun(edm::Run const &, edm::EventSetup const &)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void GenHadronAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenHadronAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenHadronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   // The following says we do not know what parameters are allowed so do no validation
   //  Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(GenHadronAnalyzer);