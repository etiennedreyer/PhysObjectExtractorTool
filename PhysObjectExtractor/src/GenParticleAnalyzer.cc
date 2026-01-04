// -*- C++ -*-
//
// Package:    GenParticleAnalyzer
// Class:      GenParticleAnalyzer
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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// classes to save data
#include "TTree.h"
#include "TFile.h"
#include <vector>


class GenFlavorLabeler {
public:
    enum Label {
        kPileup = 0,
        kFake = 1,
        kPrimary = 2,
        kFromB = 3,
        kFromBC = 4,
        kFromC = 5,
        kFromTau = 6,
        kOtherSecondary = 7
    };

    static int getLabel(const reco::GenParticle& particle) {
        // --- 0. Pile-up & 1. Fake ---
        // GenParticles in the standard collection are real (not fakes).
        // Standard MiniAOD/AOD GenParticles are typically signal-only (Hard Process).
        // If you are running on a collection that includes Mixing (Pileup), 
        // you can check the sub-event ID.
        // For standard analysis, these will rarely trigger, but here is the check:
        
        // Use a small epsilon for vertex checks if needed, but strict ancestry is safer.
        
        // --- Ancestry Analysis ---
        return determineAncestry(&particle);
    }

private:
   // Helper to check for B-Hadron (ID 500-599, 5000-5999)
   static bool isBHadron(int pdgId) {
      int aid = std::abs(pdgId);
      return (aid / 100) % 10 == 5 || (aid / 1000) % 10 == 5;
   }

    // Helper to check for C-Hadron (ID 400-499, 4000-4999)
    static bool isCHadron(int pdgId) {
        int aid = std::abs(pdgId);
        return (aid / 100) % 10 == 4 || (aid / 1000) % 10 == 4;
   }

   static int determineAncestry(const reco::Candidate* p) {
      if (!p) return kFake;

      bool hasB = false;
      bool hasC = false;
      bool hasTau = false;

      bool closestIsB = false;
      bool closestIsC = false;

      const reco::Candidate* mom = p->mother();
      
      // Loop up the chain to find Heavy Flavor ancestors
      while (mom) {
         int pid = mom->pdgId();

         if (isBHadron(pid)) {
               hasB = true;
               if (!closestIsC) closestIsB = true; 
         }
         else if (isCHadron(pid)) {
               hasC = true;
               if (!closestIsB) closestIsC = true; 
         }
         else if (std::abs(pid) == 15) {
               hasTau = true;
         }

         if (mom->numberOfMothers() > 0) mom = mom->mother(0);
         else break;
      }

      // --- Categorization ---
      if (hasB) {
         if (closestIsC) return kFromBC; // Label 4
         return kFromB;                  // Label 3
      }
      if (hasC) return kFromC;            // Label 5
      if (hasTau) return kFromTau;        // Label 6

      // --- CMSSW_5_3_X COMPATIBILITY FIX ---
      // We cannot use isPromptFinalState(). We use status codes.
      // Status 3 = Hard Process (Matrix Element)
      // Status 1 = Stable Final State
      // Status 2 = Decayed
      
      int st = p->status();

      if (st == 1 || st == 3) {
         // Distinguish "Primary" (prompt) from "OtherSecondary" (light decays like K->mu, pi->mu)
         // If the direct mother is a Light Hadron (Meson/Baryon) and NOT a heavy flavor, 
         // then this is a secondary decay.
         
         const reco::Candidate* m = p->mother();
         if (m) {
               int mpid = std::abs(m->pdgId());
               
               // Check if mother is a hadron (ID > 100) but not B, C, or Proton (2212)
               // We exclude Protons because they can be beam remnants (Primary/UE)
               bool isLightHadron = (mpid > 100 && !isBHadron(mpid) && !isCHadron(mpid) && mpid != 2212);
               
               // If it comes from a light hadron decay (e.g. Pion, Kaon), it is Secondary
               if (isLightHadron) return kOtherSecondary; // Label 7
         }
         
         // Otherwise, it implies it came from String fragmentation, Gluons, or Quarks
         return kPrimary; // Label 2
      }

      // Default fall-through for status 2 or other oddities
      return kOtherSecondary; // Label 7
   }
};

//
// class declaration
//

class GenParticleAnalyzer : public edm::EDAnalyzer
{
public:
   explicit GenParticleAnalyzer(const edm::ParameterSet &);
   ~GenParticleAnalyzer();

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
   virtual void beginJob();
   virtual void analyze(const edm::Event &, const edm::EventSetup &);
   virtual void endJob();
   virtual void beginRun(edm::Run const &, edm::EventSetup const &);
   virtual void endRun(edm::Run const &, edm::EventSetup const &);
   virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
   virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

   std::vector<std::string> particle;
   edm::InputTag jetInput;

   // ----------member data ---------------------------

   TTree *mtree;

   int numGenPart;
   std::vector<int> GenPart_status;
   std::vector<int> GenPart_label;
   std::vector<float> GenPart_pt;
   std::vector<float> GenPart_eta;
   std::vector<float> GenPart_mass;
   std::vector<int> GenPart_pdgId;
   std::vector<float> GenPart_phi;
   std::vector<float> GenPart_px;
   std::vector<float> GenPart_py;
   std::vector<float> GenPart_pz;
   std::vector<float> GenPart_vx;
   std::vector<float> GenPart_vy;
   std::vector<float> GenPart_vz;

   std::vector<int> GenPart_jetIdx;
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

GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet &iConfig) : particle(iConfig.getParameter<std::vector<std::string>>("input_particle"))

{
   // now do what ever initialization is needed
   jetInput = iConfig.getParameter<edm::InputTag>("InputJetCollection");

   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

   mtree->Branch("numGenPart", &numGenPart);
   mtree->GetBranch("numGenPart")->SetTitle("number of generator particles");
   mtree->Branch("GenPart_pt", &GenPart_pt);
   mtree->GetBranch("GenPart_pt")->SetTitle("generator particle transverse momentum");
   mtree->Branch("GenPart_eta", &GenPart_eta);
   mtree->GetBranch("GenPart_eta")->SetTitle("generator particle pseudorapidity");
   mtree->Branch("GenPart_mass", &GenPart_mass);
   mtree->GetBranch("GenPart_mass")->SetTitle("generator particle mass");
   mtree->Branch("GenPart_pdgId", &GenPart_pdgId);
   mtree->GetBranch("GenPart_pdgId")->SetTitle("generator particle PDG id");
   mtree->Branch("GenPart_phi", &GenPart_phi);
   mtree->GetBranch("GenPart_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
   // mtree->Branch("GenPart_px", &GenPart_px);
   // mtree->GetBranch("GenPart_px")->SetTitle("generator particle x coordinate of momentum vector");
   // mtree->Branch("GenPart_py", &GenPart_py);
   // mtree->GetBranch("GenPart_py")->SetTitle("generator particle y coordinate of momentum vector");
   // mtree->Branch("GenPart_pz", &GenPart_pz);
   // mtree->GetBranch("GenPart_pz")->SetTitle("generator particle z coordinate of momentum vector");
   mtree->Branch("GenPart_status", &GenPart_status);
   mtree->GetBranch("GenPart_status")->SetTitle("Particle status. 1=stable");
   
   mtree->Branch("GenPart_label", &GenPart_label);
   mtree->GetBranch("GenPart_label")->SetTitle("Particle label.");

   mtree->Branch("GenPart_vx", &GenPart_vx);
   mtree->GetBranch("GenPart_vx")->SetTitle("generator particle x coordinate its vertex");
   mtree->Branch("GenPart_vy", &GenPart_vy);
   mtree->GetBranch("GenPart_vy")->SetTitle("generator particle y coordinate its vertex");
   mtree->Branch("GenPart_vz", &GenPart_vz);
   mtree->GetBranch("GenPart_vz")->SetTitle("generator particle z coordinate its vertex");

   // mtree->Branch("GenPart_jetIdx", &GenPart_jetIdx);
   // mtree->GetBranch("GenPart_jetIdx")->SetTitle("Index of the jet the particle is in. -1 if not in a jet.");
}

GenParticleAnalyzer::~GenParticleAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void GenParticleAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;
   using namespace std;

   numGenPart = 0;
   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_phi.clear();
   GenPart_px.clear();
   GenPart_py.clear();
   GenPart_pz.clear();
   GenPart_status.clear();
   GenPart_label.clear();
   GenPart_vx.clear();
   GenPart_vy.clear();
   GenPart_vz.clear();

   GenPart_jetIdx.clear();

   Handle<reco::GenParticleCollection> gens;
   iEvent.getByLabel("genParticles", gens);

   Handle<reco::GenJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);

   unsigned int i, j;
   string s1, s2;
   std::vector<int> status_parsed;
   std::vector<int> pdgId_parsed;
   std::string delimiter = ":";

   for (i = 0; i < particle.size(); i++)
   {
      // get status and pgdId from configuration
      s1 = particle[i].substr(0, particle[i].find(delimiter));
      s2 = particle[i].substr(particle[i].find(delimiter) + 1, particle[i].size());
      // parse string to int
      status_parsed.push_back(stoi(s1));
      pdgId_parsed.push_back(stoi(s2));
   }

   if (gens.isValid())
   {
      numGenPart = gens->size();
      std::vector<const reco::GenParticle *> genJetParticles;
      std::vector<int> genJetIndices;
      if (myjets.isValid())
      {
         for (i = 0; i < myjets->size(); i++)
         {
            reco::GenJet jet = myjets->at(i);
            if (jet.pt() < 20)
            {
               continue;
            }
            std::vector<const reco::GenParticle *> genJetParticles_temp = jet.getGenConstituents();
            for (j = 0; j < genJetParticles_temp.size(); j++)
            {
               genJetParticles.push_back(genJetParticles_temp[j]);
               genJetIndices.push_back(i);
            }
         }
      }
      for (reco::GenParticleCollection::const_iterator itGenPart = gens->begin(); itGenPart != gens->end(); ++itGenPart)
      {
         // loop trough all particles selected in configuration
         for (i = 0; i < particle.size(); i++)
         {
            if (
                (status_parsed[i] == itGenPart->status() && pdgId_parsed[i] == itGenPart->pdgId()) ||
                (status_parsed[i] == 0 && pdgId_parsed[i] == 0) ||
                (status_parsed[i] == 0 && pdgId_parsed[i] == itGenPart->pdgId()) ||
                (status_parsed[i] == itGenPart->status() && pdgId_parsed[i] == 0))
            {
               GenPart_pt.push_back(itGenPart->pt());
               GenPart_eta.push_back(itGenPart->eta());
               GenPart_mass.push_back(itGenPart->mass());
               GenPart_pdgId.push_back(itGenPart->pdgId());
               GenPart_phi.push_back(itGenPart->phi());
               GenPart_status.push_back(itGenPart->status());
               int label = GenFlavorLabeler::getLabel(*itGenPart);
               GenPart_label.push_back(label);
               GenPart_px.push_back(itGenPart->px());
               GenPart_py.push_back(itGenPart->py());
               GenPart_pz.push_back(itGenPart->pz());

               GenPart_vx.push_back(itGenPart->vx());
               GenPart_vy.push_back(itGenPart->vy());
               GenPart_vz.push_back(itGenPart->vz());
               
               bool flag = false;
               for (j = 0; j < genJetParticles.size(); j++)
               {
                  const reco::GenParticle *genPart = &(*itGenPart);
                  if (genPart == genJetParticles[j])
                  {
                     GenPart_jetIdx.push_back(genJetIndices[j]);
                     // std::cout << "Particle " << itGenPart->pdgId() << " is in jet " << genJetIndices[j] << std::endl;
                     flag = true;
                     break;
                  }
               }
               if (!flag)
               {
                  GenPart_jetIdx.push_back(-1);
               }
            }
         }
      }
   }

   mtree->Fill();
   return;
}

// ------------ method called once each job just before starting event loop  ------------
void GenParticleAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenParticleAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void GenParticleAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a run  ------------
void GenParticleAnalyzer::endRun(edm::Run const &, edm::EventSetup const &)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void GenParticleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenParticleAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   // The following says we do not know what parameters are allowed so do no validation
   //  Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);