// -*- C++ -*-
//
// Package:    ParticleFlowAnalyzer
// Class:      ParticleFlowAnalyzer
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

// classes to extract ParticleFlowCandidate information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
// #include "DataFormats/Candidate/interface/OverlapChecker.h"

// classes to save data
#include "TTree.h"
#include "TFile.h"
#include <vector>

//
// class declaration
//

class ParticleFlowAnalyzer : public edm::EDAnalyzer
{
public:
   explicit ParticleFlowAnalyzer(const edm::ParameterSet &);
   ~ParticleFlowAnalyzer();

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

   int numPFCand;
   std::vector<float> PFCand_pt;
   std::vector<float> PFCand_eta;
   std::vector<float> PFCand_mass;
   std::vector<int> PFCand_pdgId;
   std::vector<float> PFCand_phi;
   std::vector<float> PFCand_px;
   std::vector<float> PFCand_py;
   std::vector<float> PFCand_pz;

   std::vector<float> PFCand_vx;
   std::vector<float> PFCand_vy;
   std::vector<float> PFCand_vz;

   std::vector<int> PFCand_jetIdx;
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

ParticleFlowAnalyzer::ParticleFlowAnalyzer(const edm::ParameterSet &iConfig) : particle(iConfig.getParameter<std::vector<std::string>>("input_particle"))

{
   // now do what ever initialization is needed
   jetInput = iConfig.getParameter<edm::InputTag>("InputJetCollection");
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");

   mtree->Branch("numPFCand", &numPFCand);
   mtree->GetBranch("numPFCand")->SetTitle("number of pfc particles");
   mtree->Branch("PFCand_pt", &PFCand_pt);
   mtree->GetBranch("PFCand_pt")->SetTitle("pflow candidate transverse momentum");
   mtree->Branch("PFCand_eta", &PFCand_eta);
   mtree->GetBranch("PFCand_eta")->SetTitle("pflow candidate pseudorapidity");
   mtree->Branch("PFCand_mass", &PFCand_mass);
   mtree->GetBranch("PFCand_mass")->SetTitle("pflow candidate mass");
   mtree->Branch("PFCand_pdgId", &PFCand_pdgId);
   mtree->GetBranch("PFCand_pdgId")->SetTitle("pflow candidate PDG id");
   mtree->Branch("PFCand_phi", &PFCand_phi);
   mtree->GetBranch("PFCand_phi")->SetTitle("pflow candidate azimuthal angle of momentum vector");
   mtree->Branch("PFCand_px", &PFCand_px);
   mtree->GetBranch("PFCand_px")->SetTitle("pflow candidate x coordinate of momentum vector");
   mtree->Branch("PFCand_py", &PFCand_py);
   mtree->GetBranch("PFCand_py")->SetTitle("pflow candidate y coordinate of momentum vector");
   mtree->Branch("PFCand_pz", &PFCand_pz);
   mtree->GetBranch("PFCand_pz")->SetTitle("pflow candidate z coordinate of momentum vector");

   mtree->Branch("PFCand_vx", &PFCand_vx);
   mtree->GetBranch("PFCand_vx")->SetTitle("pflow candidate x coordinate of its vertex");
   mtree->Branch("PFCand_vy", &PFCand_vy);
   mtree->GetBranch("PFCand_vy")->SetTitle("pflow candidate y coordinate of its vertex");
   mtree->Branch("PFCand_vz", &PFCand_vz);
   mtree->GetBranch("PFCand_vz")->SetTitle("pflow candidate z coordinate of its vertex");
   mtree->Branch("PFCand_jetIdx", &PFCand_jetIdx);
   mtree->GetBranch("PFCand_jetIdx")->SetTitle("Index of the jet the particle is in. -1 if not in a jet.");
}

ParticleFlowAnalyzer::~ParticleFlowAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void ParticleFlowAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;
   using namespace std;

   numPFCand = 0;
   PFCand_pt.clear();
   PFCand_eta.clear();
   PFCand_mass.clear();
   PFCand_pdgId.clear();
   PFCand_phi.clear();
   PFCand_px.clear();
   PFCand_py.clear();
   PFCand_pz.clear();

   PFCand_vx.clear();
   PFCand_vy.clear();
   PFCand_vz.clear();

   PFCand_jetIdx.clear();

   Handle<reco::PFCandidateCollection> pfcs;
   iEvent.getByLabel("particleFlow", pfcs);

   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);

   // OverlapChecker overlap = OverlapChecker();

   unsigned int i, j;

   if (pfcs.isValid())
   {
      numPFCand = pfcs->size();

      std::vector<reco::PFCandidatePtr> pfJetParticles;
      std::vector<int> pfJetIndices;
      if (myjets.isValid())
      {
         for (i = 0; i < myjets->size(); i++)
         {
            reco::PFJet jet = myjets->at(i);
            if (jet.pt() < 20)
            {
               continue;
            }
            std::vector<reco::PFCandidatePtr> pfJetParticles_temp = jet.getPFConstituents();
            for (j = 0; j < pfJetParticles_temp.size(); j++)
            {
               pfJetParticles.push_back(pfJetParticles_temp[j]);
               pfJetIndices.push_back(i);
            }
         }
      }
      for (reco::PFCandidateCollection::const_iterator itPFCand = pfcs->begin(); itPFCand != pfcs->end(); ++itPFCand)
      {
         // loop trough all particles selected in configuration
         for (i = 0; i < particle.size(); i++)
         {
            PFCand_pt.push_back(itPFCand->pt());
            PFCand_eta.push_back(itPFCand->eta());
            PFCand_mass.push_back(itPFCand->mass());
            PFCand_pdgId.push_back(itPFCand->pdgId());
            PFCand_phi.push_back(itPFCand->phi());
            PFCand_px.push_back(itPFCand->px());
            PFCand_py.push_back(itPFCand->py());
            PFCand_pz.push_back(itPFCand->pz());

            PFCand_vx.push_back(itPFCand->vx());
            PFCand_vy.push_back(itPFCand->vy());
            PFCand_vz.push_back(itPFCand->vz());

            bool flag = false;
            for (j = 0; j < pfJetParticles.size(); j++)
            {
               if (&(*itPFCand) == pfJetParticles[j].get())
               {
                  PFCand_jetIdx.push_back(pfJetIndices[j]);
                  flag = true;
                  break;
               }
            }
            if (!flag)
            {
               PFCand_jetIdx.push_back(-1);
            }
         }
      }
   }

   mtree->Fill();
   return;
}

// ------------ method called once each job just before starting event loop  ------------
void ParticleFlowAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void ParticleFlowAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void ParticleFlowAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a run  ------------
void ParticleFlowAnalyzer::endRun(edm::Run const &, edm::EventSetup const &)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
void ParticleFlowAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void ParticleFlowAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ParticleFlowAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   // The following says we do not know what parameters are allowed so do no validation
   //  Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(ParticleFlowAnalyzer);