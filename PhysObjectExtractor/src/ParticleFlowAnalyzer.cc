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

//classes to extract ParticleFlowCandidate information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>

//
// class declaration
//

class ParticleFlowAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ParticleFlowAnalyzer(const edm::ParameterSet&);
      ~ParticleFlowAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

     std::vector<std::string>  particle;

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

ParticleFlowAnalyzer::ParticleFlowAnalyzer(const edm::ParameterSet& iConfig):
particle(iConfig.getParameter<std::vector<std::string> >("input_particle"))

{
//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	mtree = fs->make<TTree>("Events", "Events");
    
    mtree->Branch("numPFCand",&numPFCand);
    mtree->GetBranch("numPFCand")->SetTitle("number of pfc particles");
    mtree->Branch("PFCand_pt",&PFCand_pt);
    mtree->GetBranch("PFCand_pt")->SetTitle("pflow candidate transverse momentum");
    mtree->Branch("PFCand_eta",&PFCand_eta);
    mtree->GetBranch("PFCand_eta")->SetTitle("pflow candidate pseudorapidity");
    mtree->Branch("PFCand_mass",&PFCand_mass);
    mtree->GetBranch("PFCand_mass")->SetTitle("pflow candidate mass");
    mtree->Branch("PFCand_pdgId",&PFCand_pdgId);
    mtree->GetBranch("PFCand_pdgId")->SetTitle("pflow candidate PDG id");
    mtree->Branch("PFCand_phi",&PFCand_phi);
    mtree->GetBranch("PFCand_phi")->SetTitle("pflow candidate azimuthal angle of momentum vector");
    mtree->Branch("PFCand_px",&PFCand_px);
    mtree->GetBranch("PFCand_px")->SetTitle("pflow candidate x coordinate of momentum vector");
    mtree->Branch("PFCand_py",&PFCand_py);
    mtree->GetBranch("PFCand_py")->SetTitle("pflow candidate y coordinate of momentum vector");
    mtree->Branch("PFCand_pz",&PFCand_pz);
    mtree->GetBranch("PFCand_pz")->SetTitle("pflow candidate z coordinate of momentum vector");

    mtree->Branch("PFCand_vx",&PFCand_vx);
    mtree->GetBranch("PFCand_vx")->SetTitle("pflow candidate x coordinate of its vertex");
    mtree->Branch("PFCand_vy",&PFCand_vy);
    mtree->GetBranch("PFCand_vy")->SetTitle("pflow candidate y coordinate of its vertex");
    mtree->Branch("PFCand_vz",&PFCand_vz);
    mtree->GetBranch("PFCand_vz")->SetTitle("pflow candidate z coordinate of its vertex");

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
void
ParticleFlowAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   
   numPFCand=0;
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

   Handle<reco::PFCandidateCollection> pfcs;
   iEvent.getByLabel("particleFlow", pfcs);
   
   unsigned int i;
   // string s1,s2;
   // std::vector<int> status_parsed;
   // std::vector<int> pdgId_parsed;
   // std::string delimiter = ":";

   // for(i=0;i<particle.size();i++)
   // {
   //     //get status and pgdId from configuration
   //     s1=particle[i].substr(0,particle[i].find(delimiter));
   //     s2=particle[i].substr(particle[i].find(delimiter)+1,particle[i].size());
   //     //parse string to int
   //     status_parsed.push_back(stoi(s1));
   //     pdgId_parsed.push_back(stoi(s2));
   // }
  
  if(pfcs.isValid())
  {
      numPFCand=pfcs->size();
      for (reco::PFCandidateCollection::const_iterator itPFCand=pfcs->begin(); itPFCand!=pfcs->end(); ++itPFCand)
      {
         //loop trough all particles selected in configuration
         for(i = 0; i < particle.size(); i++)
         {
         //    if((status_parsed[i]==itPFCand->status() && pdgId_parsed[i]==itPFCand->pdgId())||(status_parsed[i]==0 && pdgId_parsed[i]==0))
         //  {
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
         //  }
         }               
      }
  }
	
  mtree->Fill();
  return;

}

// ------------ method called once each job just before starting event loop  ------------
void
ParticleFlowAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
ParticleFlowAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
ParticleFlowAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
ParticleFlowAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
ParticleFlowAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
ParticleFlowAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ParticleFlowAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ParticleFlowAnalyzer);