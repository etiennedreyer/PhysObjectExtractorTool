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
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// GSFTrack 
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
// Muon
#include "DataFormats/MuonReco/interface/Muon.h"
// classes to save data
#include "TTree.h"
#include "TFile.h"
#include <vector>

//TransientTrack and IPTools for impact parameter
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// GenParticle for flavor labeling
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <cmath>

class PFFlavorLabeler {
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

    static int getLabel(const reco::PFCandidate& pfCand, const reco::GenParticleCollection* genParticles) {
        if (!genParticles) return kFake;
        
        // Find the closest GenParticle by deltaR
        const reco::GenParticle* bestMatch = nullptr;
        double minDeltaR = 0.1; // matching cone
        
        for (reco::GenParticleCollection::const_iterator itGen = genParticles->begin(); 
             itGen != genParticles->end(); ++itGen) {
            
            double deta = pfCand.eta() - itGen->eta();
            double dphi = deltaPhi(pfCand.phi(), itGen->phi());
            double dr = std::sqrt(deta*deta + dphi*dphi);
            
            if (dr < minDeltaR) {
                minDeltaR = dr;
                bestMatch = &(*itGen);
            }
        }
        
        // If no match found, label as fake
        if (!bestMatch) return kFake;
        
        // Apply GenParticle-based labeling logic
        return determineAncestry(bestMatch);
    }

private:
    static double deltaPhi(double phi1, double phi2) {
        double result = phi1 - phi2;
        while (result > M_PI) result -= 2*M_PI;
        while (result <= -M_PI) result += 2*M_PI;
        return result;
    }

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

        // Status-based categorization
        int st = p->status();

        if (st == 1 || st == 3) {
            const reco::Candidate* m = p->mother();
            if (m) {
                int mpid = std::abs(m->pdgId());
                bool isLightHadron = (mpid > 100 && !isBHadron(mpid) && !isCHadron(mpid) && mpid != 2212);
                if (isLightHadron) return kOtherSecondary; // Label 7
            }
            return kPrimary; // Label 2
        }

        return kOtherSecondary; // Label 7
    }
};

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
   std::vector<int> PFCand_label;
   
   std::vector<float> PFCand_d0;
   std::vector<float> PFCand_d0Error;
   std::vector<float> PFCand_z0;
   std::vector<float> PFCand_z0Error;

   std::vector<float> PFCand_ip3d;
   std::vector<float> PFCand_ip3dError;
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
   // mtree->Branch("PFCand_px", &PFCand_px);
   // mtree->GetBranch("PFCand_px")->SetTitle("pflow candidate x coordinate of momentum vector");
   // mtree->Branch("PFCand_py", &PFCand_py);
   // mtree->GetBranch("PFCand_py")->SetTitle("pflow candidate y coordinate of momentum vector");
   // mtree->Branch("PFCand_pz", &PFCand_pz);
   // mtree->GetBranch("PFCand_pz")->SetTitle("pflow candidate z coordinate of momentum vector");

   mtree->Branch("PFCand_vx", &PFCand_vx);
   mtree->GetBranch("PFCand_vx")->SetTitle("pflow candidate x coordinate of its vertex");
   mtree->Branch("PFCand_vy", &PFCand_vy);
   mtree->GetBranch("PFCand_vy")->SetTitle("pflow candidate y coordinate of its vertex");
   mtree->Branch("PFCand_vz", &PFCand_vz);
   mtree->GetBranch("PFCand_vz")->SetTitle("pflow candidate z coordinate of its vertex");
   // mtree->Branch("PFCand_jetIdx", &PFCand_jetIdx);
   // mtree->GetBranch("PFCand_jetIdx")->SetTitle("Index of the jet the particle is in. -1 if not in a jet.");
   
   
   // Impact parameter variables
   mtree->Branch("PFCand_d0", &PFCand_d0);
   mtree->GetBranch("PFCand_d0")->SetTitle("pflow candidate d0");
   mtree->Branch("PFCand_d0Error", &PFCand_d0Error);
   mtree->GetBranch("PFCand_d0Error")->SetTitle("pflow candidate d0Error");

   mtree->Branch("PFCand_z0", &PFCand_z0);
   mtree->GetBranch("PFCand_z0")->SetTitle("pflow candidate z0");
   mtree->Branch("PFCand_z0Error", &PFCand_z0Error);
   mtree->GetBranch("PFCand_z0Error")->SetTitle("pflow candidate z0Error");
   
   mtree->Branch("PFCand_ip3d",&PFCand_ip3d);
   mtree->GetBranch("PFCand_ip3d")->SetTitle("PFCand ip3d");
   mtree->Branch("PFCand_ip3dError",&PFCand_ip3dError);
   mtree->GetBranch("PFCand_ip3dError")->SetTitle("PFCand ip3dError");
   
   mtree->Branch("PFCand_label", &PFCand_label);
   mtree->GetBranch("PFCand_label")->SetTitle("PFCand flavor label");
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

   PFCand_d0.clear();
   PFCand_d0Error.clear();
   PFCand_z0.clear();
   PFCand_z0Error.clear();

   PFCand_ip3d.clear();
   PFCand_ip3dError.clear();

   PFCand_jetIdx.clear();
   PFCand_label.clear();

   Handle<reco::PFCandidateCollection> pfcs;
   iEvent.getByLabel("particleFlow", pfcs);

   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);

   Handle<reco::PFJetCollection> myjets;
   iEvent.getByLabel(jetInput, myjets);

   // OverlapChecker overlap = OverlapChecker();
   Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
   math::XYZPoint pv(vertices->begin()->position());
   const reco::Vertex &PV = vertices->front();
   // unsigned int i, j;
   unsigned int i;

   if (pfcs.isValid())
   {
      // numPFCand = pfcs->size();
      numPFCand = 0;

      // std::vector<reco::PFCandidatePtr> pfJetParticles;
      // std::vector<int> pfJetIndices;
      // if (myjets.isValid())
      // {
      //    for (i = 0; i < myjets->size(); i++)
      //    {
      //       reco::PFJet jet = myjets->at(i);
      //       if (jet.pt() < 20)
      //       {
      //          continue;
      //       }
      //       std::vector<reco::PFCandidatePtr> pfJetParticles_temp = jet.getPFConstituents();
      //       for (j = 0; j < pfJetParticles_temp.size(); j++)
      //       {
      //          pfJetParticles.push_back(pfJetParticles_temp[j]);
      //          pfJetIndices.push_back(i);
      //       }
      //    }
      // }
      for (reco::PFCandidateCollection::const_iterator itPFCand = pfcs->begin(); itPFCand != pfcs->end(); ++itPFCand)
      {
         // loop trough all particles selected in configuration
         for (i = 0; i < particle.size(); i++)
         {
            if ((itPFCand->pdgId() == 0) || (itPFCand->pt() < 1))
            {
               continue;
            }
            numPFCand++;
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

            // bool flag = false;
            // for (j = 0; j < pfJetParticles.size(); j++)
            // {
            //    if (&(*itPFCand) == pfJetParticles[j].get())
            //    {
            //       PFCand_jetIdx.push_back(pfJetIndices[j]);
            //       flag = true;
            //       break;
            //    }
            // }
            // if (!flag)
            // {
            //    PFCand_jetIdx.push_back(-1);
            // }
            
            // Compute flavor label
            int label = PFFlavorLabeler::kFake;
            if (genParticles.isValid()) {
                label = PFFlavorLabeler::getLabel(*itPFCand, genParticles.product());
            }
            PFCand_label.push_back(label);
            
            bool isElectron = itPFCand->particleId() == reco::PFCandidate::e;
            bool isMuon = itPFCand->particleId() == reco::PFCandidate::mu;
            bool ipCalculated = false;

            if (isMuon){
               reco::MuonRef muon = itPFCand->muonRef();
               if (muon.isNonnull()){
                  // reco::TrackRef muonTrack = muon->innerTrack();
                  reco::TrackRef muonTrack = muon->muonBestTrack();
                  if (muonTrack.isNonnull()){
                     PFCand_d0.push_back(muonTrack->dxy(PV.position()));
                     PFCand_z0.push_back(muonTrack->dz(PV.position()));
                     PFCand_d0Error.push_back(muonTrack->dxyError());
                     PFCand_z0Error.push_back(muonTrack->dzError());

                     edm::ESHandle<TransientTrackBuilder> trackBuilder;
                     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
                     reco::TransientTrack tt = trackBuilder->build(muonTrack);
                     std::pair<bool,Measurement1D> ip3dpv = IPTools::absoluteImpactParameter3D(tt, PV);
                     PFCand_ip3d.push_back(ip3dpv.second.value());
                     PFCand_ip3dError.push_back(ip3dpv.second.significance());
                     ipCalculated = true;
                  }
               }
            }
            else if (isElectron){
               reco::GsfElectronRef gsfElectron = itPFCand->gsfElectronRef();
               if (gsfElectron.isNonnull()){
                  reco::GsfTrackRef gsfTrack = gsfElectron->gsfTrack();
                  if (gsfTrack.isNonnull()){
                     PFCand_d0.push_back(gsfTrack->dxy(PV.position()));
                     PFCand_z0.push_back(gsfTrack->dz(PV.position()));
                     PFCand_d0Error.push_back(gsfTrack->dxyError());
                     PFCand_z0Error.push_back(gsfTrack->dzError());

                     edm::ESHandle<TransientTrackBuilder> trackBuilder;
                     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
                     reco::TransientTrack tt = trackBuilder->build(gsfTrack);
                     std::pair<bool,Measurement1D> ip3dpv = IPTools::absoluteImpactParameter3D(tt, PV);
                     PFCand_ip3d.push_back(ip3dpv.second.value());
                     PFCand_ip3dError.push_back(ip3dpv.second.significance());
                     ipCalculated = true;
                  }
               }
            }
            if (ipCalculated){
               continue;
            }
            reco::TrackRef track = itPFCand->trackRef();
            reco::GsfTrackRef trackGsf = itPFCand->gsfTrackRef();
            if (track.isNonnull()){
               PFCand_d0.push_back(track->dxy(PV.position()));
               PFCand_z0.push_back(track->dz(PV.position()));
               PFCand_d0Error.push_back(track->dxyError());
               PFCand_z0Error.push_back(track->dzError());

               edm::ESHandle<TransientTrackBuilder> trackBuilder;
               iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
               reco::TransientTrack tt = trackBuilder->build(track);
               std::pair<bool,Measurement1D> ip3dpv = IPTools::absoluteImpactParameter3D(tt, PV);
               PFCand_ip3d.push_back(ip3dpv.second.value());
               PFCand_ip3dError.push_back(ip3dpv.second.significance());
            }
            else if (trackGsf.isNonnull()){
               PFCand_d0.push_back(trackGsf->dxy(PV.position()));
               PFCand_z0.push_back(trackGsf->dz(PV.position()));
               PFCand_d0Error.push_back(trackGsf->dxyError());
               PFCand_z0Error.push_back(trackGsf->dzError());

               edm::ESHandle<TransientTrackBuilder> trackBuilder;
               iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
               reco::TransientTrack tt = trackBuilder->build(trackGsf);
               std::pair<bool,Measurement1D> ip3dpv = IPTools::absoluteImpactParameter3D(tt, PV);
               PFCand_ip3d.push_back(ip3dpv.second.value());
               PFCand_ip3dError.push_back(ip3dpv.second.significance());
            }
            else {
               PFCand_d0.push_back(-1000);
               PFCand_z0.push_back(-1000);
               PFCand_d0Error.push_back(-1000);
               PFCand_z0Error.push_back(-1000);
               PFCand_ip3d.push_back(-1000);
               PFCand_ip3dError.push_back(-1000);
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
