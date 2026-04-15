// -*- C++ -*-
//
// Package:    VertexAnalyzer
// Class:      VertexAnalyzer
//
/**\class VertexAnalyzer VertexAnalyzer.cc
 Vertex/VertexAnalyzer/src/VertexAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Sat Jun 12 11:03:58 CEST 2021
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TFile.h"
#include "TTree.h"
//
// class declaration
//

class GenVertexAnalyzer : public edm::EDAnalyzer {
  public:
    explicit GenVertexAnalyzer(const edm::ParameterSet &);
    ~GenVertexAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event &, const edm::EventSetup &);
    virtual void endJob();

    virtual void beginRun(edm::Run const &, edm::EventSetup const &);
    virtual void endRun(edm::Run const &, edm::EventSetup const &);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const &,
                                      edm::EventSetup const &);
    virtual void endLuminosityBlock(edm::LuminosityBlock const &,
                                    edm::EventSetup const &);

    // ----------member data ---------------------------

    TTree *mtree;
    double Sim_PV_x;
    double Sim_PV_y;
    double Sim_PV_z;
    // Reconstruction
    double Rec_PV_x;
    double Rec_PV_y;
    double Rec_PV_z;
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
GenVertexAnalyzer::GenVertexAnalyzer(const edm::ParameterSet &iConfig)

{
    // now do what ever initialization is needed

    edm::Service<TFileService> fs;
    mtree = fs->make<TTree>("Events", "Events");
    mtree->Branch("Sim_PV_x", &Sim_PV_x);
    mtree->GetBranch("Sim_PV_x")
        ->SetTitle("Simulation primary vertex x coordinate");
    mtree->Branch("Sim_PV_y", &Sim_PV_y);
    mtree->GetBranch("Sim_PV_y")
        ->SetTitle("Simulation primary vertex y coordinate");
    mtree->Branch("Sim_PV_z", &Sim_PV_z);
    mtree->GetBranch("Sim_PV_z")
        ->SetTitle("Simulation primary vertex z coordinate");
    // Reconstruction
    mtree->Branch("Rec_PV_x", &Rec_PV_x);
    mtree->GetBranch("Rec_PV_x")
        ->SetTitle("Reconstructed primary vertex x coordinate");
    mtree->Branch("Rec_PV_y", &Rec_PV_y);
    mtree->GetBranch("Rec_PV_y")
        ->SetTitle("Reconstructed primary vertex y coordinate");
    mtree->Branch("Rec_PV_z", &Rec_PV_z);
    mtree->GetBranch("Rec_PV_z")
        ->SetTitle("Reconstructed primary vertex z coordinate");
}

GenVertexAnalyzer::~GenVertexAnalyzer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void GenVertexAnalyzer::analyze(const edm::Event &iEvent,
                                const edm::EventSetup &iSetup) {
    using namespace edm;

    Sim_PV_x = -1;
    Sim_PV_y = -1;
    Sim_PV_z = -1;

    Handle<reco::GenParticleCollection> gens;
    iEvent.getByLabel("genParticles", gens);

    Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);

    math::XYZPoint pv(vertices->begin()->position());
    const reco::Vertex &PV = vertices->front();
    Rec_PV_x = PV.x();
    Rec_PV_y = PV.y();
    Rec_PV_z = PV.z();
    if (gens.isValid() && !gens->empty()) {
        for (size_t i = 0; i < gens->size(); ++i) {
            const reco::GenParticle &candidate = (*gens)[i];

            if (candidate.status() == 1) {
                Sim_PV_x = candidate.vx();
                Sim_PV_y = candidate.vy();
                Sim_PV_z = candidate.vz();
                break; // Found it, stop looping
            }
        }
    }
    mtree->Fill();
}

// ------------ method called once each job just before starting event loop
// ------------
void GenVertexAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void GenVertexAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------
void GenVertexAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void GenVertexAnalyzer::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block
// ------------
void GenVertexAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &,
                                             edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block
// ------------
void GenVertexAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &,
                                           edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void GenVertexAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions) {
    // The following says we do not know what parameters are allowed so do no
    // validation
    //  Please change this to state exactly what you do use, even if it is no
    //  parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(GenVertexAnalyzer);
