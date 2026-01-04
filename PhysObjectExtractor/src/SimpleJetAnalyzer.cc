// -*- C++ -*-
//
// Package:    JetAnalyzer
// Class:      JetAnalyzer
//

// system include files
#include <TMath.h>
#include <memory>
// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

// classes to extract PFJet information
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

// classes to save data
#include "TFile.h"
#include "TTree.h"
#include <vector>

#include "TRandom3.h"

//
// class declaration
//

class SimpleJetAnalyzer : public edm::EDAnalyzer {
public:
  explicit SimpleJetAnalyzer(const edm::ParameterSet &);
  ~SimpleJetAnalyzer();

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
                                  edm::EventSetup const &) override;

  // declare the input tag for PFJetCollection
  edm::InputTag jetInput;

  // ----------member data ---------------------------
  // jec variables
  bool isData;
  float min_pt;

  int numjet; // number of jets in the event
  TTree *mtree;
  std::vector<float> jet_e;
  std::vector<float> jet_pt;
  std::vector<float> jet_px;
  std::vector<float> jet_py;
  std::vector<float> jet_pz;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_ch;
  std::vector<float> jet_mass;
  std::vector<double> jet_btag;
  std::vector<int> jet_flavour;
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

SimpleJetAnalyzer::SimpleJetAnalyzer(const edm::ParameterSet &iConfig) {
  // now do what ever initialization is needed
  jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");

  isData = iConfig.getParameter<bool>("isData");
  min_pt = iConfig.getParameter<double>("minPt");

  mtree->Branch("numberjet", &numjet);
  mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
  mtree->Branch("jet_e", &jet_e);
  mtree->GetBranch("jet_e")->SetTitle("Uncorrected Jet Energy");
  mtree->Branch("jet_pt", &jet_pt);
  mtree->GetBranch("jet_pt")->SetTitle("Uncorrected Transverse Jet Momentum");
  mtree->Branch("jet_px", &jet_px);
  mtree->GetBranch("jet_px")->SetTitle("X-Component of Jet Momentum");
  mtree->Branch("jet_py", &jet_py);
  mtree->GetBranch("jet_py")->SetTitle("Y-Component of Jet Momentum");
  mtree->Branch("jet_pz", &jet_pz);
  mtree->GetBranch("jet_pz")->SetTitle("Z-Component of Jet Momentum");
  mtree->Branch("jet_eta", &jet_eta);
  mtree->GetBranch("jet_eta")->SetTitle("Jet Eta");
  mtree->Branch("jet_phi", &jet_phi);
  mtree->GetBranch("jet_phi")->SetTitle("Jet Phi");
  mtree->Branch("jet_ch", &jet_ch);
  mtree->GetBranch("jet_ch")->SetTitle("Jet Charge");
  mtree->Branch("jet_mass", &jet_mass);
  mtree->GetBranch("jet_mass")->SetTitle("Jet Mass");
  mtree->Branch("jet_btag", &jet_btag);
  mtree->GetBranch("jet_btag")->SetTitle("Jet Btagging Discriminant (CSV)");
  mtree->Branch("jet_flavour", &jet_flavour);
  mtree->GetBranch("jet_flavour")->SetTitle("Jet Parton Flavour");
}

SimpleJetAnalyzer::~SimpleJetAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void SimpleJetAnalyzer::analyze(const edm::Event &iEvent,
                                const edm::EventSetup &iSetup) {

  using namespace edm;
  using namespace std;

  Handle<reco::PFJetCollection> myjets;
  iEvent.getByLabel(jetInput, myjets);
  Handle<reco::JetTagCollection> btags;
  iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);
  Handle<double> rhoHandle;
  iEvent.getByLabel(InputTag("kt6PFJets:rho"), rhoHandle);
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  Handle<reco::JetFlavourInfoMatchingCollection> injets;
  if (!isData) {
    iEvent.getByLabel(InputTag("jetFlavourInfosAK5PFJets"), injets);
  }

  numjet = 0;
  jet_e.clear();
  jet_pt.clear();
  jet_px.clear();
  jet_py.clear();
  jet_pz.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_ch.clear();
  jet_mass.clear();
  jet_btag.clear();
  jet_flavour.clear();

  if (myjets.isValid()) {

    int hadronFlavour;

    for (reco::PFJetCollection::const_iterator itjet = myjets->begin();
         itjet != myjets->end(); ++itjet) {
      if (!isData) {
        if (itjet->pt() >= min_pt) {

          jet_e.push_back(itjet->energy());
          jet_pt.push_back(itjet->pt());
          jet_px.push_back(itjet->px());
          jet_py.push_back(itjet->py());
          jet_pz.push_back(itjet->pz());
          jet_eta.push_back(itjet->eta());
          jet_phi.push_back(itjet->phi());
          jet_ch.push_back(itjet->charge());
          jet_mass.push_back(itjet->mass());
          if (btags.isValid() && (itjet - myjets->begin()) < btags->size()) {
            jet_btag.push_back(
                btags->operator[](itjet - myjets->begin()).second);
          } else
            jet_btag.push_back(-999);

          reco::JetFlavourInfo aInfo =
              injets->operator[](itjet - myjets->begin()).second;
          hadronFlavour = aInfo.getPartonFlavour();
          jet_flavour.push_back(hadronFlavour);

          ++numjet;
        }
      }
    }
  }

  mtree->Fill();
  return;
}

// ------------ method called once each job just before starting event loop
// ------------
void SimpleJetAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void SimpleJetAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------
void SimpleJetAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run ------------
void SimpleJetAnalyzer::endRun(edm::Run const &, edm::EventSetup const &) {}
// ------------ method called when starting to processes a luminosity block
// ------------
void SimpleJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &,
                                             edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block
// ------------
void SimpleJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &,
                                           edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for
// the module  ------------
void SimpleJetAnalyzer::fillDescriptions(
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
DEFINE_FWK_MODULE(SimpleJetAnalyzer);
