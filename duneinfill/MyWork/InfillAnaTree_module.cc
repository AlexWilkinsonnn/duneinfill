////////////////////////////////////////////////////////////////////////
// Class:       InfillAnaTree
// Plugin Type: analyzer (Unknown Unknown)
// File:        InfillAnaTree_module.cc
//
// Wed Feb 16 2022 Alex Wilkinson
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Wire.h"

#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace Infill {
  class InfillAnaTree;
}


class Infill::InfillAnaTree : public art::EDAnalyzer {
public:
  explicit InfillAnaTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  InfillAnaTree(InfillAnaTree const&) = delete;
  InfillAnaTree(InfillAnaTree&&) = delete;
  InfillAnaTree& operator=(InfillAnaTree const&) = delete;
  InfillAnaTree& operator=(InfillAnaTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  std::string fPerfectParticleIDLabel;
  std::string fRealParticleIDLabel;
  std::string fInfillParticleIDLabel;
  std::string fPerfectPFParticleLabel;
  std::string fRealPFParticleLabel;
  std::string fInfillPFParticleLabel;
  std::string fTruthLabel;
  std::string fPerfectTrackLabel;
  std::string fRealTrackLabel;
  std::string fInfillTrackLabel;
  
  TTree*              fTreePID;
  std::vector<double> fPerfectMuonChi2s;
  std::vector<double> fPerfectPionChi2s;
  std::vector<double> fPerfectKaonChi2s;
  std::vector<double> fPerfectProtonChi2s;
  std::vector<double> fRealMuonChi2s;
  std::vector<double> fRealPionChi2s;
  std::vector<double> fRealKaonChi2s;
  std::vector<double> fRealProtonChi2s;
  std::vector<double> fInfillMuonChi2s;
  std::vector<double> fInfillPionChi2s;
  std::vector<double> fInfillKaonChi2s;
  std::vector<double> fInfillProtonChi2s;
  int                 fPerfectNumPIDs;
  int                 fRealNumPIDs;
  int                 fInfillNumPIDs;
  std::vector<int>    fTrueParticles;
  int                 fTrueNumParticles;

  TTree*           fTreePFParticle;
  std::vector<int> fPerfectPdgs;
  std::vector<int> fRealPdgs;
  std::vector<int> fInfillPdgs;
  std::vector<int> fPerfectNDaughters;
  std::vector<int> fRealNDaughters;
  std::vector<int> fInfillNDaughters;
  int              fPerfectNPFParticles;      
  int              fRealNPFParticles;
  int              fInfillNPFParticles;

  TTree*              fTreeTrack;
  std::vector<double> fPerfectLengths;
  std::vector<double> fRealLengths;
  std::vector<double> fInfillLengths;

  int fRun;
  int fSubRun;
  int fEventNum;
};


Infill::InfillAnaTree::InfillAnaTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fPerfectParticleIDLabel (p.get<std::string> ("PerfectParticleIDLabel")),
    fRealParticleIDLabel    (p.get<std::string> ("RealParticleIDLabel")),
    fInfillParticleIDLabel  (p.get<std::string> ("InfillParticleIDLabel")),
    fPerfectPFParticleLabel (p.get<std::string> ("PerfectPFParticleLabel")),
    fRealPFParticleLabel    (p.get<std::string> ("RealPFParticleLabel")),
    fInfillPFParticleLabel  (p.get<std::string> ("InfillPFParticleLabel")),
    fTruthLabel             (p.get<std::string> ("TruthLabel")),
    fPerfectTrackLabel      (p.get<std::string> ("PerfectTrackLabel")),
    fRealTrackLabel         (p.get<std::string> ("RealTrackLabel")),
    fInfillTrackLabel       (p.get<std::string> ("InfillTrackLabel"))
{
  consumes<std::vector<anab::ParticleID>>(fPerfectParticleIDLabel);
  consumes<std::vector<anab::ParticleID>>(fRealParticleIDLabel);
  consumes<std::vector<anab::ParticleID>>(fInfillParticleIDLabel);

  consumes<std::vector<simb::MCTruth>>(fTruthLabel);

  consumes<std::vector<recob::PFParticle>>(fPerfectPFParticleLabel);
  consumes<std::vector<recob::PFParticle>>(fRealPFParticleLabel);
  consumes<std::vector<recob::PFParticle>>(fInfillPFParticleLabel);

  consumes<std::vector<recob::Track>>(fPerfectTrackLabel);
  consumes<std::vector<recob::Track>>(fRealTrackLabel);
  consumes<std::vector<recob::Track>>(fInfillTrackLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreePID = tfs->make<TTree>("InfillAnaTree", "InfillAnaTree");
  fTreePID->Branch("Run", &fRun, "run/I");
  fTreePID->Branch("SubRun", &fSubRun, "subrun/I");
  fTreePID->Branch("EventNum", &fEventNum, "eventnum/I");
  fTreePID->Branch("PerfectPionChi2s", &fPerfectPionChi2s);
  fTreePID->Branch("PerfectMuonChi2s", &fPerfectMuonChi2s);
  fTreePID->Branch("PerfectKaonChi2s", &fPerfectKaonChi2s);
  fTreePID->Branch("PerfectProtonChi2s", &fPerfectProtonChi2s);
  fTreePID->Branch("RealPionChi2s", &fRealPionChi2s);
  fTreePID->Branch("RealMuonChi2s", &fRealMuonChi2s);
  fTreePID->Branch("RealKaonChi2s", &fRealKaonChi2s);
  fTreePID->Branch("RealProtonChi2s", &fRealProtonChi2s);
  fTreePID->Branch("InfillPionChi2s", &fInfillPionChi2s);
  fTreePID->Branch("InfillMuonChi2s", &fInfillMuonChi2s);
  fTreePID->Branch("InfillKaonChi2s", &fInfillKaonChi2s);
  fTreePID->Branch("InfillProtonChi2s", &fInfillProtonChi2s);
  fTreePID->Branch("PerfectNumPIDs", &fPerfectNumPIDs, "perfectnumpids/I");
  fTreePID->Branch("RealNumPIDs", &fRealNumPIDs, "realnumpids/I");
  fTreePID->Branch("InfillNumPIDs", &fInfillNumPIDs, "infillnumpids/I");
  fTreePID->Branch("TrueParticles", &fTrueParticles);
  fTreePID->Branch("TrueNumParticles", &fTrueNumParticles, "truenumparticles/I");

  fTreePFParticle = tfs->make<TTree>("InfillAnaTreePFParticle", "InfillAnaTreePFParticle");
  fTreePFParticle->Branch("Run", &fRun, "run/I");
  fTreePFParticle->Branch("SubRun", &fSubRun, "subrun/I");
  fTreePFParticle->Branch("EventNum", &fEventNum, "eventnum/I");
  fTreePFParticle->Branch("PerfectPdgs", &fPerfectPdgs);
  fTreePFParticle->Branch("RealPdgs", &fRealPdgs);
  fTreePFParticle->Branch("InfillPdgs", &fInfillPdgs);
  fTreePFParticle->Branch("PerfectNDaughters", &fPerfectNDaughters);
  fTreePFParticle->Branch("RealNDaughters", &fRealNDaughters);
  fTreePFParticle->Branch("InfillNDaughters", &fInfillNDaughters);
  fTreePFParticle->Branch("PerfectNPFParticles", &fPerfectNPFParticles, "perfectnpfparticles/I");
  fTreePFParticle->Branch("RealNPFParticles", &fRealNPFParticles, "realnpfparticles/I");
  fTreePFParticle->Branch("InfillNPFParticles", &fInfillNPFParticles, "infillnpfparticles/I");

  fTreeTrack = tfs->make<TTree>("InfillAnaTreeTrack", "InfillAnaTreeTrack");
  fTreeTrack->Branch("Run", &fRun, "run/I");
  fTreeTrack->Branch("SubRun", &fSubRun, "subrun/I");
  fTreeTrack->Branch("EventNum", &fEventNum, "eventnum/I");
  fTreeTrack->Branch("PerfectLengths", &fPerfectLengths);
  fTreeTrack->Branch("RealLengths", &fRealLengths);
  fTreeTrack->Branch("InfillLengths", &fInfillLengths);
}

void Infill::InfillAnaTree::analyze(art::Event const& e)
{
  this->reset();

  // Get event id info
  fRun = e.id().run();
  fSubRun = e.id().subRun(); 
  fEventNum = e.id().event();

  // Dump anab::ParticleID data
  const auto particlePIDsPerfect = e.getValidHandle<std::vector<anab::ParticleID>>(fPerfectParticleIDLabel);
  const auto particlePIDsReal = e.getValidHandle<std::vector<anab::ParticleID>>(fRealParticleIDLabel);
  const auto particlePIDsInfill = e.getValidHandle<std::vector<anab::ParticleID>>(fInfillParticleIDLabel);

  for (auto pID : *particlePIDsPerfect) { 
    for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
      if (pIDAlgScore.fAlgName == "Chi2") {
        if (pIDAlgScore.fAssumedPdg == 13) {
          fPerfectMuonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 321) {
          fPerfectKaonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 211) {
          fPerfectPionChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 2212) {
          fPerfectProtonChi2s.push_back(pIDAlgScore.fValue);
        }
      }
    }
  }

  for (auto pID : *particlePIDsReal) { 
    for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
      if (pIDAlgScore.fAlgName == "Chi2") {
        if (pIDAlgScore.fAssumedPdg == 13) {
          fRealMuonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 321) {
          fRealKaonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 211) {
          fRealPionChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 2212) {
          fRealProtonChi2s.push_back(pIDAlgScore.fValue);
        }
      }
    }
  }

  for (auto pID : *particlePIDsInfill) { 
    for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
      if (pIDAlgScore.fAlgName == "Chi2") {
        if (pIDAlgScore.fAssumedPdg == 13) {
          fInfillMuonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 321) {
          fInfillKaonChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 211) {
          fInfillPionChi2s.push_back(pIDAlgScore.fValue);
        }
        else if (pIDAlgScore.fAssumedPdg == 2212) {
          fInfillProtonChi2s.push_back(pIDAlgScore.fValue);
        }
      }
    }
  }

  fPerfectNumPIDs = particlePIDsPerfect->size();
  fRealNumPIDs = particlePIDsReal->size();
  fInfillNumPIDs = particlePIDsInfill->size();

  const auto trueGenParticles = e.getValidHandle<std::vector<simb::MCTruth>>(fTruthLabel);

  for (auto truth : *trueGenParticles) {
    fTrueNumParticles += truth.NParticles();
    for (int i = 0; i < truth.NParticles(); i++) {
      const simb::MCParticle truthParticle = truth.GetParticle(i);
      fTrueParticles.push_back(truthParticle.PdgCode());
    }
  }

  fTreePID->Fill();

  // Dump recob::PFparticle data
  const auto PFParticlesPerfect = e.getValidHandle<std::vector<recob::PFParticle>>(fPerfectPFParticleLabel);
  const auto PFParticlesReal = e.getValidHandle<std::vector<recob::PFParticle>>(fRealPFParticleLabel);
  const auto PFParticlesInfill = e.getValidHandle<std::vector<recob::PFParticle>>(fInfillPFParticleLabel);

  for (auto PFParticle : *PFParticlesPerfect) {
    fPerfectPdgs.push_back(PFParticle.PdgCode());
    fPerfectNDaughters.push_back(PFParticle.NumDaughters());
  }
  for (auto PFParticle : *PFParticlesReal) {
    fRealPdgs.push_back(PFParticle.PdgCode());
    fRealNDaughters.push_back(PFParticle.NumDaughters());
  }   
  for (auto PFParticle : *PFParticlesInfill) {
    fInfillPdgs.push_back(PFParticle.PdgCode());
    fInfillNDaughters.push_back(PFParticle.NumDaughters());
  }

  fPerfectNPFParticles = PFParticlesPerfect->size();
  fRealNPFParticles = PFParticlesReal->size();
  fInfillNPFParticles = PFParticlesInfill->size();

  fTreePFParticle->Fill();

  // Dump recob::Track data
  const auto TracksPerfect = e.getValidHandle<std::vector<recob::Track>>(fPerfectTrackLabel);
  const auto TracksReal = e.getValidHandle<std::vector<recob::Track>>(fRealTrackLabel);
  const auto TracksInfill = e.getValidHandle<std::vector<recob::Track>>(fInfillTrackLabel);

  for (auto track : *TracksPerfect) {
    fPerfectLengths.push_back(track.Length());
  }
  for (auto track : *TracksReal) {
    fRealLengths.push_back(track.Length());
  }
  for (auto track : *TracksInfill) {
    fInfillLengths.push_back(track.Length());
  }

  fTreeTrack->Fill();
}

void Infill::InfillAnaTree::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void Infill::InfillAnaTree::endJob()
{
}

void Infill::InfillAnaTree::reset() 
{
  fRun = -999;
  fSubRun = -999;
  fEventNum = -999;

  fPerfectPionChi2s.clear();
  fPerfectMuonChi2s.clear();
  fPerfectKaonChi2s.clear();
  fPerfectProtonChi2s.clear();
  fRealPionChi2s.clear();
  fRealMuonChi2s.clear();
  fRealKaonChi2s.clear();
  fRealProtonChi2s.clear();
  fInfillPionChi2s.clear();
  fInfillMuonChi2s.clear();
  fInfillKaonChi2s.clear();
  fInfillProtonChi2s.clear();
  fPerfectNumPIDs = 0;
  fRealNumPIDs = 0;
  fInfillNumPIDs = 0;
  fTrueParticles.clear();
  fTrueNumParticles = 0;

  fPerfectPdgs.clear();
  fRealPdgs.clear();
  fInfillPdgs.clear();
  fPerfectNDaughters.clear();
  fRealNDaughters.clear();
  fInfillNDaughters.clear();
  fPerfectNPFParticles = 0;
  fRealNPFParticles = 0;
  fInfillNPFParticles = 0;
}

DEFINE_ART_MODULE(Infill::InfillAnaTree)
