////////////////////////////////////////////////////////////////////////
// Class:       InfillAnaTreePID
// Plugin Type: analyzer (Unknown Unknown)
// File:        InfillAnaTreePID_module.cc
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
  class InfillAnaTreePID;
}


class Infill::InfillAnaTreePID : public art::EDAnalyzer {
public:
  explicit InfillAnaTreePID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  InfillAnaTreePID(InfillAnaTreePID const&) = delete;
  InfillAnaTreePID(InfillAnaTreePID&&) = delete;
  InfillAnaTreePID& operator=(InfillAnaTreePID const&) = delete;
  InfillAnaTreePID& operator=(InfillAnaTreePID&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My function.
  void reset();

private:
  const geo::GeometryCore* fGeom;

  TTree*              fTree;
  int                 fRun;
  int                 fSubRun;
  int                 fEventNum;
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
};


Infill::InfillAnaTreePID::InfillAnaTreePID(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{
  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChPerfect"));
  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019"));
  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019Infill"));

  consumes<std::vector<simb::MCTruth>>(art::InputTag("generator", "", "SinglesGen"));
}

void Infill::InfillAnaTreePID::analyze(art::Event const& e)
{
  this->reset();

  fRun = e.id().run();
  fSubRun = e.id().subRun(); 
  fEventNum = e.id().event();
  std::cout << fRun << "," << fSubRun << "," << fEventNum << "\n";

  const auto particlePIDsPerfect = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChPerfect"));
  const auto particlePIDsReal = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019"));
  const auto particlePIDsInfill = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019Infill"));

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

  const auto trueGenParticles = e.getValidHandle<std::vector<simb::MCTruth>>(art::InputTag("generator", "", "SinglesGen"));

  for (auto truth : *trueGenParticles) {
    fTrueNumParticles += truth.NParticles();
    for (int i = 0; i < truth.NParticles(); i++) {
      const simb::MCParticle truthParticle = truth.GetParticle(i);
      fTrueParticles.push_back(truthParticle.PdgCode());
    }
  }

  fTree->Fill();
}

void Infill::InfillAnaTreePID::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("InfillAnaTreePID", "InfillAnaTreePID");
  fTree->Branch("Run", &fRun, "run/I");
  fTree->Branch("SubRun", &fSubRun, "subrun/I");
  fTree->Branch("EventNum", &fEventNum, "eventnum/I");
  fTree->Branch("PerfectPionChi2s", &fPerfectPionChi2s);
  fTree->Branch("PerfectMuonChi2s", &fPerfectMuonChi2s);
  fTree->Branch("PerfectKaonChi2s", &fPerfectKaonChi2s);
  fTree->Branch("PerfectProtonChi2s", &fPerfectProtonChi2s);
  fTree->Branch("RealPionChi2s", &fRealPionChi2s);
  fTree->Branch("RealMuonChi2s", &fRealMuonChi2s);
  fTree->Branch("RealKaonChi2s", &fRealKaonChi2s);
  fTree->Branch("RealProtonChi2s", &fRealProtonChi2s);
  fTree->Branch("InfillPionChi2s", &fInfillPionChi2s);
  fTree->Branch("InfillMuonChi2s", &fInfillMuonChi2s);
  fTree->Branch("InfillKaonChi2s", &fInfillKaonChi2s);
  fTree->Branch("InfillProtonChi2s", &fInfillProtonChi2s);
  fTree->Branch("PerfectNumPIDs", &fPerfectNumPIDs, "perfectnumpids/I");
  fTree->Branch("RealNumPIDs", &fRealNumPIDs, "realnumpids/I");
  fTree->Branch("InfillNumPIDs", &fInfillNumPIDs, "infillnumpids/I");
  fTree->Branch("TrueParticles", &fTrueParticles);
  fTree->Branch("TrueNumParticles", &fTrueNumParticles, "truenumparticles/I");
}

void Infill::InfillAnaTreePID::endJob()
{
}

void Infill::InfillAnaTreePID::reset() 
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
}

DEFINE_ART_MODULE(Infill::InfillAnaTreePID)
