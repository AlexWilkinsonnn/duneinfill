////////////////////////////////////////////////////////////////////////
// Class:       RecoAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        RecoAna_module.cc
//
// Generated at Thu Feb 10 15:31:48 2022 by Alex Wilkinson using cetskelgen
// from  version .
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

#include <map>
#include <string>


namespace Infill {
  class RecoAna;
}


class Infill::RecoAna : public art::EDAnalyzer {
public:
  explicit RecoAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RecoAna(RecoAna const&) = delete;
  RecoAna(RecoAna&&) = delete;
  RecoAna& operator=(RecoAna const&) = delete;
  RecoAna& operator=(RecoAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

};


Infill::RecoAna::RecoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{
  // consumes<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChPerfect"));
  // consumes<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChnov2019"));

  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChPerfect"));
  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019"));
  consumes<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019Infill"));

  consumes<std::vector<recob::PFParticle>>(art::InputTag("pandora", "", "RecoChPerfect"));
  consumes<std::vector<recob::PFParticle>>(art::InputTag("pandora", "", "RecoChnov2019"));

  // consumes<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChPerfect"));
  // consumes<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChnov2019"));

  consumes<std::vector<recob::Track>>(art::InputTag("pandoraTrack", "", "RecoChPerfect"));
  consumes<std::vector<recob::Track>>(art::InputTag("pandoraTrack", "", "RecoChnov2019"));

  // consumes<std::vector<anab::Calorimetry>>(art::InputTag("pandoracalo", "", "RecoChPerfect"));
  // consumes<std::vector<anab::Calorimetry>>(art::InputTag("pandoracalo", "", "RecoChnov2019"));

  // consumes<std::vector<anab::Calorimetry>>(art::InputTag("pandoraShowercalo", "", "RecoChPerfect"));
  // consumes<std::vector<anab::Calorimetry>>(art::InputTag("pandoraShowercalo", "", "RecoChnov2019"));

  consumes<std::vector<recob::Shower>>(art::InputTag("pandoraShower", "", "RecoChPerfect"));
  consumes<std::vector<recob::Shower>>(art::InputTag("pandoraShower", "", "RecoChnov2019"));

  consumes<std::vector<simb::MCTruth>>(art::InputTag("generator", "", "SinglesGen"));
}

void Infill::RecoAna::analyze(art::Event const& e)
{
  const auto pIDsPerfect = e.getValidHandle<std::vector<recob::PFParticle>>(art::InputTag("pandora", "", "RecoChPerfect"));
  const auto pIDsReal = e.getValidHandle<std::vector<recob::PFParticle>>(art::InputTag("pandora", "", "RecoChnov2019"));

  std::map<int, int> pIDCntrPerfect;
  std::map<int, int> pIDCntrReal;
  for (auto pID : *pIDsPerfect) {
    pIDCntrPerfect[pID.PdgCode()]++;
  }
  for (auto pID : *pIDsReal) {
    pIDCntrReal[pID.PdgCode()]++;
  }
  
  for (auto pIDCnt : pIDCntrPerfect) {
    std::cout << pIDCnt.first << ": " << pIDCnt.second << "\n";
  }
  std::cout << "---\n";
  for (auto pIDCnt : pIDCntrReal) {
    std::cout << pIDCnt.first << ": " << pIDCnt.second << "\n";
  }
  std::cout << "\n\n";

  const auto tracksPerfect = e.getValidHandle<std::vector<recob::Track>>(art::InputTag("pandoraTrack", "", "RecoChPerfect"));
  const auto tracksReal = e.getValidHandle<std::vector<recob::Track>>(art::InputTag("pandoraTrack", "", "RecoChnov2019"));

  std::map<int, int> trackLengthCntrPerfect;
  std::map<int, int> trackLengthCntrReal;
  for(auto track : *tracksPerfect) {
    trackLengthCntrPerfect[(int)(track.Length() + 0.5)/10]++;
  }
  for(auto track : *tracksReal) {
    trackLengthCntrReal[(int)(track.Length() + 0.5)/10]++;
  }

  for (auto trackLengthCnt : trackLengthCntrPerfect) {
    std::cout << trackLengthCnt.first << ": " << trackLengthCnt.second << "\n";
  }
  std::cout << "---\n";
  for (auto trackLengthCnt : trackLengthCntrReal) {
    std::cout << trackLengthCnt.first << ": " << trackLengthCnt.second << "\n";
  } 
  std::cout << "\n\n";

  const auto showersPerfect = e.getValidHandle<std::vector<recob::Shower>>(art::InputTag("pandoraShower", "", "RecoChPerfect"));
  const auto showersReal = e.getValidHandle<std::vector<recob::Shower>>(art::InputTag("pandoraShower", "", "RecoChnov2019"));

  std::map<float, int> showerOpenAngleCntrPerfect;
  std::map<float, int> showerOpenAngleCntrReal;
  for (auto shower : *showersPerfect) {
    showerOpenAngleCntrPerfect[(float)(((int)((shower.OpenAngle() + 0.05)*10)))/10]++;
  }
  for (auto shower : *showersReal) {
    showerOpenAngleCntrReal[(float)(((int)((shower.OpenAngle() + 0.05)*10)))/10]++;
  }

  for (auto showerOpenAngleCnt : showerOpenAngleCntrPerfect) {
    std::cout << showerOpenAngleCnt.first << ": " << showerOpenAngleCnt.second << "\n";
  }
  std::cout << "---\n";
  for (auto showerOpenAngleCnt : showerOpenAngleCntrReal) {
    std::cout << showerOpenAngleCnt.first << ": " << showerOpenAngleCnt.second << "\n";
  } 
  std::cout << "\n\n"; 

  std::map<int, int> showerLengthCntrPerfect;
  std::map<int, int> showerLengthCntrReal;
  for (auto shower : *showersPerfect) {
    showerLengthCntrPerfect[(int)(shower.Length() + 0.5)/10]++;
  }
  for (auto shower : *showersReal) {
    showerLengthCntrReal[(int)(shower.Length() + 0.5)/10]++;
  }

  for (auto showerLengthCnt : showerLengthCntrPerfect) {
    std::cout << showerLengthCnt.first << ": " << showerLengthCnt.second << "\n";
  }
  std::cout << "---\n";
  for (auto showerLengthCnt : showerLengthCntrReal) {
    std::cout << showerLengthCnt.first << ": " << showerLengthCnt.second << "\n";
  } 
  std::cout << "\n\n"; 

  const auto particlePIDsPerfect = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChPerfect"));
  const auto particlePIDsReal = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019"));
  const auto particlePIDsInfill = e.getValidHandle<std::vector<anab::ParticleID>>(art::InputTag("pandorapid", "", "RecoChnov2019Infill"));

  std::map<int, int> particlePIDBestCntrPerfect;
  std::map<int, int> particlePIDBestCntrReal;
  std::map<int, int> particlePIDBestCntrInfill;
  for (auto pID : *particlePIDsPerfect) {
    if (pID.ParticleIDAlgScores().size() != 0) {
      int pdg = 0;
      int score = 0;
      for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
        if (pIDAlgScore.fAlgName == "Chi2") {
          if (pIDAlgScore.fValue > score) {
            score = pIDAlgScore.fValue;
            pdg = pIDAlgScore.fAssumedPdg;
          }
        }
      }
      particlePIDBestCntrPerfect[pdg]++;
    }
  }
  for (auto pID : *particlePIDsReal) {
    if (pID.ParticleIDAlgScores().size() != 0) {
      int pdg = 0;
      int score = 0;
      for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
        if (pIDAlgScore.fAlgName == "Chi2") {
          if (pIDAlgScore.fValue > score) {
            score = pIDAlgScore.fValue;
            pdg = pIDAlgScore.fAssumedPdg;
          }
        }
      }
      particlePIDBestCntrReal[pdg]++;
    }
  }
  for (auto pID : *particlePIDsInfill) {
    if (pID.ParticleIDAlgScores().size() != 0) {
      int pdg = 0;
      int score = 0;
      for (auto pIDAlgScore : pID.ParticleIDAlgScores()) {
        if (pIDAlgScore.fAlgName == "Chi2") {
          if (pIDAlgScore.fValue > score) {
            score = pIDAlgScore.fValue;
            pdg = pIDAlgScore.fAssumedPdg;
          }
        }
      }
      particlePIDBestCntrInfill[pdg]++;
    }
  }

  std::cout << "--Perfect--\n";
  for (auto particlePIDBestCnt : particlePIDBestCntrPerfect) {
    std::cout << particlePIDBestCnt.first << ": " << particlePIDBestCnt.second << "\n";
  }
  std::cout << "--Real--\n";
  for (auto particlePIDBestCnt : particlePIDBestCntrReal) {
    std::cout << particlePIDBestCnt.first << ": " << particlePIDBestCnt.second << "\n";
  } 
  std::cout << "--Infill--\n";
  for (auto particlePIDBestCnt : particlePIDBestCntrInfill) {
    std::cout << particlePIDBestCnt.first << ": " << particlePIDBestCnt.second << "\n";
  } 
  std::cout << "--True--\n";

  const auto trueGenParticles = e.getValidHandle<std::vector<simb::MCTruth>>(art::InputTag("generator", "", "SinglesGen"));

  std::map<int, int> trueParticleCntr;
  for (auto truth : *trueGenParticles) {
    for (int i = 0; i < truth.NParticles(); i++) { 
      const simb::MCParticle truthParticle = truth.GetParticle(i);
      trueParticleCntr[truthParticle.PdgCode()]++;
    }
  }

  for (auto trueParticleCnt : trueParticleCntr) {
    std::cout << trueParticleCnt.first << ": " << trueParticleCnt.second << "\n";
  }
}

void Infill::RecoAna::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();
}

void Infill::RecoAna::endJob()
{
}

DEFINE_ART_MODULE(Infill::RecoAna)
