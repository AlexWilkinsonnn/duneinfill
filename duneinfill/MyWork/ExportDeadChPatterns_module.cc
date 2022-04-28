////////////////////////////////////////////////////////////////////////
// Class:       Interactive
// Plugin Type: analyzer (Unknown Unknown)
// File:        Interactive_module.cc
//
// Copied 28 Apr 2022
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

#include <vector>
#include <set>
#include <map>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

namespace Infill {
  class ExportDeadChPatterns;
}

class Infill::ExportDeadChPatterns : public art::EDAnalyzer {
public:
  explicit ExportDeadChPatterns(fhicl::ParameterSet const& p);

  ExportDeadChPatterns(ExportDeadChPatterns const&) = delete;
  ExportDeadChPatterns(ExportDeadChPatterns&&) = delete;
  ExportDeadChPatterns& operator=(ExportDeadChPatterns const&) = delete;
  ExportDeadChPatterns& operator=(ExportDeadChPatterns&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:
  const geo::GeometryCore* fGeom;

  TTree*                        fTreeDeadChPatterns;
  std::vector<std::vector<int>> fDeadChPatternsZ;
  std::vector<std::vector<int>> fDeadChPatternsU;
  std::vector<std::vector<int>> fDeadChPatternsV;
};


Infill::ExportDeadChPatterns::ExportDeadChPatterns(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  art::ServiceHandle<art::TFileService> tfs;
  
  fTreeDeadChPatterns = tfs->make<TTree>("deadchs", "deadchs");
  fTreeDeadChPatterns->Branch("deadchsZ", &fDeadChPatternsZ);
  fTreeDeadChPatterns->Branch("deadchsU", &fDeadChPatternsU);
  fTreeDeadChPatterns->Branch("deadchsV", &fDeadChPatternsV);
}

void Infill::ExportDeadChPatterns::analyze(art::Event const& e)
{
}

void Infill::ExportDeadChPatterns::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const std::set<raw::ChannelID_t> badChs = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  const std::set<raw::ChannelID_t> noisyChs = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().NoisyChannels();

  std::set<raw::ChannelID_t> deadChs;
  deadChs.insert(badChs.begin(), badChs.end());
  deadChs.insert(noisyChs.begin(), noisyChs.end());

  std::cout << badChs.size() << " bad channels, " << noisyChs.size() << "noisy channels -> " << deadChs.size() << " dead\n";

  std::map<readout::ROPID, std::vector<int>> deadChPatterns;

  for (auto ch : deadChs) {
    deadChPatterns[fGeom->ChannelToROP(ch)].push_back(
      (int)(ch - fGeom->FirstChannelInROP(fGeom->ChannelToROP(ch))));
  }

  for (auto rIDDeadChs : deadChPatterns) {
    if (fGeom->View(rIDDeadChs.first) == geo::kZ) {
      fDeadChPatternsZ.push_back(rIDDeadChs.second);
    }
    else if (fGeom->View(rIDDeadChs.first) == geo::kU) {
      fDeadChPatternsU.push_back(rIDDeadChs.second);
    }
    else if (fGeom->View(rIDDeadChs.first) == geo::kV) {
      fDeadChPatternsV.push_back(rIDDeadChs.second);
    }
  } 

  fTreeDeadChPatterns->Fill();
}

void Infill::ExportDeadChPatterns::endJob()
{
}

DEFINE_ART_MODULE(Infill::ExportDeadChPatterns)
