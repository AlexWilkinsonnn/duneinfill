////////////////////////////////////////////////////////////////////////
// Class:       ExportDigits
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExportDigits_module.cc
//
// Copied on 28 Apr 2022.
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

#include <memory>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace Infill {
  class ExportDigits;
}

class Infill::ExportDigits : public art::EDAnalyzer {
public:
  explicit ExportDigits(fhicl::ParameterSet const& p);

  ExportDigits(ExportDigits const&) = delete;
  ExportDigits(ExportDigits&&) = delete;
  ExportDigits& operator=(ExportDigits const&) = delete;
  ExportDigits& operator=(ExportDigits&&) = delete;

  struct TPCPlaneInfo;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();

private:
  const geo::GeometryCore* fGeom;

  TTree*                        fTreeDigits;
  // Need to use int instead of short to avoid providing a definition file for ROOT
  std::vector<std::vector<int>> fDigitsZ;
  std::vector<std::vector<int>> fDigitsU;
  std::vector<std::vector<int>> fDigitsV;

  unsigned int              fCIndex;
  std::vector<unsigned int> fTIndices;
  std::vector<TPCPlaneInfo> fTIDPlanesVec;

  std::string fDigitsLabel;
};


struct Infill::ExportDigits::TPCPlaneInfo {
  TPCPlaneInfo(geo::TPCID _tID, const geo::GeometryCore* geom)
    : tID(_tID)
  {
    for (geo::PlaneID pID : geom->IteratePlaneIDs(_tID)) {
      if (geom->View(pID) == geo::kZ) {
        std::cout << "Z plane: " << pID << "\n";
        rIDZ = geom->WirePlaneToROP(pID);
      }
      else if (geom->View(pID) == geo::kU) {
        std::cout << "U plane: " << pID << "\n";
        rIDU = geom->WirePlaneToROP(pID);
      }
      else if (geom->View(pID) == geo::kV) {
        std::cout << "V plane: " << pID << "\n";
        rIDV = geom->WirePlaneToROP(pID);
      }
    }
  }

  geo::TPCID     tID;
  readout::ROPID rIDZ;
  readout::ROPID rIDU;
  readout::ROPID rIDV;
};

Infill::ExportDigits::ExportDigits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fCIndex      (p.get<unsigned int>("CryoIndex")),
    fTIndices    (p.get<std::vector<unsigned int>>("TpcIndices")),
    fDigitsLabel (p.get<std::string>("DigitsLabel"))
{
  consumes<std::vector<raw::RawDigit>>(fDigitsLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeDigits = tfs->make<TTree>("digits", "digits");
  fTreeDigits->Branch("digit_vecsZ", &fDigitsZ);
  fTreeDigits->Branch("digit_vecsU", &fDigitsU);
  fTreeDigits->Branch("digit_vecsV", &fDigitsV);
}

void Infill::ExportDigits::analyze(art::Event const& e)
{
  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel(fDigitsLabel, digs);

  for (auto tIDPlanes : fTIDPlanesVec) {
    this->reset();

    readout::ROPID rIDZ = tIDPlanes.rIDZ;
    readout::ROPID rIDU = tIDPlanes.rIDU;
    readout::ROPID rIDV = tIDPlanes.rIDV;

    fDigitsZ = std::vector<std::vector<int>>(480);
    fDigitsU = std::vector<std::vector<int>>(800);
    fDigitsV = std::vector<std::vector<int>>(800);
    raw::ChannelID_t firstChZ = fGeom->FirstChannelInROP(rIDZ);
    raw::ChannelID_t firstChU = fGeom->FirstChannelInROP(rIDU);
    raw::ChannelID_t firstChV = fGeom->FirstChannelInROP(rIDV);

    for (const raw::RawDigit& dig : *digs) {
      if (fGeom->ChannelToROP(dig.Channel()) == rIDZ) {
        raw::RawDigit::ADCvector_t adcs(dig.Samples());
        raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

        fDigitsZ[dig.Channel() - firstChZ] = std::vector<int>(6000);
        for (unsigned int tick = 0; tick < 6000; tick++) {
          const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

          fDigitsZ[dig.Channel() - firstChZ][tick] = (int)adc;
        }
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == rIDU) {
        raw::RawDigit::ADCvector_t adcs(dig.Samples());
        raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

        fDigitsU[dig.Channel() - firstChU] = std::vector<int>(6000);
        for (unsigned int tick = 0; tick < 6000; tick++) {
          const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

          fDigitsU[dig.Channel() - firstChU][tick] = (int)adc;
        }
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == rIDV) {
        raw::RawDigit::ADCvector_t adcs(dig.Samples());
        raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

        fDigitsV[dig.Channel() - firstChV] = std::vector<int>(6000);
        for (unsigned int tick = 0; tick < 6000; tick++) {
          const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

          fDigitsV[dig.Channel() - firstChV][tick] = (int)adc;
        }
      }
    }

    fTreeDigits->Fill();
  }
}

void Infill::ExportDigits::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  const geo::CryostatID cID(fCIndex);
  for (auto tIndex : fTIndices) { 
    std::cout << "TPC " << tIndex << "\n";
    const geo::TPCID tID(cID, tIndex);
    fTIDPlanesVec.push_back(TPCPlaneInfo(tID, fGeom));
  }
}

void Infill::ExportDigits::endJob()
{
}

void Infill::ExportDigits::reset()
{
  fDigitsZ.clear();
  fDigitsU.clear();
  fDigitsV.clear();
}

DEFINE_ART_MODULE(Infill::ExportDigits)

