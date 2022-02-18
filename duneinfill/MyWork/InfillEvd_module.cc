////////////////////////////////////////////////////////////////////////
// Class:       InfillEvd
// Plugin Type: analyzer (Unknown Unknown)
// File:        InfillEvd_module.cc
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
#include "lardataobj/RecoBase/Wire.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include <string>
#include <vector>
#include <map>
#include "TH2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace Infill {
  class InfillEvd;
}


class Infill::InfillEvd : public art::EDAnalyzer {
public:
  explicit InfillEvd(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  InfillEvd(InfillEvd const&) = delete;
  InfillEvd(InfillEvd&&) = delete;
  InfillEvd& operator=(InfillEvd const&) = delete;
  InfillEvd& operator=(InfillEvd&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // My functions.
  void reset();
  void fillRawDigit(TH2D*, const raw::RawDigit&, const raw::ChannelID_t&);
  void fillWire(TH2D*, const recob::Wire&, const raw::ChannelID_t&);
  void fillHit(TH2D*, const recob::Hit&, const raw::ChannelID_t&);

private:
  const geo::GeometryCore* fGeom;

  std::set<raw::ChannelID_t> fBadChs;
  readout::ROPID             fRIDZ;
  readout::ROPID             fRIDUV;
  raw::ChannelID_t           fFirstChZ;
  raw::ChannelID_t           fFirstChUV;

  TH2D*  fHistWiresPerfectZ;
  TH2D*  fHistWiresPerfectUV;
  TH2D*  fHistWiresnov2019Z;
  TH2D*  fHistWiresnov2019UV;
  TH2D*  fHistWiresnov2019InfillZ;
  TH2D*  fHistWiresnov2019InfillUV;
  TH2D*  fHistWiresCaldataPerfectZ;
  TH2D*  fHistWiresCaldataPerfectUV;
  TH2D*  fHistWiresCaldatanov2019Z;
  TH2D*  fHistWiresCaldatanov2019UV;
  TH2D*  fHistWiresCaldatanov2019InfillZ;
  TH2D*  fHistWiresCaldatanov2019InfillUV;
  TH2D*  fHistHitsPerfectZ;
  TH2D*  fHistHitsPerfectUV;
  TH2D*  fHistHitsnov2019Z;
  TH2D*  fHistHitsnov2019UV;
  TH2D*  fHistHitsnov2019InfillZ;
  TH2D*  fHistHitsnov2019InfillUV;
  TH2D*  fHistDigsPerfectZ;
  TH2D*  fHistDigsPerfectUV;
  TH2D*  fHistDigsnov2019Z;
  TH2D*  fHistDigsnov2019UV;
  TH2D*  fHistDigsnov2019InfillZ;
  TH2D*  fHistDigsnov2019InfillUV;

  TTree*                    fTreeBadChs;
  int                       fROP;
  std::vector<unsigned int> fROPBadChs;
};


Infill::InfillEvd::InfillEvd(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChPerfect"));
  consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChnov2019"));
  consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChnov2019Infill"));

  consumes<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChPerfect"));
  consumes<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChnov2019"));
  consumes<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChnov2019Infill"));

  consumes<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChPerfect"));
  consumes<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChnov2019"));
  consumes<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChnov2019Infill"));

  consumes<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq", "DetsimStage1"));
  consumes<std::vector<raw::RawDigit>>(art::InputTag("infill", "", "InfillChannelsPD"));
}

void Infill::InfillEvd::analyze(art::Event const& e)
{
  // Write RawDigits data to TH2s.
  const auto digsPerfect =       e.getValidHandle<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq", "DetsimStage1"));
  const auto digsnov2019Infill = e.getValidHandle<std::vector<raw::RawDigit>>(art::InputTag("infill", "", "InfillChannelsPD"));

  for (const raw::RawDigit &dig : *digsPerfect) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      this->fillRawDigit(fHistDigsPerfectZ, dig, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
      this->fillRawDigit(fHistDigsPerfectUV, dig, fFirstChUV);
    }
  }

  for (const raw::RawDigit &dig : *digsPerfect) {
    if (!fBadChs.count(dig.Channel())) {
      if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
        this->fillRawDigit(fHistDigsnov2019Z, dig, fFirstChZ);
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
        this->fillRawDigit(fHistDigsnov2019UV, dig, fFirstChUV);
      }
    }
  }

  for (const raw::RawDigit &dig : *digsnov2019Infill) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      this->fillRawDigit(fHistDigsnov2019InfillZ, dig, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
      this->fillRawDigit(fHistDigsnov2019InfillUV, dig, fFirstChUV);
    }
  }

  std::cout << "RawDigits written.\n";

  // Write Wire data to TH2s.
  const auto wiresPerfect =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChPerfect"));
  const auto wiresnov2019 =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChnov2019"));
  const auto wiresnov2019Infill = e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("wclsdatasp", "gauss", "RecoChnov2019Infill"));

  for (const recob::Wire &wire : *wiresPerfect) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresPerfectZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresPerfectUV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresnov2019) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresnov2019Z, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresnov2019UV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresnov2019Infill) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresnov2019InfillZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresnov2019InfillUV, wire, fFirstChUV);
    }
  }

  std::cout << "Wires written.\n";

  // Write caldata Wire data to TH2s (want to know if wirecell has bad channels hardcoded).
  const auto wiresCaldataPerfect =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChPerfect"));
  const auto wiresCaldatanov2019 =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChnov2019"));
  const auto wiresCaldatanov2019Infill = e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChnov2019Infill"));

  for (const recob::Wire &wire : *wiresCaldataPerfect) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresCaldataPerfectZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresCaldataPerfectUV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresCaldatanov2019) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresCaldatanov2019Z, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresCaldatanov2019UV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresCaldatanov2019Infill) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresCaldatanov2019InfillZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresCaldatanov2019InfillUV, wire, fFirstChUV);
    }
  }

  std::cout << "caldata wires written.\n";

  // Write Hit data to TH2s
  const auto hitsPerfect =       e.getValidHandle<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChPerfect"));
  const auto hitsnov2019 =       e.getValidHandle<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChnov2019"));
  const auto hitsnov2019Infill = e.getValidHandle<std::vector<recob::Hit>>(art::InputTag("gaushit", "", "RecoChnov2019Infill"));

  for (const recob::Hit &hit : *hitsPerfect) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsPerfectZ, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsPerfectUV, hit, fFirstChUV);
    }
  }

  for (const recob::Hit &hit : *hitsnov2019) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsnov2019Z, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsnov2019UV, hit, fFirstChUV);
    }
  }

  for (const recob::Hit &hit : *hitsnov2019Infill) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsnov2019InfillZ, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsnov2019InfillUV, hit, fFirstChUV);
    }
  }

  std::cout << "Hits written.\n";
}

void Infill::InfillEvd::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  art::ServiceHandle<art::TFileService> tfs;

  fHistWiresPerfectZ =        tfs->make<TH2D>("wclsdatasp_gauss_RecoChPerfectZ", "Perfect WCSP Gauss;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresPerfectUV =       tfs->make<TH2D>("wclsdatasp_gauss_RecoChPerfectUV", "Perfect WCSP Gauss;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistWiresnov2019Z =        tfs->make<TH2D>("wclsdatasp_gauss_RecoChnov2019Z", "nov2019 WCSP Gauss;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresnov2019UV =       tfs->make<TH2D>("wclsdatasp_gauss_RecoChnov2019UV", "nov2019 WCSP Gauss;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistWiresnov2019InfillZ =  tfs->make<TH2D>("wclsdatasp_gauss_RecoChnov2019InfillZ", "nov2019Infill WCSP Gauss;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresnov2019InfillUV = tfs->make<TH2D>("wclsdatasp_gauss_RecoChnov2019InfillUV", "nov2019Infill WCSP Gauss;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fHistWiresCaldataPerfectZ =        tfs->make<TH2D>("caldata_dataprep_RecoChPerfectZ", "Perfect dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresCaldataPerfectUV =       tfs->make<TH2D>("caldata_dataprep_RecoChPerfectUV", "Perfect dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistWiresCaldatanov2019Z =        tfs->make<TH2D>("caldata_dataprep_RecoChnov2019Z", "nov2019 dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresCaldatanov2019UV =       tfs->make<TH2D>("caldata_dataprep_RecoChnov2019UV", "nov2019 dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistWiresCaldatanov2019InfillZ =  tfs->make<TH2D>("caldata_dataprep_RecoChnov2019InfillZ", "nov2019Infill dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresCaldatanov2019InfillUV = tfs->make<TH2D>("caldata_dataprep_RecoChnov2019InfillUV", "nov2019Infill dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fHistHitsPerfectZ =        tfs->make<TH2D>("gaushit_RecoChPerfectZ", "Perfect Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsPerfectUV =       tfs->make<TH2D>("gaushit_RecoChPerfectUV", "Perfect Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistHitsnov2019Z =        tfs->make<TH2D>("gaushit_RecoChnov2019Z", "nov2019 Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsnov2019UV =       tfs->make<TH2D>("gaushit_RecoChnov2019UV", "nov2019 Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistHitsnov2019InfillZ =  tfs->make<TH2D>("gaushit_RecoChnov2019InfillZ", "nov2019Infill Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsnov2019InfillUV = tfs->make<TH2D>("gaushit_RecoChnov2019InfillUV", "nov2019Infill Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fHistDigsPerfectZ =        tfs->make<TH2D>("rawdigits_RecoChPerfectZ", "Perfect Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsPerfectUV =       tfs->make<TH2D>("rawdigits_RecoChPerfectUV", "Perfect Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistDigsnov2019Z =        tfs->make<TH2D>("rawdigits_RecoChnov2019Z", "nov2019 Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsnov2019UV =       tfs->make<TH2D>("rawdigits_RecoChnov2019UV", "nov2019 Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistDigsnov2019InfillZ =  tfs->make<TH2D>("rawdigits_RecoChnov2019InfillZ", "nov2019Infill Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsnov2019InfillUV = tfs->make<TH2D>("rawdigits_RecoChnov2019InfillUV", "nov2019Infill Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fTreeBadChs = tfs->make<TTree>("nov2019_badchannels", "nov2019 Bad Channels");
  fTreeBadChs->Branch("ROP", &fROP, "rop/I");
  fTreeBadChs->Branch("BadChannels", &fROPBadChs);

  // Map bad channels to ROP if ROP is not facing wall.
  fBadChs = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>> activeRIDBadChsZ;
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>> activeRIDBadChsUV;

  for (const raw::ChannelID_t &ch : fBadChs) {
    const readout::ROPID rID = fGeom->ChannelToROP(ch);

    for (const geo::TPCID &tID : fGeom->ROPtoTPCs(rID)) { 
      // Width is drift coordinate x. DetHalfWidth: Normal is 178.813, wall facing is 0.736875.
      if (fGeom->DetHalfWidth(tID) > 170) { // Not facing wall.
        if (fGeom->View(ch) == geo::kZ) {
          activeRIDBadChsZ[fGeom->ChannelToROP(ch)].push_back(ch);
        }
        else if (fGeom->View(ch) == geo::kU || fGeom->View(ch) == geo::kV) {
          activeRIDBadChsUV[fGeom->ChannelToROP(ch)].push_back(ch);
        }
      }
    }
  }

  // Get collection and induction ROPs with the most bad channels.
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>>::iterator mostBadChsZ =
    std::max_element(activeRIDBadChsZ.begin(), activeRIDBadChsZ.end(),
    [] (const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &a, const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &b)
    -> bool{ return a.second.size() < b.second.size(); });
  fRIDZ = mostBadChsZ->first;
  fFirstChZ = fGeom->FirstChannelInROP(fRIDZ);

  std::map<readout::ROPID, std::vector<raw::ChannelID_t>>::iterator mostBadChsUV =
    std::max_element(activeRIDBadChsUV.begin(), activeRIDBadChsUV.end(),
    [] (const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &a, const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &b)
    -> bool{ return a.second.size() < b.second.size(); });
  fRIDUV = mostBadChsUV->first;
  fFirstChUV = fGeom->FirstChannelInROP(fRIDUV);

  // Dump bad channel info
  this->reset();
  fROP = fRIDZ.ROP;
  for (auto &ch : mostBadChsZ->second) { 
    ch -= fFirstChZ;
  }
  fROPBadChs = mostBadChsZ->second;
  fTreeBadChs->Fill();

  this->reset();
  fROP = fRIDUV.ROP;
  for (auto &ch : mostBadChsUV->second) {
    ch -= fFirstChUV;
  }
  fROPBadChs = mostBadChsUV->second;
  fTreeBadChs->Fill();
}

void Infill::InfillEvd::endJob()
{
}

void Infill::InfillEvd::fillRawDigit(TH2D *hist, const raw::RawDigit &dig, const raw::ChannelID_t &firstCh)
{
  raw::RawDigit::ADCvector_t adcs(dig.Samples());
  raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

  for(unsigned int tick = 0; tick < adcs.size(); ++tick) {
    const int adc = adcs[tick] ? int(adcs[tick]) - dig.GetPedestal() : 0;

    hist->Fill(dig.Channel() - firstCh, tick, adc);    
  }
}

void Infill::InfillEvd::fillWire(TH2D *hist, const recob::Wire &wire, const raw::ChannelID_t &firstCh)
{
  for (unsigned int tick = 0; tick < wire.Signal().size(); ++tick) {
    hist->Fill(wire.Channel() - firstCh, tick, wire.Signal()[tick]);
  }
}

void Infill::InfillEvd::fillHit(TH2D *hist, const recob::Hit &hit, const raw::ChannelID_t &firstCh)
{
  for (int tick = hit.StartTick(); tick <= hit.EndTick(); ++tick) {
    hist->Fill(hit.Channel() - firstCh, tick, hit.SummedADC());
  }
}

void Infill::InfillEvd::reset()
{
  fROP = -1;
  fROPBadChs.clear();
}

DEFINE_ART_MODULE(Infill::InfillEvd)
