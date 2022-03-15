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

  std::set<raw::ChannelID_t> fDeadChs;
  readout::ROPID             fRIDZ;
  readout::ROPID             fRIDUV;
  raw::ChannelID_t           fFirstChZ;
  raw::ChannelID_t           fFirstChUV;

  TH2D*  fHistWiresPerfectZ;
  TH2D*  fHistWiresPerfectUV;
  TH2D*  fHistWiresDeadZ;
  TH2D*  fHistWiresDeadUV;
  TH2D*  fHistWiresDeadInfillZ;
  TH2D*  fHistWiresDeadInfillUV;
  // TH2D*  fHistWiresCaldataPerfectZ;
  // TH2D*  fHistWiresCaldataPerfectUV;
  // TH2D*  fHistWiresCaldataDeadZ;
  // TH2D*  fHistWiresCaldataDeadUV;
  // TH2D*  fHistWiresCaldataDeadInfillZ;
  // TH2D*  fHistWiresCaldataDeadInfillUV;
  TH2D*  fHistHitsPerfectZ;
  TH2D*  fHistHitsPerfectUV;
  TH2D*  fHistHitsDeadZ;
  TH2D*  fHistHitsDeadUV;
  TH2D*  fHistHitsDeadInfillZ;
  TH2D*  fHistHitsDeadInfillUV;
  TH2D*  fHistDigsPerfectZ;
  TH2D*  fHistDigsPerfectUV;
  TH2D*  fHistDigsDeadZ;
  TH2D*  fHistDigsDeadUV;
  TH2D*  fHistDigsDeadInfillZ;
  TH2D*  fHistDigsDeadInfillUV;

  TTree*                    fTreeDeadChs;
  int                       fROP;
  std::vector<unsigned int> fROPDeadChs;

  bool fBadChsInfilled;
  bool fNoisyChsInfilled;
  std::string fPerfectWCSPLabel;
  std::string fRealWCSPLabel;
  std::string fInfillWCSPLabel;
  std::string fPerfectHitLabel;
  std::string fRealHitLabel;
  std::string fInfillHitLabel;
  std::string fPerfectDigitLabel;
  std::string fInfillDigitLabel;
};


Infill::InfillEvd::InfillEvd(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fBadChsInfilled    (p.get<bool>        ("BadChsInfilled")),
    fNoisyChsInfilled  (p.get<bool>        ("NoisyChsInfilled")),
    fPerfectWCSPLabel  (p.get<std::string> ("PerfectWCSPLabel")),
    fRealWCSPLabel     (p.get<std::string> ("RealWCSPLabel")),
    fInfillWCSPLabel   (p.get<std::string> ("InfillWCSPLabel")),
    fPerfectHitLabel   (p.get<std::string> ("PerfectHitLabel")),
    fRealHitLabel      (p.get<std::string> ("RealHitLabel")),
    fInfillHitLabel    (p.get<std::string> ("InfillHitLabel")),
    fPerfectDigitLabel (p.get<std::string> ("PerfectDigitLabel")),
    fInfillDigitLabel  (p.get<std::string> ("InfillDigitLabel"))
{
  // consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChPerfect"));
  // consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChDead"));
  // consumes<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChDeadInfill"));

  consumes<std::vector<recob::Wire>>(fPerfectWCSPLabel);
  consumes<std::vector<recob::Wire>>(fRealWCSPLabel);
  consumes<std::vector<recob::Wire>>(fInfillWCSPLabel);

  consumes<std::vector<recob::Hit>>(fPerfectHitLabel);
  consumes<std::vector<recob::Hit>>(fRealHitLabel);
  consumes<std::vector<recob::Hit>>(fInfillHitLabel);

  consumes<std::vector<raw::RawDigit>>(fPerfectDigitLabel);
  consumes<std::vector<raw::RawDigit>>(fInfillDigitLabel);
}

void Infill::InfillEvd::analyze(art::Event const& e)
{
  // Write RawDigits data to TH2s.
  const auto digsPerfect =       e.getValidHandle<std::vector<raw::RawDigit>>(fPerfectDigitLabel);
  const auto digsDeadInfill = e.getValidHandle<std::vector<raw::RawDigit>>(fInfillDigitLabel);

  for (const raw::RawDigit &dig : *digsPerfect) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      this->fillRawDigit(fHistDigsPerfectZ, dig, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
      this->fillRawDigit(fHistDigsPerfectUV, dig, fFirstChUV);
    }
  }

  for (const raw::RawDigit &dig : *digsPerfect) {
    if (!fDeadChs.count(dig.Channel())) {
      if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
        this->fillRawDigit(fHistDigsDeadZ, dig, fFirstChZ);
      }
      else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
        this->fillRawDigit(fHistDigsDeadUV, dig, fFirstChUV);
      }
    }
  }

  for (const raw::RawDigit &dig : *digsDeadInfill) {
    if (fGeom->ChannelToROP(dig.Channel()) == fRIDZ) {
      this->fillRawDigit(fHistDigsDeadInfillZ, dig, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(dig.Channel()) == fRIDUV) {
      this->fillRawDigit(fHistDigsDeadInfillUV, dig, fFirstChUV);
    }
  }

  std::cout << "RawDigits written.\n";

  // Write Wire data to TH2s.
  const auto wiresPerfect =       e.getValidHandle<std::vector<recob::Wire>>(fPerfectWCSPLabel);
  const auto wiresDead =       e.getValidHandle<std::vector<recob::Wire>>(fRealWCSPLabel);
  const auto wiresDeadInfill = e.getValidHandle<std::vector<recob::Wire>>(fInfillWCSPLabel);

  for (const recob::Wire &wire : *wiresPerfect) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresPerfectZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresPerfectUV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresDead) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresDeadZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresDeadUV, wire, fFirstChUV);
    }
  }

  for (const recob::Wire &wire : *wiresDeadInfill) {
    if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
      this->fillWire(fHistWiresDeadInfillZ, wire, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
      this->fillWire(fHistWiresDeadInfillUV, wire, fFirstChUV);
    }
  }

  std::cout << "Wires written.\n";

  // Write caldata Wire data to TH2s (want to know if wirecell has bad channels hardcoded).
  // const auto wiresCaldataPerfect =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChPerfect"));
  // const auto wiresCaldataDead =       e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChDead"));
  // const auto wiresCaldataDeadInfill = e.getValidHandle<std::vector<recob::Wire>>(art::InputTag("caldata", "dataprep", "RecoChDeadInfill"));

  // for (const recob::Wire &wire : *wiresCaldataPerfect) {
  //   if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
  //     this->fillWire(fHistWiresCaldataPerfectZ, wire, fFirstChZ);
  //   }
  //   else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
  //     this->fillWire(fHistWiresCaldataPerfectUV, wire, fFirstChUV);
  //   }
  // }

  // for (const recob::Wire &wire : *wiresCaldataDead) {
  //   if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
  //     this->fillWire(fHistWiresCaldataDeadZ, wire, fFirstChZ);
  //   }
  //   else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
  //     this->fillWire(fHistWiresCaldataDeadUV, wire, fFirstChUV);
  //   }
  // }

  // for (const recob::Wire &wire : *wiresCaldataDeadInfill) {
  //   if (fGeom->ChannelToROP(wire.Channel()) == fRIDZ) {
  //     this->fillWire(fHistWiresCaldataDeadInfillZ, wire, fFirstChZ);
  //   }
  //   else if (fGeom->ChannelToROP(wire.Channel()) == fRIDUV) {
  //     this->fillWire(fHistWiresCaldataDeadInfillUV, wire, fFirstChUV);
  //   }
  // }

  // std::cout << "caldata wires written.\n";

  // Write Hit data to TH2s
  const auto hitsPerfect =       e.getValidHandle<std::vector<recob::Hit>>(fPerfectHitLabel);
  const auto hitsDead =       e.getValidHandle<std::vector<recob::Hit>>(fRealHitLabel);
  const auto hitsDeadInfill = e.getValidHandle<std::vector<recob::Hit>>(fInfillHitLabel);

  for (const recob::Hit &hit : *hitsPerfect) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsPerfectZ, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsPerfectUV, hit, fFirstChUV);
    }
  }

  for (const recob::Hit &hit : *hitsDead) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsDeadZ, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsDeadUV, hit, fFirstChUV);
    }
  }

  for (const recob::Hit &hit : *hitsDeadInfill) {
    if (fGeom->ChannelToROP(hit.Channel()) == fRIDZ) {
      this->fillHit(fHistHitsDeadInfillZ, hit, fFirstChZ);
    }
    else if (fGeom->ChannelToROP(hit.Channel()) == fRIDUV) {
      this->fillHit(fHistHitsDeadInfillUV, hit, fFirstChUV);
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
  fHistWiresDeadZ =        tfs->make<TH2D>("wclsdatasp_gauss_RecoChDeadZ", "Dead WCSP Gauss;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresDeadUV =       tfs->make<TH2D>("wclsdatasp_gauss_RecoChDeadUV", "Dead WCSP Gauss;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistWiresDeadInfillZ =  tfs->make<TH2D>("wclsdatasp_gauss_RecoChDeadInfillZ", "DeadInfill WCSP Gauss;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistWiresDeadInfillUV = tfs->make<TH2D>("wclsdatasp_gauss_RecoChDeadInfillUV", "DeadInfill WCSP Gauss;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  // fHistWiresCaldataPerfectZ =        tfs->make<TH2D>("caldata_dataprep_RecoChPerfectZ", "Perfect dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  // fHistWiresCaldataPerfectUV =       tfs->make<TH2D>("caldata_dataprep_RecoChPerfectUV", "Perfect dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  // fHistWiresCaldataDeadZ =        tfs->make<TH2D>("caldata_dataprep_RecoChDeadZ", "Dead dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  // fHistWiresCaldataDeadUV =       tfs->make<TH2D>("caldata_dataprep_RecoChDeadUV", "Dead dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  // fHistWiresCaldataDeadInfillZ =  tfs->make<TH2D>("caldata_dataprep_RecoChDeadInfillZ", "DeadInfill dataprep;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  // fHistWiresCaldataDeadInfillUV = tfs->make<TH2D>("caldata_dataprep_RecoChDeadInfillUV", "DeadInfill dataprep;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fHistHitsPerfectZ =        tfs->make<TH2D>("gaushit_RecoChPerfectZ", "Perfect Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsPerfectUV =       tfs->make<TH2D>("gaushit_RecoChPerfectUV", "Perfect Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistHitsDeadZ =        tfs->make<TH2D>("gaushit_RecoChDeadZ", "Dead Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsDeadUV =       tfs->make<TH2D>("gaushit_RecoChDeadUV", "Dead Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistHitsDeadInfillZ =  tfs->make<TH2D>("gaushit_RecoChDeadInfillZ", "DeadInfill Gauss Hit;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistHitsDeadInfillUV = tfs->make<TH2D>("gaushit_RecoChDeadInfillUV", "DeadInfill Gauss Hit;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fHistDigsPerfectZ =        tfs->make<TH2D>("rawdigits_RecoChPerfectZ", "Perfect Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsPerfectUV =       tfs->make<TH2D>("rawdigits_RecoChPerfectUV", "Perfect Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistDigsDeadZ =        tfs->make<TH2D>("rawdigits_RecoChDeadZ", "Dead Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsDeadUV =       tfs->make<TH2D>("rawdigits_RecoChDeadUV", "Dead Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);
  fHistDigsDeadInfillZ =  tfs->make<TH2D>("rawdigits_RecoChDeadInfillZ", "DeadInfill Raw Digits;Channel;Tick", 480, 0, 480, 6000, 0, 6000);
  fHistDigsDeadInfillUV = tfs->make<TH2D>("rawdigits_RecoChDeadInfillUV", "DeadInfill Raw Digits;Channel;Tick", 800, 0, 800, 6000, 0, 6000);

  fTreeDeadChs = tfs->make<TTree>("Dead_deadchannels", "Dead Dead Channels");
  fTreeDeadChs->Branch("ROP", &fROP, "rop/I");
  fTreeDeadChs->Branch("DeadChannels", &fROPDeadChs);

  // Get channels that have been infilled.
  const std::set<raw::ChannelID_t> badChs =   art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  const std::set<raw::ChannelID_t> noisyChs = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().NoisyChannels();

  if (fBadChsInfilled && fNoisyChsInfilled) {
    std::merge(
      badChs.begin(), badChs.end(), noisyChs.begin(), noisyChs.end(), 
      std::inserter(fDeadChs, fDeadChs.begin())
    );
  }
  else if (fBadChsInfilled) {
    fDeadChs.insert(badChs.begin(), badChs.end());
  }
  else if (fNoisyChsInfilled) {
    fDeadChs.insert(noisyChs.begin(), noisyChs.end());
  }

  // Map dead channels to ROP if ROP is not facing wall.
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>> activeRIDDeadChsZ;
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>> activeRIDDeadChsUV;

  for (const raw::ChannelID_t &ch : fDeadChs) {
    const readout::ROPID rID = fGeom->ChannelToROP(ch);
    // std::cout << ch - fGeom->FirstChannelInROP(rID) << " ";
    for (const geo::TPCID &tID : fGeom->ROPtoTPCs(rID)) { 
      // Width is drift coordinate x. DetHalfWidth: Normal is 178.813, wall facing is 0.736875.
      if (fGeom->DetHalfWidth(tID) > 170) { // Not facing wall.
        if (fGeom->View(ch) == geo::kZ) {
          activeRIDDeadChsZ[fGeom->ChannelToROP(ch)].push_back(ch);
        }
        else if (fGeom->View(ch) == geo::kU || fGeom->View(ch) == geo::kV) {
          activeRIDDeadChsUV[fGeom->ChannelToROP(ch)].push_back(ch);
        }
      }
    }
  }
  // std::cout << "\n";

  // Get collection and induction ROPs with the most dead channels.
  std::map<readout::ROPID, std::vector<raw::ChannelID_t>>::iterator mostDeadChsZ =
    std::max_element(activeRIDDeadChsZ.begin(), activeRIDDeadChsZ.end(),
    [] (const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &a, const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &b)
    -> bool{ return a.second.size() < b.second.size(); });
  fRIDZ = mostDeadChsZ->first;
  fFirstChZ = fGeom->FirstChannelInROP(fRIDZ);

  std::map<readout::ROPID, std::vector<raw::ChannelID_t>>::iterator mostDeadChsUV =
    std::max_element(activeRIDDeadChsUV.begin(), activeRIDDeadChsUV.end(),
    [] (const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &a, const std::pair<readout::ROPID, std::vector<raw::ChannelID_t>> &b)
    -> bool{ return a.second.size() < b.second.size(); });
  fRIDUV = mostDeadChsUV->first;
  fFirstChUV = fGeom->FirstChannelInROP(fRIDUV);

  // Dump dead channel info
  this->reset();
  fROP = fRIDZ.ROP;
  for (auto &ch : mostDeadChsZ->second) { 
    ch -= fFirstChZ;
  }
  fROPDeadChs = mostDeadChsZ->second;
  fTreeDeadChs->Fill();

  this->reset();
  fROP = fRIDUV.ROP;
  for (auto &ch : mostDeadChsUV->second) {
    ch -= fFirstChUV;
  }
  fROPDeadChs = mostDeadChsUV->second;
  fTreeDeadChs->Fill();
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
  fROPDeadChs.clear();
}

DEFINE_ART_MODULE(Infill::InfillEvd)
