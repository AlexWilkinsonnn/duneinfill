#ifndef MY_PROTODUNE_TRACK_UTILS_H
#define MY_PROTODUNE_TRACK_UTILS_H

///////////////////////////////////////////////////////////////
// MyProtoDUNETrackUtils
//  - Class to help analysers access useful track information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
//
// Copied from protoduneana/protoduneana/Utilities on 3 Mar 2022.
// Removing hardcoded emtrkmichelid tags.
///////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "art/Framework/Principal/Event.h"
#include "TProfile.h"

namespace protoana {

  class MyProtoDUNETrackUtils {

  public:

    MyProtoDUNETrackUtils();
    ~MyProtoDUNETrackUtils();

  std::pair<double, int> GetVertexMichelScore(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    const std::string emtrkmichelidModule,
    double min_length = 5., double min_x = -200.,
    double max_x = 0., double min_y = 200., double max_y = 500.,
    double min_z = 25., bool check_wire = true, double check_x = 0,
    double check_y = 0., double check_z = 0.);
  std::pair<double, double> GetVertexMichelScore_weight_by_charge(
      const recob::Track & track, const art::Event & evt,
      const std::string trackModule, const std::string hitModule,
      const std::string emtrkmichelidModule,
      double min_length = 5., double min_x = -200.,
      double max_x = 0., double min_y = 200., double max_y = 500.,
      double min_z = 25., bool check_wire = true, double check_x = 0,
      double check_y = 0., double check_z = 0.);
  std::pair<double, int> GetVertexMichelScoreAlt(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    const std::string emtrkmichelidModule,
    double min_length = 5., double min_x = -200.,
    double max_x = 0., double min_y = 200., double max_y = 500.,
    double min_z = 25., bool check_wire = true, double check_x = 0,
    double check_y = 0., double check_z = 0.);

  private:


  };

}

#endif

