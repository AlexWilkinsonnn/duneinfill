#include "duneinfill/MyPduneUtils/MyProtoDUNETrackUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "TFile.h"
#include "TH1F.h"

#include <string>

protoana::MyProtoDUNETrackUtils::MyProtoDUNETrackUtils(){

}

protoana::MyProtoDUNETrackUtils::~MyProtoDUNETrackUtils(){

}

//Jake's implementation
std::pair<double, int> protoana::MyProtoDUNETrackUtils::GetVertexMichelScore(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    const std::string emtrkmichelidModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  art::ServiceHandle<geo::Geometry> geom;
  anab::MVAReader<recob::Hit, 4> hitResults(evt, emtrkmichelidModule);

  //Skip short tracks
  //Skip tracks that start/end outside of interesting volume
  auto start = track.Vertex();
  auto end = track.End();
  if ((TMath::Max(start.X(), end.X()) > max_x) ||
      (TMath::Min(start.X(), end.X()) < min_x) ||
      (TMath::Max(start.Y(), end.Y()) > max_y) ||
      (TMath::Min(start.Y(), end.Y()) < min_y) ||
      (TMath::Min(start.Z(), end.Z()) < min_z) ||
      (track.Length() < min_length)) {
    return {-1., 0};
  }


  //Get the hits from the TrackHitMetas and only for view 2 (collection plane)
  protoana::ProtoDUNETrackUtils trackUtil;
  std::map<size_t, const recob::Hit *>
      hits_from_traj = trackUtil.GetRecoHitsFromTrajPoints(track, evt, trackModule);
  std::vector<const recob::Hit *> hits_from_traj_view2;
  std::vector<size_t> index_from_traj_view2;

  for (auto it = hits_from_traj.begin(); it != hits_from_traj.end(); ++it) {
    if (it->second->View() != 2) continue;
    hits_from_traj_view2.push_back(it->second); 
    index_from_traj_view2.push_back(it->first);
  }

  //Find the vertex hit & info to compare to later
  double highest_z = -100.; 
  int vertex_tpc = -1;
  int vertex_wire = -1;
  float vertex_peak_time = -1.;

  if (check_z) {
    for (const auto * hit : hits_from_traj_view2) {
      double wire_z = geom->Wire(hit->WireID()).GetCenter().Z();
      if (wire_z > highest_z) {
        highest_z = wire_z;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  else {
    const recob::TrackTrajectory & traj = track.Trajectory();
    double highest_diff = -1.;
    for (size_t i = 0; i < hits_from_traj_view2.size(); ++i) {
      const recob::Hit * hit = hits_from_traj_view2[i];
      size_t traj_index = index_from_traj_view2[i];
      
      double traj_x = traj.LocationAtPoint(traj_index).X();
      double traj_y = traj.LocationAtPoint(traj_index).Y();
      double traj_z = traj.LocationAtPoint(traj_index).Z();

      double diff = sqrt(std::pow((traj_x - check_x), 2) + 
                         std::pow((traj_y - check_y), 2) + 
                         std::pow((traj_z - check_z), 2));
      if (diff > highest_diff) {
        highest_diff = diff;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  
  std::pair<double, int> results = {0., 0};

  //Go through all hits in the event.
  auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(hitModule);
  std::vector<art::Ptr<recob::Hit>> hit_vector;
  art::fill_ptr_vector(hit_vector, allHits);
  art::FindManyP<recob::Track> tracks_from_all_hits(allHits, evt, trackModule);
  for (size_t i = 0; i < hit_vector.size(); ++i) {

    //If this hit is in the trajectory hits vector, skip
    const recob::Hit * theHit = hit_vector[i].get();
    if (std::find(hits_from_traj_view2.begin(),
        hits_from_traj_view2.end(),
        theHit) != hits_from_traj_view2.end()) {
      continue;
    }

    //Skip hits that are outside of our TPC/plane or window of interest
    int wire = theHit->WireID().Wire;
    float peak_time = theHit->PeakTime();
    int tpc = theHit->WireID().TPC; 
    int plane = theHit->View();
    if ((abs(wire - vertex_wire) > 15) ||
        (abs(peak_time - vertex_peak_time) > 100.) ||
        (tpc != vertex_tpc) || (plane != 2)) {
      continue;
    }

    //It's ok if the hits don't come from any track
    //or if that track is the primary one, because sometimes the michel hits
    //are associated to it.
    //
    //It's not ok if the hits come from another track that is long
    //(i.e. an actual daughter). Skip these
    auto & tracks_from_hit = tracks_from_all_hits.at(hit_vector[i].key());
    if (!tracks_from_hit.empty() &&
        (tracks_from_hit[0].get()->ID() != track.ID()) &&
        (tracks_from_hit[0].get()->Length() > 25.))
      continue;

    //add up the CNN results 
    std::array<float, 4> cnn_out = hitResults.getOutput(hit_vector[i]);
    results.first += cnn_out[hitResults.getIndex("michel")];
    results.second += 1;
  }

  return results;
}
//Yinrui's modification
std::pair<double, double> protoana::MyProtoDUNETrackUtils::GetVertexMichelScore_weight_by_charge(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    const std::string emtrkmichelidModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  art::ServiceHandle<geo::Geometry> geom;
  anab::MVAReader<recob::Hit, 4> hitResults(evt, emtrkmichelidModule);

  //Skip short tracks
  //Skip tracks that start/end outside of interesting volume
  auto start = track.Vertex();
  auto end = track.End();
  if ((TMath::Max(start.X(), end.X()) > max_x) ||
      (TMath::Min(start.X(), end.X()) < min_x) ||
      (TMath::Max(start.Y(), end.Y()) > max_y) ||
      (TMath::Min(start.Y(), end.Y()) < min_y) ||
      (TMath::Min(start.Z(), end.Z()) < min_z) ||
      (track.Length() < min_length)) {
    return {-1., 0};
  }


  //Get the hits from the TrackHitMetas and only for view 2 (collection plane)
  protoana::ProtoDUNETrackUtils trackUtil;
  std::map<size_t, const recob::Hit *>
      hits_from_traj = trackUtil.GetRecoHitsFromTrajPoints(track, evt, trackModule);
  std::vector<const recob::Hit *> hits_from_traj_view2;
  std::vector<size_t> index_from_traj_view2;

  for (auto it = hits_from_traj.begin(); it != hits_from_traj.end(); ++it) {
    if (it->second->View() != 2) continue;
    hits_from_traj_view2.push_back(it->second);
    index_from_traj_view2.push_back(it->first);
  }

  //Find the vertex hit & info to compare to later
  double highest_z = -100.;
  int vertex_tpc = -1;
  int vertex_wire = -1;
  float vertex_peak_time = -1.;

  if (check_z) {
    for (const auto * hit : hits_from_traj_view2) {
      double wire_z = geom->Wire(hit->WireID()).GetCenter().Z();
      if (wire_z > highest_z) {
        highest_z = wire_z;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  else {
    const recob::TrackTrajectory & traj = track.Trajectory();
    double highest_diff = -1.;
    for (size_t i = 0; i < hits_from_traj_view2.size(); ++i) {
      const recob::Hit * hit = hits_from_traj_view2[i];
      size_t traj_index = index_from_traj_view2[i];
      
      double traj_x = traj.LocationAtPoint(traj_index).X();
      double traj_y = traj.LocationAtPoint(traj_index).Y();
      double traj_z = traj.LocationAtPoint(traj_index).Z();

      double diff = sqrt(std::pow((traj_x - check_x), 2) +
                         std::pow((traj_y - check_y), 2) +
                         std::pow((traj_z - check_z), 2));
      if (diff > highest_diff) {
        highest_diff = diff;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  
  std::pair<double, double> results = {0., 0.};

  //Go through all hits in the event.
  auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(hitModule);
  std::vector<art::Ptr<recob::Hit>> hit_vector;
  art::fill_ptr_vector(hit_vector, allHits);
  art::FindManyP<recob::Track> tracks_from_all_hits(allHits, evt, trackModule);
  for (size_t i = 0; i < hit_vector.size(); ++i) {

    //If this hit is in the trajectory hits vector, skip
    const recob::Hit * theHit = hit_vector[i].get();
    if (std::find(hits_from_traj_view2.begin(),
        hits_from_traj_view2.end(),
        theHit) != hits_from_traj_view2.end()) {
      continue;
    }

    //Skip hits that are outside of our TPC/plane or window of interest
    int wire = theHit->WireID().Wire;
    float peak_time = theHit->PeakTime();
    int tpc = theHit->WireID().TPC;
    int plane = theHit->View();
    if ((abs(wire - vertex_wire) > 15) ||
        (abs(peak_time - vertex_peak_time) > 100.) ||
        (tpc != vertex_tpc) || (plane != 2)) {
      continue;
    }

    //It's ok if the hits don't come from any track
    //or if that track is the primary one, because sometimes the michel hits
    //are associated to it.
    //
    //It's not ok if the hits come from another track that is long
    //(i.e. an actual daughter). Skip these
    auto & tracks_from_hit = tracks_from_all_hits.at(hit_vector[i].key());
    if (!tracks_from_hit.empty() &&
        (tracks_from_hit[0].get()->ID() != track.ID()) &&
        (tracks_from_hit[0].get()->Length() > 25.))
      continue;

    //add up the CNN results
    std::array<float, 4> cnn_out = hitResults.getOutput(hit_vector[i]);
    double hitcharge = hit_vector[i]->Integral();
    results.first += hitcharge*cnn_out[hitResults.getIndex("michel")];
    results.second += hitcharge;
  }

  return results;
}
//Ajib's implementation
std::pair<double, int> protoana::MyProtoDUNETrackUtils::GetVertexMichelScoreAlt(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    const std::string emtrkmichelidModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  // Get all tracks
  std::vector < art::Ptr < recob::Track > > trkList;
  auto trkListHandle = evt.getHandle < std::vector < recob::Track > >(trackModule);
  if (trkListHandle) {
    art::fill_ptr_vector(trkList, trkListHandle);
  }

  // Get all hits
  std::vector < art::Ptr < recob::Hit > > hitList;
  auto hitListHandle = evt.getHandle < std::vector < recob::Hit > >(hitModule);
  if (hitListHandle) {
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  // Get track-hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkListHandle, evt,trackModule); // to associate tracks and hits

  // Get hit-track association
  art::FindManyP<recob::Track> thass(hitListHandle, evt, trackModule); //to associate hit just trying

  // Get CNN scores
  anab::MVAReader<recob::Hit, 4> hitResults(evt, emtrkmichelidModule);

  int endwire = -1;
  int endtpc = -1;
  double endpeakt = -1;
  std::vector<int> wirekeys;

  if (fmthm.isValid()){
    float zlast0=-99999;
    auto vhit=fmthm.at(track.ID());
    auto vmeta=fmthm.data(track.ID());
    for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
      bool fBadhit = false;
      if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
        fBadhit = true;
        continue;
      }
      if (vmeta[ii]->Index()>=track.NumberTrajectoryPoints()){
        throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<track.NumberTrajectoryPoints()<<" for track index "<<track.ID()<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
      }
      if (!track.HasValidPoint(vmeta[ii]->Index())){
        fBadhit = true;
        continue;
      }
      auto loc = track.LocationAtPoint(vmeta[ii]->Index());
      if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
      if (loc.Z()<-100) continue; //hit not on track
      if(vhit[ii]->WireID().Plane==2){
        wirekeys.push_back(vhit[ii].key());
        float zlast=loc.Z();
        if(zlast>zlast0){
          zlast0=zlast;
          endwire=vhit[ii]->WireID().Wire;
          endpeakt=vhit[ii]->PeakTime();
          endtpc=vhit[ii]->WireID().TPC;
        }
      }
    }
  }

  int ndaughterhits = 0;
  double average_daughter_score_mic = 0;
  
  for(size_t hitl=0;hitl<hitList.size();hitl++){
    std::array<float,4> cnn_out=hitResults.getOutput(hitList[hitl]);
    auto & tracks = thass.at(hitList[hitl].key());
    // hit not on the track
    if (std::find(wirekeys.begin(), wirekeys.end(), hitl) != wirekeys.end()) continue;
    // hit not on a long track
    if (!tracks.empty() && int(tracks[0].key()) != track.ID() && trkList[tracks[0].key()]->Length()>25) continue;
    int planeid=hitList[hitl]->WireID().Plane;
    if (planeid!=2) continue;
    int tpcid=hitList[hitl]->WireID().TPC;
    if (tpcid!=endtpc) continue;
    float peakth1=hitList[hitl]->PeakTime();
    int wireh1=hitList[hitl]->WireID().Wire;
    if(std::abs(wireh1-endwire)<8 && std::abs(peakth1-endpeakt)<150 && tpcid==endtpc){
      ++ndaughterhits;
      average_daughter_score_mic += cnn_out[hitResults.getIndex("michel")];
      //std::cout<<hitList[hitl]->WireID().Wire<<" "<<hitList[hitl]->PeakTime()<<" "<<hitList[hitl]->Integral()<<" "<<cnn_out[hitResults.getIndex("michel")]<<std::endl;
    }
  }

  if (ndaughterhits) average_daughter_score_mic /= ndaughterhits;
  
  return {average_daughter_score_mic, ndaughterhits};
}
