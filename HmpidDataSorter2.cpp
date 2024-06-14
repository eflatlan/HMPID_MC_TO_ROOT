#ifndef HMPID_DATA_SORTER2_H
#define HMPID_DATA_SORTER2_H
#include "Tracks.cpp"

#include <cstdlib>  // For free
#include <cxxabi.h> // For abi::__cxa_demangle
#include <iomanip>
#include <typeinfo> // For typeid

#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Trigger.h"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <algorithm>
#include <map>
#include <vector>
// #include <omp.h>
#include "Steer/MCKinematicsReader.h"
#include <iostream>

#include <set>

struct Entry {
  int eidMip;
  int pdgHitMc;
  int tidMip;
};

class HmpidDataSorter2 {

  McTruth mcTruth;

  TFile *tFile = new TFile("data.root", "RECREATE");
  
  // Tree holding the probability of each cluster in event (sum of all species and all tracks)
  TTree *tSumProballTracks = new TTree(
      "SumProballTracks",
      "Tree probability across all tracks and species for chamber-event");

  // tree holding track-informatiom about the actual track to be reconstructed
  TTree *tTrack =
      new TTree("ThisTrack",
                "Tree containing information about track to be reconstructed");

  // tree holding track-informatiom about all the other tracks in same event and chamber
  TTree *tOtherTracks = new TTree(
      "OtherTracks",
      "Tree containing information about other tracks in chamber-event");

  // tree holding all information about all clusters in chamber and event that exceeds a certain charge-threshold
  // used for FeedBack photons
  TTree *tHighChargeClu = new TTree("HighChargeClusters",
                                    "Tree containing clusters of high charge");

  // tree holding all information about the clusters in the given event and chamber 
  TTree *tClusters =
      new TTree("ClusterCandidates", "Tree containing cluster-candidates");


  // The MC truth values
  TTree *tMcTruth = new TTree("McTruth", "Tree containing MC truth");

  // struct holding fields about all info to be stored in trees [this track = (tTrack, <field>ThisTrack), and other tracks = (tOtherTracks, <field>OtherTracks)]
  Tracks::TrackAttributes trackAttributes;

  // struct holding information about all cluster information to be stored (tClusters, )
  Tracks::ClusterData clusterData;


  std::vector<float> sumProbabilityAllTracks;

  Tracks::HighChargeClu highChargeClu;


  // count matching information
  int numPdgMipTrackOk = 0;          // distance not ok, PDG/MIP ok
  int numDist_PdgMipTrackOk = 0;    // PDG and distane ok
  int numPdgMipTrackNotOk = 0;      // PDG/MIP is ok, distance not ok 
  int numDist_PdgMipTrackNotOk = 0; // neither distane or MIP is ok

  // o2::steer::MCKinematicsReader mcReader; 
  std::unique_ptr<o2::steer::MCKinematicsReader> mcReader; // reader of MC information

  // aliases
  using MatchInfoHMP = o2::dataformats::MatchInfoHMP;
  using MCLabel = o2::MCCompLabel;
  using Cluster = o2::hmpid::Cluster;

  //"3D" structure for Clusters
  using ChamberClustersMap = std::map<int, std::vector<o2::hmpid::Cluster>>;
  using EventChamberClustersMap = std::map<int, ChamberClustersMap>;

  //"3D" structure for mathcinfo
  using ChamberMatchInfoMap =
      std::map<int, std::vector<o2::dataformats::MatchInfoHMP>>;
  using EventChamberMatchInfoMap = std::map<int, ChamberMatchInfoMap>;

  //"3D" structure for mathcinfo-MC-truth
  using ChamberMCLabelMap = std::map<int, std::vector<o2::MCCompLabel>>;
  using EventChamberMCLabelMap = std::map<int, ChamberMCLabelMap>;

  std::vector<o2::MCTrack> mcTrackArr, *mcTrackArrPtr = &mcTrackArr;

  using MCTruthContainerLabel =
      o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

public:
  // HMPIDDataSorter2(/*const std::vector<o2::hmpid::Cluster>& allClusters,
  // const std::vector<o2::dataformats::MatchInfoHMP>& allMatchInfo, const
  // std::vector<o2::MCCompLabel>& allMcLabels*/) {
  HmpidDataSorter2() {
    mcReader = std::make_unique<o2::steer::MCKinematicsReader>(
        "collisioncontext.root");
  }

  void fillTrees() {
    tSumProballTracks->Fill();
    tTrack->Fill();
    tOtherTracks->Fill();
    tHighChargeClu->Fill();
    tClusters->Fill();
    tMcTruth->Fill();
  }

  void writeTrees() {
    tSumProballTracks->Write();
    tTrack->Write();
    tOtherTracks->Write();
    tHighChargeClu->Write();
    tClusters->Write();
    tMcTruth->Write();

    std::cout << "Number of entries in tSumProballTracks: "
              << tSumProballTracks->GetEntries() << std::endl;
    std::cout << "Number of entries in tTrack: " << tTrack->GetEntries()
              << std::endl;
    std::cout << "Number of entries in tOtherTracks: "
              << tOtherTracks->GetEntries() << std::endl;
    std::cout << "Number of entries in tHighChargeClu: "
              << tHighChargeClu->GetEntries() << std::endl;
    std::cout << "Number of entries in tClusters: " << tClusters->GetEntries()
              << std::endl;
    std::cout << "Number of entries in tMcTruth: " << tMcTruth->GetEntries()
              << std::endl;
  }


  // to compare teh eventID and TID for two MCCompLabels, 
  // fx to compare a MIP to the matched track
  static bool compareMCLabels(const o2::MCCompLabel &a,
                              const o2::MCCompLabel &b) {
    if (a.getEventID() != b.getEventID()) {
      return a.getEventID() < b.getEventID();
    } else {
      return a.getTrackID() < b.getTrackID();
    }
  }

  using TrackMap =
      std::map<int, std::vector<std::pair<MCLabel, int>>>; // Inner map: trackID
                                                           // -> Labels
  using EventMap = std::map<int, TrackMap>; // Outer map: eventID -> TrackMap
  EventMap mEventMap;

  void setClusterMcTruth(
      o2::dataformats::MCTruthContainer<o2::MCCompLabel> cluLabels) {
    cluLblArr = cluLabels;

    for (int i = 0; i < cluLabels.getIndexedSize(); i++) {
      auto labels = cluLabels.getLabels(i);
      for (auto label : labels) {
        mEventMap[label.getEventID()][label.getTrackID()].push_back(
            {label, i}); // set the Cluster MC-truth position (should then
                         // corespond to cluter-index also)
      }
    }
  }

  // if no matching is found between MIP-track:
  // find the cluster-index of the eid-tid combination,
  // given by mcTrack : to find which cluster it is
  void printLabelsMcFromCluIndex(int eid, int tid,
                                 std::vector<int> &cluMipindices) {

    std::vector<std::pair<o2::MCCompLabel, int>> cluLabelsFromTrack =
        lookupLabels(eid, tid);

    LOGP(info, "cluLabelsFromTrack size {}", cluLabelsFromTrack.size());
    for (const auto &cluLabelPair : cluLabelsFromTrack) {

      const auto &cluLabel = cluLabelPair.first;
      LOGP(info,
           "        From cluLabel | Event: {}, track: {}, source: {} || "
           "Cluindex {}",
           cluLabel.getEventID(), cluLabel.getTrackID(), cluLabel.getSourceID(),
           cluLabelPair.second);

      auto cluIndex = cluLabelPair.second;

      // store all indices for cluster-MIPs (kan det i det hele tatt vare fler?)
      if (cluMipindices.size() > 1 && cluMipindices.back() != cluIndex) {
        cluMipindices.push_back(cluIndex);
      } else if (cluMipindices.size() == 0) {
        cluMipindices.push_back(cluIndex);
      }

      if (mcReader->getTrack(cluLabel)) {
        const auto &mcCluFromTrack = mcReader->getTrack(cluLabel);
        if (mcCluFromTrack) {
          int pdgCluFromTrack = mcCluFromTrack->GetPdgCode();

          LOGP(info, "        pdgCluFromTrack pdgCode = {} ", pdgCluFromTrack);
        }
      }
    }
  }

  std::vector<std::pair<MCLabel, int>> lookupLabels(int eventID, int trackID) {
    std::vector<std::pair<MCLabel, int>> pairedLabels;

    // Check if the eventID exists in the eventMap
    auto eventIt = mEventMap.find(eventID);
    if (eventIt != mEventMap.end()) {
      // Event ID found, now look for the trackID in the TrackMap
      const TrackMap &trackMap = eventIt->second;
      auto trackIt = trackMap.find(trackID);
      if (trackIt != trackMap.end()) {
        // Track ID found, return the vector of MCLabel objects
        pairedLabels = trackIt->second;
      }
    }
    // Sort pairedLabels based on the int values
    std::sort(
        pairedLabels.begin(), pairedLabels.end(),
        [](const std::pair<MCLabel, int> &a, const std::pair<MCLabel, int> &b) {
          return a.second <
                 b.second; // Compare based on the int part of the pair
        });

    // Extract MCLabel objects from the sorted pairs into the result vector
    std::vector<MCLabel> result;
    for (const auto &labelPair : pairedLabels) {
      result.push_back(
          labelPair.first); // Add only the MCLabel part of each sorted pair
    }

    if (pairedLabels.size() == 0) {
      LOGP(info, "Did not find the Cluster associated with the track!");
    }

    return pairedLabels;
  }

  void iterateOverMatchedTracks() {



    // Set branches and memberfields for the different TTRees
    tClusters->Branch("ClusterData_xValues", &clusterData.xValues);
    tClusters->Branch("ClusterData_yValues", &clusterData.yValues);
    tClusters->Branch("ClusterData_qValues", &clusterData.qValues);
    tClusters->Branch("ClusterData_sizeValues", &clusterData.sizeValues);
    tClusters->Branch("ClusterData_thetaCerValues",
                      &clusterData.thetaCerValues);
    tClusters->Branch("ClusterData_phiCerValues", &clusterData.phiCerValues);
    tClusters->Branch("ClusterData_sigmaRingValues",
                      &clusterData.sigmaRingValues);
    tClusters->Branch("ClusterData_pionProbs", &clusterData.pionProbs);
    tClusters->Branch("ClusterData_kaonProbs", &clusterData.kaonProbs);
    tClusters->Branch("ClusterData_protonProbs", &clusterData.protonProbs);
    tClusters->Branch("ClusterData_protonProbsNorm",
                      &clusterData.protonProbsNorm);
    tClusters->Branch("ClusterData_kaonProbsNorm", &clusterData.kaonProbsNorm);
    tClusters->Branch("ClusterData_pionProbsNorm", &clusterData.pionProbsNorm);
    tClusters->Branch("ClusterData_sumProbabilityTrack",
                      &clusterData.sumProbabilityTrack);
    tClusters->Branch("ClusterData_rawSizeValues", &clusterData.rawSizeValues);
    tClusters->Branch("ClusterData_numRawClustersValues",
                      &clusterData.numRawClustersValues);

    tTrack->Branch("TrackAttributes_xMipThisTrack",
                   &trackAttributes.xMipThisTrack);
    tTrack->Branch("TrackAttributes_yMipThisTrack",
                   &trackAttributes.yMipThisTrack);
    tTrack->Branch("TrackAttributes_xRadThisTrack",
                   &trackAttributes.xRadThisTrack);
    tTrack->Branch("TrackAttributes_yRadThisTrack",
                   &trackAttributes.yRadThisTrack);
    tTrack->Branch("TrackAttributes_xPCThisTrack",
                   &trackAttributes.xPCThisTrack);
    tTrack->Branch("TrackAttributes_yPCThisTrack",
                   &trackAttributes.yPCThisTrack);
    tTrack->Branch("TrackAttributes_thetaPThisTrack",
                   &trackAttributes.thetaPThisTrack);
    tTrack->Branch("TrackAttributes_phiPThisTrack",
                   &trackAttributes.phiPThisTrack);
    tTrack->Branch("TrackAttributes_momentumThisTrack",
                   &trackAttributes.momentumThisTrack);
    tTrack->Branch("TrackAttributes_qMipThisTrack",
                   &trackAttributes.qMipThisTrack);
    tTrack->Branch("TrackAttributes_sizeMipThisTrack",
                   &trackAttributes.sizeMipThisTrack);
    tTrack->Branch("TrackAttributes_mipPcDistThisTrack",
                   &trackAttributes.mipPcDistThisTrack);
    tTrack->Branch("TrackAttributes_ckovThPionThisTrack",
                   &trackAttributes.ckovThPionThisTrack);
    tTrack->Branch("TrackAttributes_ckovThKaonThisTrack",
                   &trackAttributes.ckovThKaonThisTrack);
    tTrack->Branch("TrackAttributes_ckovThProtonThisTrack",
                   &trackAttributes.ckovThProtonThisTrack);
    tTrack->Branch("TrackAttributes_refIndexThisTrack",
                   &trackAttributes.refIndexThisTrack);

    tTrack->Branch("TrackAttributes_ckovReconThisTrack",
                   &trackAttributes.ckovReconThisTrack);
    tTrack->Branch("TrackAttributes_ckovReconMassHypThisTrack",
                   &trackAttributes.ckovReconMassHypThisTrack);

    tTrack->Branch("TrackAttributes_numCkovHough",
                   &trackAttributes.numCkovHough);
    tTrack->Branch("TrackAttributes_numCkovHoughMH",
                   &trackAttributes.numCkovHoughMH);



    tOtherTracks->Branch("TrackAttributes_xMipsOtherTracks",
                         &trackAttributes.xMipsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_yMipsOtherTracks",
                         &trackAttributes.yMipsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_xRadsOtherTracks",
                         &trackAttributes.xRadsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_yRadsOtherTracks",
                         &trackAttributes.yRadsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_xPCsOtherTracks",
                         &trackAttributes.xPCsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_yPCsOtherTracks",
                         &trackAttributes.yPCsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_thetaPsOtherTracks",
                         &trackAttributes.thetaPsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_phiPsOtherTracks",
                         &trackAttributes.phiPsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_momentumsOtherTracks",
                         &trackAttributes.momentumsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_qMipsOtherTracks",
                         &trackAttributes.qMipsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_sizeMipsOtherTracks",
                         &trackAttributes.sizeMipsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_mipPcDistsOtherTracks",
                         &trackAttributes.mipPcDistsOtherTracks);
    tOtherTracks->Branch("TrackAttributes_ckovThPionOtherTracks",
                         &trackAttributes.ckovThPionOtherTracks);
    tOtherTracks->Branch("TrackAttributes_ckovThKaonOtherTracks",
                         &trackAttributes.ckovThKaonOtherTracks);
    tOtherTracks->Branch("TrackAttributes_ckovThProtonOtherTracks",
                         &trackAttributes.ckovThProtonOtherTracks);
    tOtherTracks->Branch("TrackAttributes_refIndexesOtherTracks",
                         &trackAttributes.refIndexesOtherTracks);

    tOtherTracks->Branch("TrackAttributes_ckovReconOtherTracks",
                         &trackAttributes.ckovReconOtherTracks);
    tOtherTracks->Branch("TrackAttributes_ckovReconMassHypOtherTracks",
                         &trackAttributes.ckovReconMassHypOtherTracks);


    tHighChargeClu->Branch("highChargeClu_x", &highChargeClu.x);
    tHighChargeClu->Branch("highChargeClu_y", &highChargeClu.y);
    tHighChargeClu->Branch("highChargeClu_q", &highChargeClu.q);
    tHighChargeClu->Branch("highChargeClu_size", &highChargeClu.size);


    tSumProballTracks->Branch("sumProbabilityAllTracks",
                              &sumProbabilityAllTracks);


    tMcTruth->Branch(
        "mcTruth_isTrackToReconKnownPdg",
        &mcTruth.isTrackToReconKnownPdg); // > if we take into account known pdg
                                          // for requirements
    tMcTruth->Branch("mcTruth_isMipMatchedCorrectly",
                     &mcTruth.isMipMatchedCorrectly);
    tMcTruth->Branch("mcTruth_pdgCodeTrack;", &mcTruth.pdgCodeTrack);
    tMcTruth->Branch("mcTruth_pdgCodeClu", &mcTruth.pdgCodeClu);
    tMcTruth->Branch("mcTruth_numCountedPhotonsPC",
                     &mcTruth.numCountedPhotonsPC);

    // MC-truth
    // int pdgCodeTrack, pdgCodeClu;
    // bool isMipMatchedCorrectly;

    LOGP(info, "\n=======================================");
    LOGP(info, "iterateOverMatchedTracks");

    // to get sigmaRing
    const double refIndexTmp = 1.2904;
    int trigNum = 0;
    for (const auto &eventEntry : matchInfoByEventChamber) {
      LOGP(info, "\n=======================================");
      LOGP(info, "  Event {}", eventEntry.first);
      LOGP(info, "=======================================");


      const auto &trig = triggers[eventEntry.first]; 

      std::vector<o2::hmpid::Cluster> clustersInEvent;
      // std::vector<o2::hmpid::Cluster> clusterLabelsInEvent;

      // std::ordered_map<int, std::vector<o2::DataFormatsHMP::cluster>>
      // clusterMaps;
      std::array<std::vector<o2::hmpid::Cluster>, 7> clusterArray;

      // ef : do more elegant, but aon this
      // to map mipIndex which has just the index of the mip
      // in one event and one chamber of clusters
      std::array<std::vector<int>, 7> cluIndexArray;

      int tnum = 0;
      int cluNumPre = 0;

      for (int cluNum = trig.getFirstEntry(); cluNum <= trig.getLastEntry();
           cluNum++) {

        if (trig.getNumberOfObjects() < 1) {
          continue;
        }

        if (cluNum < mClusters.size()) {
          const auto &clu = mClusters[cluNum];

          const int chNum = clu.ch();
          if (chNum >= o2::hmpid::Param::EChamberData::kMinCh &&
              chNum <= o2::hmpid::Param::EChamberData::kMaxCh) {

            clustersInEvent.emplace_back(clu);
            clusterArray[chNum].emplace_back(clu);

            // sok som
            // indexOfMipGlobal = cluIndexArray[chNum][index];
            cluIndexArray[chNum].emplace_back(cluNum);
          }
          // cluNumPre = evNum;
        }
        tnum++;
      }


      // Looping over all tracks for the given event
      for (const auto &chamberEntry : eventEntry.second) {

        // here all the matched and unmatched tracks for a given event and
        // chamber

        // get clusters for the given event and chamber
        // const auto& clusters =
        // clustersByEventChamber[eventEntry.first][chamberEntry.first]; const
        // auto& clusters = getClusters(eventEntry.first, chamberEntry.first);
        const auto &clustersInChamber = clusterArray[chamberEntry.first];

        const auto &mcMatchInfoArr =
            mcMatchInfoByEventChamber[eventEntry.first][chamberEntry.first];



        LOGP(info,
             "clustersInEvent size {} : in chamber-- clustersInChamber size {}",
             clustersInEvent.size(), clustersInChamber.size());

        int trackNum = 0;
        LOGP(info, "=======================================");
        LOGP(info, "    Event {}, chamber {} : num Clusters = {}",
             eventEntry.first, chamberEntry.first, clustersInChamber.size());

        const o2::MCTrack *mcTrack = nullptr;
        std::vector<MCRecon> mcReconObjects;


        // looping over all the chambers
        // chamberEntry.second is the array of MatchInfo for the given chamber number 
        // chamber-number given by chamberEntry.first
        for (const auto &matchInfo : chamberEntry.second) {
          int hardCodedMipIndex = matchInfo.getMipclusIndex();

          float xMip, yMip;
          int qMip, nph;

          float xPcUnc, yPcUnc, xRad, yRad, thetaP, phiP;

          matchInfo.getHMPIDtrk(xPcUnc, yPcUnc, thetaP, phiP);

          matchInfo.getHMPIDmip(xMip, yMip, qMip, nph);


          // ef > TODO > this field should not exist in MatchInfoHMP > replace it with the linear extrapolation
          matchInfo.getHMPIDrad(xRad, yRad);

          const auto dist = TMath::Sqrt((xPcUnc - xMip) * (xPcUnc - xMip) +
                                        (yPcUnc - yMip) * (yPcUnc - yMip));


          if (dist > 6.) {
            Printf("               Too large distance %.0f : Skip", dist);
            trackNum++;
            continue;
          }
          /*
          if (!matchInfo.getMatchStatus()) {
            Printf("               was not matched");
            continue;
          }*/

          LOGP(info,
               "        "
               "\n\n=========================================================");
          LOGP(info, "        Track number {} : matched Status {}", trackNum,
               matchInfo.getMatchStatus());


          // ef > verbose pringing
          /*
          Printf("              xMip %.1f xPcUnc %.1f ", xMip, xPcUnc);
          Printf("               yMip %.1f yPcUnc %.1f ", yMip, yPcUnc);

          // LOGP(info, "        Dist PC CONC-UNC: deltaX {} deltaY {}", xPcCon
          // - xPcUnc, yPcCon - yPcUnc);
          Printf("               PC_CONC-MIP : deltaX %.1f deltaY %.1f",
                 xPcUnc - xMip, yPcUnc - yMip);
          Printf("               dist %.1f", dist);
          Printf("               xPc %.1f, mip %.1f, ra %.1f | yPc %.1f mip "
                 "%.1f ra %.1f",
                 yPcUnc, xMip, xRad, yPcUnc, yMip, yRad);

          */

          // get MC truth
          int trackIdKine = -1;
          int eventIdKine = -1;
          int sourceIdKine = -1;

          const auto &mcMatchInfo = mcMatchInfoArr[trackNum];

          if (!mcMatchInfo.isValid()) {
            LOGP(info, "        mcMatchInfo not valid");
          }

          trackIdKine = mcMatchInfo.getTrackID();
          eventIdKine = mcMatchInfo.getEventID();
          sourceIdKine = mcMatchInfo.getSourceID();

          // int eventID, int trackID

          // ef : get cluMC truths for the given track


          // ef : if we want to find the MIP of the clusters indicated by the
          // track
          bool printCluLabels = false;
          if (printCluLabels) {
            std::vector<std::pair<o2::MCCompLabel, int>> cluLabelsFromTrack =
                lookupLabels(eventIdKine, trackIdKine);

            LOGP(info, "cluLabelsFromTrack size {}", cluLabelsFromTrack.size());
            for (const auto &cluLabelPair : cluLabelsFromTrack) {

              const auto &cluLabel = cluLabelPair.first;
              LOGP(info,
                   "        From cluLabel | Event: {}, track: {}, source: {} "
                   "|| Cluindex {}",
                   cluLabel.getEventID(), cluLabel.getTrackID(),
                   cluLabel.getSourceID(), cluLabelPair.second);

              if (mcReader->getTrack(cluLabel)) {
                const auto &mcCluFromTrack = mcReader->getTrack(cluLabel);

                if (mcCluFromTrack) {
                  int pdgCluFromTrack = mcCluFromTrack->GetPdgCode();

                  LOGP(info, "        pdgCluFromTrack pdgCode = {} ",
                       pdgCluFromTrack);
                }
              }
            }
          }

          LOGP(info,
               "        From mcMatchInfo | Event: {}, track: {}, source: {}",
               eventIdKine, trackIdKine, sourceIdKine);


          if (mcReader == nullptr) {
            LOGP(info, "      mcReader == nullptr");
          }

          try {

            /*LOGP(info, "        try mcReader w IP types: trackIdKine: {} {},
             * eventIdKine: {} {}, sourceIdKine: {} {}", trackIdKine,
             * typeid(trackIdKine).name(), eventIdKine,
             * typeid(eventIdKine).name(), sourceIdKine,
             * typeid(sourceIdKine).name());*/

            mcTrack = mcReader->getTrack(mcMatchInfo);
            // mcTrack = mcReader->getTrack(sourceIdKine, eventIdKine,
            // trackIdKine);
          } catch (const std::exception &e) {
            LOGP(error,
                 "       Exception caught while trying to read MC track: %s",
                 e.what());
          } catch (...) {
            LOGP(error, "       Unknown exception caught while trying to read "
                        "MC track");
          }

          /*
          LOGP(info, "        From mcTrack | Event: {}, track: {}, source: {}",
          eventIdKine, trackIdKine, sourceIdKine);*/

          if (mcTrack == nullptr) {
            LOGP(info, "        MC track not found for event {} and track {}",
                 eventIdKine, trackIdKine);
            continue;
          }

          TParticlePDG *pPDG =
              TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode());

          auto pdgCode = mcTrack->GetPdgCode();



          // index of first and last daughter 
          const int startIdxDaughter = mcTrack->getFirstDaughterTrackId();
          const int lastIdxDaughter = mcTrack->getLastDaughterTrackId();

          LOGP(info, "firstDaugter {} lastDaughtr {}", startIdxDaughter,
               lastIdxDaughter);

          const int numCkovPhotsTotal = lastIdxDaughter - startIdxDaughter;



          // > ef : i dont thinjk this is correct approach,
          // should probably instead loop over the indices from startIdxDaughter > lastIdxDaughter and count all that has the motherTrack as trackID
          LOGP(info, "Number of Ckov Photons for track {}", numCkovPhotsTotal);

          auto pT = TMath::Abs(mcTrack->GetStartVertexMomentumX() *
                                   mcTrack->GetStartVertexMomentumX() +
                               mcTrack->GetStartVertexMomentumY() *
                                   mcTrack->GetStartVertexMomentumY());

          int pdgFromTrack = matchInfo.getParticlePdg();
          int mipPDG = -1;

          LOGP(info, "        Particle {}: mcTrack pdgCode = {}", trackNum,
               pdgCode);


          // flag stating if MIP was correctly matched with the corresponding track
          // can only be done in MC    
          bool isMipMatched = false;

          if (true /*pdgCode != pdgFromMip*/) {
            // LOGP(info, "      PDG code ulik! : num clustersInChamber = {}",
            // clustersInChamber.size());

            LOGP(info, "      **********************************");
            LOGP(info, "      Check clusters \n");

            const int mipcluSize = matchInfo.getMipClusSize();

            // void setIdxHMPClus(int ch, int idx) { mIdxHMPClus = ch * 1000000
            // + idx; } matching.setIdxHMPClus(iCh, index + 1000 * cluSize); //
            // set chamber, index of cluster + cluster size

            const int mipIndex = matchInfo.getIdxHMPClus();

            int chTrack = matchInfo.getChamber();

            int mipch = mipIndex / 1000000; // Extract chamber
            int remainder =
                mipIndex %
                1000000; // Remainder after removing chamber information
            int mipSz = remainder / 1000; // Extract cluster size
            int index = remainder % 1000; // Extract index of cluster

            const std::vector<int> &clustersIndexChamber = cluIndexArray[mipch];

            int mipIndexGobal = 0;
            if (index < clustersIndexChamber.size()) {
              mipIndexGobal = clustersIndexChamber[index];
            } else {
              LOGP(
                  info,
                  "      index {} greater than size of clustersIndexChamber {}",
                  index, clustersIndexChamber.size());
            }
            LOGP(info, "      mipIndexGobal {}", mipIndexGobal);
            // indexOfMipGlobal = cluIndexArray[chNum][index];

            const int cluTrigStartIndex = trig.getFirstEntry();
            const int cluTrigEndIndex = trig.getLastEntry();


            // total number of clusters aon
            const int numCluTotal = mClusters.size();

            // number of clusters in current HMP trigger
            const int numCluInTrig = trig.getNumberOfObjects();



            /* ef > verbose printing to debug MIP index 
            LOGP(info,
                 "      clustersInEvent size {} simpleClusters size {} "
                 "numCluInTrig {} ",
                 clustersInEvent.size(), simpleClusters.size(), numCluInTrig);

            LOGP(info,
                 "       cluTrigStartIndex {} cluTrigEndIndex {} numCluTotal "
                 "{} : numCluInTrig {} indexMIP {}   hardCodedMipIndex {}",
                 cluTrigStartIndex, cluTrigEndIndex, numCluTotal, numCluInTrig,
                 mipIndexGobal, hardCodedMipIndex);
            */

            // to check index :
            const int indexTotal = cluTrigStartIndex + index;



            Printf("              mipIndex %d mipch %d mipSz %d mipIndexGobal "
                   "(%d/%d)",
                   mipIndex, mipch, mipSz, mipIndexGobal, numCluTotal);

            // get hit-->dig-->clu MC-truth for MIP
            const auto &mipLabels =
                cluLblArr.getLabels(hardCodedMipIndex /*mipIndexGobal*/);
            LOGP(info, "      Get MC-truth for MIP : size {}",
                 mipLabels.size());

            int indexLabel = 1;

            // ef : sort labels on 1. eventID 2. trackID
            std::sort(mipLabels.begin(), mipLabels.end(), compareMCLabels);


            // loop over all hit MC labels for the MIP
            for (const auto &mipLabel : mipLabels) {
              if (isMipMatched) {
                break;
              }

              const auto &mcTruthHit = mcReader->getTrack(mipLabel);


              // receive MC infomration about MIP
              int pdgHitMc = mcTruthHit->GetPdgCode();
              const auto &tidMip = mipLabel.getTrackID();
              const auto &eidMip = mipLabel.getEventID();

              // check pdg and trackID of matched
              // MIP-track pair
              if (pdgCode == pdgHitMc && tidMip == trackIdKine) {
                // LOGP(info, "MC track-cluster matched \n pdgCode "


                // pdg code and tid for Track and matched MIP was the same
                isMipMatched = true;
                mipPDG = pdgHitMc;
                break;
              }

              if (pdgHitMc != 50000050) {
                mipPDG = pdgHitMc;
              }


              /* ef : verbose printing
              if (!isMipMatched) {
                Printf("              MC label %d: eventID = %d, trackID = %d, "
                       "sourceID = %d",
                       indexLabel, mipLabel.getEventID(), mipLabel.getTrackID(),
                       mipLabel.getSourceID());
                Printf("            Hit MC-truth %d/%zu : pdg %d", indexLabel,
                       mipLabels.size(), pdgHitMc);
              } */

              indexLabel++;
            }

            // ef :compare track-label
            // with matched mipch
            // fx check pdg code,
            // can also chekc tradkID etc
            if (isMipMatched) {
              LOGP(info, "      matched with correct cluster\n\n\n");

              if (matchInfo.getMatchStatus()) {
                numDist_PdgMipTrackOk++;
              } else {
                numPdgMipTrackOk++;
              }

            } 
            
            
            // track was not matched with correct MIP
            // meaning that another MIP-candidate was closer to the extrapolated track
            else {
              if (matchInfo.getMatchStatus()) {
                numDist_PdgMipTrackNotOk++;
              } else {
                numPdgMipTrackNotOk++;
              }

              // ef : return index of cluster, and print the pos etc

              LOGP(info, "      matched with wrong cluster\n\n\n");
              std::vector<int> cluMipindices;

              LOGP(info, "      Looking up cluster-MIP index");

              printLabelsMcFromCluIndex(eventIdKine, trackIdKine,
                                        cluMipindices);

              // ef : TODO : read cluster-MC truth index of mip
              // retrieve cluster MIP info and positons
              LOGP(info, "      Looking up cluster-MIP positions");

              for (const int &cluMipIndex : cluMipindices) {

                if (cluMipIndex < mClusters.size()) {

                  // MIP matched w track
                  const auto &clu =
                      mClusters[hardCodedMipIndex /*mipIndexGobal*/];

                  // actual MIP per indexing
                  const auto &clu2 = mClusters[cluMipIndex];

                  const auto distMipMatched2PC =
                      TMath::Sqrt((clu.x() - xPcUnc) * (clu.x() - xPcUnc) +
                                  (clu.y() - yPcUnc) * (clu.y() - yPcUnc));



                  // several print statements for comparing actual MIP to matched MIP
                  {

                    auto labels = cluLblArr.getLabels(
                        hardCodedMipIndex /*mipIndexGobal*/);
                    LOGP(info, "        checking MC-labels  size {}",
                         labels.size());

                    for (auto label : labels) {
                      LOGP(info, "        evID {} tid {}", label.getEventID(),
                           label.getTrackID());
                    }
                  }

                  {

                    auto labels = cluLblArr.getLabels(cluMipIndex);
                    LOGP(info, "        checking MC-labels  size {}",
                         labels.size());

                    for (auto label : labels) {
                      LOGP(info, "        evID {} tid {}", label.getEventID(),
                           label.getTrackID());
                    }
                  }

                  const auto distMipCorrect2PC =
                      TMath::Sqrt((clu2.x() - xPcUnc) * (clu2.x() - xPcUnc) +
                                  (clu2.y() - yPcUnc) * (clu2.y() - yPcUnc));

                  Printf("          chamber for matched MIP %d : chamber "
                         "correct MIP %d",
                         clu.ch(), clu2.ch());

                  LOGP(info,
                       "      possibly misplaced cluster ? hardCodedMipIndex "
                       "{}  mipIndexGobal {} actualIndexMip  {}",
                       hardCodedMipIndex, mipIndexGobal, cluMipIndex);

                  /*
                  {
                  for(int iLabel = 0; iLabel<cluLblArr.getIndexedSize();
                  iLabel++)
                  {
                  LOGP(info, "EVENT-info for matched cluster {}",
                  clu.getEventNumber()); auto labels =
                  cluLblArr.getLabels(mipIndexGobal);

                  LOGP(info, "checking MC-labels  size {}", labels.size());

                  for(auto label : labels ) {
                  LOGP(info, "evID {} tid {}", label.getEventID(),
                  label.getTrackID());
                  }
                  }
                  }*/

                  const int trigFirst = trig.getFirstEntry();
                  const int trigLast = trig.getLastEntry();



                  // ef: TODO: find if mipIndexGobal and cluMipIndex are in
                  // different triggers then sth wrong with logic
                  // if(mipIndexGobal)

                  // ef : verbose printing
                  if (false) {
                    LOGP(info, "      trigFirst {} trigLast  {}", trigFirst,
                       trigLast);
                    Printf("          Actual Cluster at x %.1f y %.1f : mip "
                          "Matched was xMip %.1f yMip %.1f xPc %.1f yPc %.1f",
                          clu2.x(), clu2.y(), xMip, yMip, xPcUnc, yPcUnc);

                    Printf("          mip mipIndexGobal %d :: at x %.1f y %.1f",
                          mipIndexGobal, clu.x(), clu.y());

                    Printf("          MIP-PC for matched MIP %.1f : for correct "
                          "MIP %.1f",
                          distMipMatched2PC, distMipCorrect2PC);

                    Printf("          matched MIP (q = %.1f, size %d) : correct "
                          "MIP (q = %.1f, size %d)",
                          clu.q(), clu.size(), clu2.q(), clu2.size());
                  }




                  Printf("\n\n\n");

                } else {
                  LOGP(info, "        out of range");
                }
              }

              LOGP(info, "_________");

              // did not find the cluster, possible that the cluster from the
              // track was in DEAD-time?



              // ef : verbose checking of other MIP candidates in proximity
              bool lookTroughMips = false;
              if (cluMipindices.size() == 0 && lookTroughMips) {
                LOGP(info, "        Did not find cluster for the track, even "
                           "in other triggers ");
                LOGP(info,
                     "        Look through other clusters in proximity  ");

                for (int cluNum = trig.getFirstEntry();
                     cluNum <= trig.getLastEntry(); cluNum++) {
                  const auto &c = mClusters[cluNum];

                  if (c.ch() != mipch)
                    continue;

                  if (c.q() > 100 && c.size() > 3 && c.size() < 10) {
                    const auto cluDist2PC =
                        TMath::Sqrt((c.x() - xPcUnc) * (c.x() - xPcUnc) +
                                    (c.y() - yPcUnc) * (c.y() - yPcUnc));

                    if (cluDist2PC < 10) {
                      const auto &clabels = cluLblArr.getLabels(cluNum);

                      std::set<std::pair<int, int>> uniqueCombinations;
                      std::vector<Entry> entries;

                      for (const auto &clabel : clabels) {
                        const auto &cluTruthHit = mcReader->getTrack(clabel);

                        int pdgHitMc = cluTruthHit->GetPdgCode();

                        const auto &tidMip = clabel.getTrackID();

                        const auto &eidMip = clabel.getEventID();

                        std::pair<int, int> pdgTidPair =
                            std::make_pair(pdgHitMc, tidMip);

                        if (uniqueCombinations.insert(pdgTidPair).second) {
                          entries.push_back({eidMip, pdgHitMc, tidMip});
                        }
                      }
                      for (const auto &entry : entries) {
                        std::cout << "              eid: " << entry.eidMip
                                  << ", pdg: " << entry.pdgHitMc
                                  << ", tid: " << entry.tidMip << std::endl;
                      }
                    } // if(cluDist2PC < 10) {
                  }   //  if(c.q() < 150 || c.size() < 3 )
                }     // clus in trig
              }       // cluMipindices-== 0
              // did not find the cluster, possible that the cluster from the
              // track was in DEAD-time?
            } // else MIP was not matched
          }


          LOGP(info, "isMipMatched {} pdgCodeTrack {} mipPDG {}", isMipMatched,
               pdgCode, mipPDG);

          // MC-truth
          // int pdgCodeTrack, pdgCodeClu;
          // bool isMipMatchedCorrectly;

          // matchInfo.getHMPIDtrk(xRad, yRad, xPcUnc, yPcUnc, thetaP, phiP);
          // matchInfo.getHMPIDmip(xMip, yMip, qMip, nph);

          // store TID of photons that are given from the track
          std::vector<int> tidOfPhotsFromTrk;

          // clusters that contain Ckov photon hits
          std::vector<o2::hmpid::Cluster> ckovPhotonClustersFromTrk;
          // ef > change toe clustersInChamber instead of simpleCluster
          std::vector<o2::hmpid::Cluster> cluTemps;

          std::vector<int> indexOfCkovClusters;

          int indexOfCluInCh;
          for (int cluNum = trig.getFirstEntry(); cluNum <= trig.getLastEntry();
               cluNum++) {

            if (cluNum < mClusters.size()) {
              const auto &clu = mClusters[cluNum];

              const int chNum = clu.ch();
              if (!(chNum >= o2::hmpid::Param::EChamberData::kMinCh &&
                    chNum <= o2::hmpid::Param::EChamberData::kMaxCh)) {
                continue;
              }

              // if(clu.q() > 200 && clu.size() >= 3) {
              //   continue;
              // }

              auto x = clu.x();
              auto y = clu.y();

              auto dist2mip =
                  std::sqrt((xMip - x) * (xMip - x) + (yMip - y) * (yMip - y));

              const auto &cluLabels = cluLblArr.getLabels(cluNum);

              for (const auto &cluLbl : cluLabels) {

                const auto &eid = cluLbl.getEventID();
                const auto &tid = cluLbl.getTrackID();
                const auto &sid = cluLbl.getSourceID();

                const auto &mcClu = mcReader->getTrack(cluLbl);

                const auto &motherTid = mcClu->getMotherTrackId();
                const auto &motherTid2 = mcClu->getSecondMotherTrackId();


                // ef > quite verbose printing of Ckov photons that belong to track
                if (motherTid == trackIdKine) {
                  std::cout << "This is a ckov photon of the track!\n";
                  Printf("Cluster > dist2mip %.2f", dist2mip);
                  Printf("charge %.2f size %d", clu.q(), clu.size());
                  LOGP(info,
                       "cluster-size {} localMaximums {} number of labels {}",
                       clu.size(), clu.numLocMax(), cluLabels.size());

                  LOGP(info,
                       " CluLabel | Event: {}, track: {}, source: {} >> PDG {}",
                       eid, tid, sid, mcClu->GetPdgCode());

                  // add the index of the current cluster if it is was not added
                  // already

                  if (indexOfCkovClusters.size() == 0) {
                    indexOfCkovClusters.push_back(indexOfCluInCh);
                    LOGP(info, "added ckov clu for trk size now {}",
                         indexOfCkovClusters.size());
                  } else {
                    if (indexOfCluInCh != indexOfCkovClusters.back()) {

                      LOGP(info, "indexOfCluInCh {} back {}", indexOfCluInCh,
                           indexOfCkovClusters.back());

                      indexOfCkovClusters.push_back(indexOfCluInCh);
                      LOGP(info, "added ckov clu for trk size now {}",
                           indexOfCkovClusters.size());
                    }
                  }

                  bool isAdded = false;
                  for (const auto &tidPhot : tidOfPhotsFromTrk) {
                    if (tid == tidPhot) {
                      isAdded = true;
                    }
                  }
                  if (!isAdded) {
                    tidOfPhotsFromTrk.push_back(tid);
                    ckovPhotonClustersFromTrk.push_back(clu);
                  }

                  // trackIdKine = mcMatchInfo.getTrackID();
                  // eventIdKine = mcMatchInfo.getEventID();
                  // sourceIdKine = mcMatchInfo.getSourceID();

                  // getTrack(int source, int event, int track)
                  // LOGP(info, " CluLabel | MotherTID {}, 2ndMotherTID: {} ",
                  // motherTid, motherTid2);

                  // this means the cluster-hit label should correspond to a
                  // Ckov Photon becasue the motherID of the cluster-label
                  // matches the TID from teh MC-track from ITS

                  const auto &mcCluMother2 =
                      mcReader->getTrack(sid, eid, motherTid2);
                  const auto &mcCluMother =
                      mcReader->getTrack(sid, eid, motherTid);
                  if (mcCluMother) {
                    std::cout << "  Mother > TID : " << motherTid
                              << " PDG :" << mcCluMother->GetPdgCode();
                  }

                  if (mcCluMother2 && motherTid2 != -1) {
                    std::cout << "  || 2ndMother > TID : " << motherTid2
                              << " PDG :" << mcCluMother2->GetPdgCode();
                  }
                  std::cout << std::endl << std::endl;
                }
              }

              cluTemps.push_back(clu);
              indexOfCluInCh++;
            }
          }

          std::vector<o2::hmpid::Cluster> diffCkovclusters;
          float xPrev = 0, yPrev = 0;
          // LOGP(info, "Looping over Ckov photons");
          for (const auto &ckovPhotClu : ckovPhotonClustersFromTrk) {
            auto x = ckovPhotClu.x();
            auto y = ckovPhotClu.y();
            auto dist2mip =
                std::sqrt((xMip - x) * (xMip - x) + (yMip - y) * (yMip - y));
            // Printf("Cluster > dist2mip %.2f",  dist2mip);
            // Printf("charge %.2f size %d",  ckovPhotClu.q(),
            // ckovPhotClu.size());

            if (x != xPrev && y != yPrev) {
              diffCkovclusters.push_back(ckovPhotClu);
            }
            xPrev = x;
            yPrev = y;
          }

          // Total number of Ckov photons from track
          LOGP(info, "total number of Ckov Photons from track {}",
               numCkovPhotsTotal);

          // Only those on the HMPID PC
          LOGP(info, "Number of Clusters/Photons on HMPID PC : ");
          LOGP(info, "        number of Ckov Photons from track {}",
               tidOfPhotsFromTrk.size());
          LOGP(info,
               "        number of Clusters containing ckov photons from track "
               "> {}",
               ckovPhotonClustersFromTrk.size());
          LOGP(info,
               "        number of differemt Clusters containing ckov photons "
               "from track > {}",
               diffCkovclusters.size());

          MCRecon mcReconObj(xRad, yRad, clustersInChamber, matchInfo);

          // Recon reconObj2(xRad, yRad, clustersInChamber, matchInfo);

          LOGP(info, " Recon reconOb xRad {}, yRad {}", xRad, yRad);

          mcReconObj.process(pdgCode, indexOfCkovClusters);

          McTruth mcTruthCp(isMipMatched, pdgCode, mipPDG,
                            mcReconObj.isTrackAppropriate(),
                            diffCkovclusters.size());

          // pass object by const ref,
          mcReconObj.setMcTruth(mcTruthCp);

          // MC-truth
          // int pdgCodeTrack, pdgCodeClu;
          // bool isMipMatchedCorrectly;
          mcReconObjects.push_back(mcReconObj);
          trackNum++;
        }

        bool tracksToReconstruct = false; 

        // check if there is any valid tracks according to requirements
        // can be dione in loop above also
        for (auto &mcReconObj : mcReconObjects) {
          // auto sumProbAllTracks +=
          if (mcReconObj.isTrackAppropriate()) { // isTrackAppropriate() returns bool  mIsTrackToReconstruct
            tracksToReconstruct = true;
          }
        }

        // if no valid tracks,
        if (!tracksToReconstruct) {
          continue;
        }

        // ef > maybe set it more elegatly based on q-dist to see how many %
        // typically is above 50 ADC
        auto numHighCharge = std::ceil(clustersInChamber.size() / 2);
        if (numHighCharge < 1) {
          numHighCharge = 1;
        }

        highChargeClu.setNumberOfEntries(numHighCharge);



        // adding all clusters above a certain cut to array of clusters with high charge
        // to help with feedback photons
        for (const auto &clu : clustersInChamber) {
          if (clu.q() > 50) {
            highChargeClu.addToFields(clu.x(), clu.y(), clu.q(), clu.size());
          }
        }

        // clear
        sumProbabilityAllTracks.clear();

        // resize for the cluster-properties
        // also sets all fields to zero
        clusterData.setNumberOfEntries(clustersInChamber.size());

        sumProbabilityAllTracks.resize(clustersInChamber.size(), 0.0f);

        // number of matched tracks in trigger--chamber
        const size_t numMatchedTracks = mcReconObjects.size();

        // resize track-attributes, and
        trackAttributes.setNumberOfTracks(numMatchedTracks);

        // find all matched tracks in chamber for this trigger
        if (mcReconObjects.size() > 1) {
          int reconNum = 0;
          Tracks::Tracks tracks; // fix namespace etc..
          for (auto &mcReconObj : mcReconObjects) {
            //if (mcReconObj.isTrackMatched() {

              // fill all track attributes from all tracks that is matched
              tracks.addTrackCands(&mcReconObj);

              LOGP(info, "exit addTracks");

              // calculates the normalized probability across all species within
              // the track passes sumProbabilityAllTracks by referene, and adds
              // the probability of the entire track (sumProbabilityAllTracks +=
              // sumProbTrack)
              mcReconObj.calculateNormalizedProb(sumProbabilityAllTracks);

              std::cout << "\n\n sumProbabilityAllTracks tnum" << reconNum++
                        << "/" << mcReconObjects.size() << std::endl;
              for (size_t i = 0; i < sumProbabilityAllTracks.size(); ++i) {
                if (sumProbabilityAllTracks[i] != 0.0) {
                  std::cout << "i: " << i << " " << std::fixed
                            << std::setprecision(2)
                            << sumProbabilityAllTracks[i] << " ";
                }
              }
              std::cout << "\n" << std::endl;
            // }
          }

          // find all tracks in chamber for this trigger that meets further
          // requirments
          int indexOfReconObj = 0;
          for (auto &mcReconObj : mcReconObjects) {
            // auto sumProbAllTracks +=
            if (mcReconObj.isTrackAppropriate() or true) {

              // mcReconObj.setSumProbabilityAcrossTracks(sumProb); // prob across
              // all species and tracks

              // mcReconObj.setOtherTrackAttributes(); // prob across all species
              // and tracks xValues = mcReconObj.getXvalues();

              // pass mcTruth by ref to get the struct that was sat in mcReconObj
              mcReconObj.getMcTruth(mcTruth);

              // fill all track attributes (for this and other tracks)
              writeTrackAttributes(tracks, indexOfReconObj);
              writeClusterAttributes(mcReconObj, clusterData);

              fillTrees();

              trackAttributes.clear();

              LOGP(info, "filled TTRee entry");

              // clear all vectors
              // clearFields();
            }
            indexOfReconObj++;
          }
        }

        // just one object, we dont need to normalize
        else if (mcReconObjects.size() == 1) {
          auto &mcReconObj = mcReconObjects[0];

          Tracks::Tracks tracks; // fix namespace etc..

          if (mcReconObj.isTrackAppropriate() or true) {

            LOGP(info, "only 1 reconObject >> isTrackAppropriate");

            mcReconObj.calculateNormalizedProb(sumProbabilityAllTracks);

            tracks.addTrackCands(&mcReconObj);
            writeTrackAttributes(tracks, 0);
            writeClusterAttributes(mcReconObj, clusterData);

            // pass mcTruth by ref to get the struct that was sat in mcReconObj
            mcReconObj.getMcTruth(mcTruth);

            fillTrees();
            LOGP(info, "filled TTRee entry");

            trackAttributes.clear(); // ef added

            // clear all vectors
            // clearFields();
          }
        }

        highChargeClu.clearDataFields();

        // loop over reconObject, do postprocessing ...
      }

      LOGP(info, "exit eventLoop");
    }

    LOGP(info, "number of distance ok  and match ok {}", numDist_PdgMipTrackOk);
    LOGP(info, "number of distance X and match ok {}", numPdgMipTrackOk);
    LOGP(info, "number of distance ok and match X {}",
         numDist_PdgMipTrackNotOk);
    LOGP(info, "number of distance X and match X {}", numPdgMipTrackNotOk);

    tFile->cd();

    writeTrees();

    LOGP(info, "exit iterateMC");
    tFile->Close();
    LOGP(info, "exit iterateMC");
  }

  void
  organizeAndSortClusters(const std::vector<o2::hmpid::Cluster> &allClusters) {

    LOGP(info, "organizeAndSortClusters : numClusters IN  {}; OUT = {} ",
         allClusters.size(), mClusters.size());

    for (const auto &cluster : allClusters) {
      mClusters.emplace_back(cluster);


      // clustersByEventChamber[eventNum][chamberNum].push_back(cluster);
      // LOGP(info, "organizeAndSortClusters : eventNum {} chamberNum {} : size
      // {}", eventNum, chamberNum,
      // clustersByEventChamber[eventNum][chamberNum].size()); LOGP(info,
      // "organizeAndSortClusters : eventNum {} chamberNum {} : size {}",
      // eventNum, chamberNum, getClusters(eventNum,chamberNum).size());
    }

    int tnum = 0;
    LOGP(info, "organizeAndSortClusters : nTriggers {}", triggers.size());
  }

  void organizeAndSortMatchInfo(
      const std::vector<o2::dataformats::MatchInfoHMP> &allMatchInfo,
      const std::vector<o2::MCCompLabel> &allMcLabels) {

    for (size_t i = 0; i < allMatchInfo.size(); ++i) {
      const auto &matchInfo = allMatchInfo[i];
      const auto &mcMatchInfo = allMcLabels[i];

      int eventNum = matchInfo.getEventNumberFromTrack();
      int chamberNum = matchInfo.getChamber();

      matchInfoByEventChamber[eventNum][chamberNum].push_back(matchInfo);
      mcMatchInfoByEventChamber[eventNum][chamberNum].push_back(mcMatchInfo);
      // LOGP(info, "pushing back match eventNum {} chamberNum {} ", eventNum,
      // chamberNum);
    }
  }

  std::vector<o2::hmpid::Cluster> getClusters(int eventNum, int chamberNum) {
    return clustersByEventChamber[eventNum][chamberNum];
  }

  std::vector<o2::dataformats::MatchInfoHMP> getMatchInfo(int eventNum,
                                                          int chamberNum) {
    return matchInfoByEventChamber[eventNum][chamberNum];
  }

  std::vector<o2::MCCompLabel> getMcMatchInfo(int eventNum, int chamberNum) {
    return mcMatchInfoByEventChamber[eventNum][chamberNum];
  }

  void setTriggers(std::vector<o2::hmpid::Trigger> *trigArrPtr) {
    if (!trigArrPtr) {
      throw std::runtime_error("trigArrPtr nullptr");
      return;
    }
    for (const auto trig : *trigArrPtr) {
      triggers.emplace_back(trig);
      LOGP(info, "pbck trig w size {}", trig.getNumberOfObjects());
    }
  }

  float calcMassFromCkov(float p, float n, float ckov) {
    p = std::abs(p);
    const float refIndexFreon = n;

    const float cos_ckov = std::cos(ckov);

    const float term = n * p * cos_ckov;
    float m_squared = term * term - p * p;

    // Sanity check to avoid taking the square root of a negative number
    if (m_squared < 0) {
      return 0;
    }

    return std::sqrt(m_squared);
  }

  float calcCkovFromMass(float p, float n, int pdg) {
    // Constants for particle masses (in GeV/c^2)
    const float mass_Muon = 0.10566, mass_Pion = 0.1396, mass_Kaon = 0.4937,
                mass_Proton = 0.938;
    auto particlePDG = TDatabasePDG::Instance()->GetParticle(pdg);
    double mass = particlePDG ? particlePDG->Mass() : 0.;

    float m; // variable to hold the mass
    p = std::abs(p);
    // Switch based on absolute value of PDG code
    switch (std::abs(pdg)) {
    case 13:
      m = mass_Muon;
      break;
    case 211:
      m = mass_Pion;
      break;
    case 321:
      m = mass_Kaon;
      break;
    case 2212:
      m = mass_Proton;
      break;
    default:
      return mass; // return 0 if PDG code doesn't match any known codes
    }

    const float p_sq = p * p;
    const float refIndexFreon = n; // Assuming n is the refractive index
    const float cos_ckov_denom = p * refIndexFreon;

    // Sanity check
    if (p_sq + m * m < 0) {
      return 0;
    }

    const auto cos_ckov =
        static_cast<float>(TMath::Sqrt(p_sq + m * m) / cos_ckov_denom);

    // Sanity check
    if (cos_ckov > 1 || cos_ckov < -1) {
      return 0;
    }
    const auto ckovAngle = static_cast<float>(TMath::ACos(cos_ckov));
    return ckovAngle;
  }

  void writeClusterAttributes(const MCRecon &mcReconObj,
                              Tracks::ClusterData &clusterData) {
    clusterData.pionProbs = mcReconObj.getProtonProb();
    clusterData.kaonProbs = mcReconObj.getKaonProb();
    clusterData.protonProbs = mcReconObj.getPionProb();

    clusterData.pionProbsNorm = mcReconObj.getProtonProbNorm();
    clusterData.kaonProbsNorm = mcReconObj.getKaonProbNorm();
    clusterData.protonProbsNorm = mcReconObj.getPionProbNorm();

    // Positional and cerenkov angle values
    clusterData.xValues = mcReconObj.getXValues();
    clusterData.yValues = mcReconObj.getYValues();
    clusterData.thetaCerValues = mcReconObj.getThetaCerValues();
    clusterData.phiCerValues = mcReconObj.getPhiCerValues();
    clusterData.qValues = mcReconObj.getQValues();
    clusterData.sigmaRingValues = mcReconObj.getSigmaRingValues();
    clusterData.sumProbabilityTrack = mcReconObj.getSumProbTrack();
    clusterData.sizeValues = mcReconObj.getSizeValues();

    clusterData.rawSizeValues = mcReconObj.getRawSizeValues();
    clusterData.numRawClustersValues = mcReconObj.getNumClustersValues();


    /* ef > verbose printing
    int positiveXCount = 0;
    for (int i = 0; i < clusterData.xValues.size(); i++) {
      if (clusterData.xValues[i] > 0) {
        ++positiveXCount; // Increment counter when condition is met

        printf(
            "PROBS %.1f %.1f %.1f %.1f %.1f %.1f || x,y,th,phi  %.1f %.1f %.1f "
            "%.1f  ||  q si sumProb size %.1f %.4f %.1f %d || rawSize "
            "numRawClu %d %d  ",
            clusterData.pionProbs[i], clusterData.kaonProbs[i],
            clusterData.protonProbs[i], clusterData.pionProbsNorm[i],
            clusterData.kaonProbsNorm[i], clusterData.protonProbsNorm[i],
            clusterData.xValues[i], clusterData.yValues[i],
            clusterData.thetaCerValues[i], clusterData.phiCerValues[i],
            clusterData.qValues[i], clusterData.sigmaRingValues[i],
            clusterData.sumProbabilityTrack[i], // Assuming this is an array
                                                // like the others
            clusterData
                .sizeValues[i], // Assuming this is an array like the others
            clusterData
                .rawSizeValues[i], // Assuming this is an array like the others
            clusterData.numRawClustersValues[i]); // Assuming this is an array
                                                  // like the others

        if (positiveXCount >= 3) {
          break; // Exit loop after clusterData.xValues[i] > 0 has occurred 3
                 // times
        }
      }
    } */
  }

  void writeTrackAttributes(const Tracks::Tracks &tracks,
                            int indexOfReconObject) {
    for (size_t i = 0; i < tracks.xMips.size(); ++i) {
      if (i == indexOfReconObject) {
        // Assign values to "this track" attributes
        trackAttributes.xMipThisTrack = tracks.xMips[i];
        trackAttributes.yMipThisTrack = tracks.yMips[i];
        trackAttributes.xRadThisTrack = tracks.xRads[i];
        trackAttributes.yRadThisTrack = tracks.yRads[i];
        trackAttributes.xPCThisTrack = tracks.xPCs[i];
        trackAttributes.yPCThisTrack = tracks.yPCs[i];
        trackAttributes.thetaPThisTrack = tracks.thetaPs[i];
        trackAttributes.phiPThisTrack = tracks.phiPs[i];
        trackAttributes.momentumThisTrack = tracks.momentums[i];
        trackAttributes.qMipThisTrack = tracks.qMips[i];
        trackAttributes.sizeMipThisTrack = tracks.sizeMips[i];
        trackAttributes.mipPcDistThisTrack = tracks.mipPcDists[i];
        trackAttributes.ckovThPionThisTrack = tracks.ckovThPion[i];
        trackAttributes.ckovThKaonThisTrack = tracks.ckovThKaon[i];
        trackAttributes.ckovThProtonThisTrack = tracks.ckovThProton[i];
        trackAttributes.refIndexThisTrack = tracks.refIndexes[i];

        LOGP(info, " writeTrackAttributes tracks.xRads[i] {}",
             tracks.xRads[i]); // Append values to "other tracks" vectors

        trackAttributes.ckovReconMassHypThisTrack = tracks.ckovReconMassHyp[i];
        trackAttributes.ckovReconThisTrack = tracks.ckovRecon[i];

        // number of selected photons
        trackAttributes.numCkovHough = tracks.numCkovHough[i];
        trackAttributes.numCkovHoughMH = tracks.numCkovHoughMH[i];

      } else {

        LOGP(info, " writeTrackAttributes tracks.xRads[i] {}",
             tracks.xRads[i]); // Append values to "other tracks" vectors
        trackAttributes.xMipsOtherTracks.push_back(tracks.xMips[i]);
        trackAttributes.yMipsOtherTracks.push_back(tracks.yMips[i]);
        trackAttributes.xRadsOtherTracks.push_back(tracks.xRads[i]);
        trackAttributes.yRadsOtherTracks.push_back(tracks.yRads[i]);
        trackAttributes.xPCsOtherTracks.push_back(tracks.xPCs[i]);
        trackAttributes.yPCsOtherTracks.push_back(tracks.yPCs[i]);
        trackAttributes.thetaPsOtherTracks.push_back(tracks.thetaPs[i]);
        trackAttributes.phiPsOtherTracks.push_back(tracks.phiPs[i]);
        trackAttributes.momentumsOtherTracks.push_back(tracks.momentums[i]);
        trackAttributes.qMipsOtherTracks.push_back(tracks.qMips[i]);
        trackAttributes.sizeMipsOtherTracks.push_back(tracks.sizeMips[i]);
        trackAttributes.mipPcDistsOtherTracks.push_back(tracks.mipPcDists[i]);
        trackAttributes.ckovThPionOtherTracks.push_back(tracks.ckovThPion[i]);
        trackAttributes.ckovThKaonOtherTracks.push_back(tracks.ckovThKaon[i]);
        trackAttributes.ckovThProtonOtherTracks.push_back(
            tracks.ckovThProton[i]);
        trackAttributes.refIndexesOtherTracks.push_back(tracks.refIndexes[i]);

        trackAttributes.ckovReconMassHypOtherTracks.push_back(
            tracks.ckovReconMassHyp[i]);
        trackAttributes.ckovReconOtherTracks.push_back(tracks.ckovRecon[i]);
      }
    }
  }

private:
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> cluLblArr,
      *cluLblArrPtr = &cluLblArr;

  std::vector<o2::hmpid::Cluster> mClusters;
  std::vector<o2::hmpid::Trigger> triggers;
  EventChamberMatchInfoMap matchInfoByEventChamber;
  EventChamberMCLabelMap mcMatchInfoByEventChamber;
  EventChamberClustersMap clustersByEventChamber;
};
#endif
