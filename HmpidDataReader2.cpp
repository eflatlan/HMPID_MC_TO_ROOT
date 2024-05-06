#ifndef HMPID_DATA_READER2_H
#define HMPID_DATA_READER2_H

#include "HmpidDataSorter2.cpp"

#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h> // gRoot

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <fstream>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"

// C++ header files and libraries
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <gsl/gsl>
#include <iostream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;

using Clusters = o2::hmpid::Cluster;
using Cluster = o2::hmpid::Cluster; //, o2::hmpid::Digit, o2::hmpid::Trigger,
                                    // o2::hmpid::Clusterer;

using MCLabelContainer = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;

class HmpidDataReader2 {
private:
  bool mUseMc = true;

  TFile fileKine;
  TFile fileClu;
  TFile fileMatch;

  TTree *treeMatch;
  TTree *treeClu;
  TTree *treeKine;

  // mc vector : holdeing MC truth-information from ITS-TPC
  // entries read by using label.getEventID() etc on elems in  mLabelHMP vector
  std::vector<o2::MCTrack> mcTrackArr, *mcTrackArrPtr = &mcTrackArr;

  // Match information from HMPID
  std::vector<o2::MCCompLabel> mLabelHMP, *mLabelHMPPtr = &mLabelHMP;
  std::vector<o2::dataformats::MatchInfoHMP> mMatches, *matchArrPtr = &mMatches;

  // Cluster information from HMPID
  std::vector<Cluster> cluArr, *cluArrPtr = &cluArr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> cluLblArr,
      *cluLblArrPtr = &cluLblArr;
  // unsure of implementation and usage of MC truth here, not strictly needed

  // Trigger information from HMPID
  std::vector<Trigger> *trigArr = nullptr;

public:
  HmpidDataReader2(const char *matchFileName, const char *cluFileName,
                   const char *mcFileName)
      : fileKine(mcFileName, "READ"), fileClu(cluFileName, "READ"),
        fileMatch(matchFileName, "READ") {
    // Check if files are opened correctly

    bool mUseMC = true;

    /*if (fileClu.IsZombie()) {
      LOGP(error, "Failed to open cluster File: {}", cluFileName);
      throw std::runtime_error("Failed to open cluFileName.");
    } else {
      LOGP(info, "cluster File opened successfully: {}", cluFileName);
    }*/

    if (fileMatch.IsZombie()) {
      LOGP(info, "Failed to open match file: {}", matchFileName);
      // throw std::runtime_error("Failed to open matchFileName.");
    } else {
      LOGP(info, "Match file opened successfully: {}", matchFileName);
    }

    // ok : we use mcReader instead
    /*if(mUseMC) {
      if (fileKine.IsZombie()) {
        LOGP(error, "Failed to open mc truth file: {}", mcFileName);
        throw std::runtime_error("Failed to open mcFileName.");
      } else {
        LOGP(info, "Failed to open mc truth file:  {}", mcFileName);
      }
    } else {
      LOGP(info, "Mc truth false");
    }*/

    // Initialize TTree objects
    initializeClusterTree();
    initializeMatchTree(/*0, 0, 0*/);

    /* use mcReader instead to read MC truth information from ITS-TPC
    if(mUseMc) {
      initializeMCTree(mcTrackArrPtr);
    } */

    if (!treeMatch) {
      LOGP(error, "Failed to initialize treeMatch TTree object.");
      throw std::runtime_error("Failed to initialize treeMatch TTree object.");
    } else {
      LOGP(info, "treeMatch TTree object initialized successfully.");
    }

    if (!treeClu) {
      LOGP(error, "Failed to initialize treeClu TTree object.");
      throw std::runtime_error("Failed to initialize treeClu TTree object.");
    } else {
      LOGP(info, "treeClu TTree object initialized successfully.");
    }

    /*
    if (mUseMc) {
      if (!treeKine) {
        LOGP(error, "Failed to initialize treeKine TTree object.");
        throw std::runtime_error("Failed to initialize treeKine TTree object.");
      } else {
        LOGP(info, "treeKine TTree object initialized successfully.");
      }
    } else {
      LOGP(info, "Tree MC truth not initialized because Mc truth false");
    }*/

    // mcTrackArrPtr
    // bool mUseMc = true;
    if (mUseMC) {
      LOGP(info, "HmpidDataReader2 : useMC {} ; call hmpidDataSorter2", mUseMC);
      HmpidDataSorter2 hmpidDataSorter2;

      hmpidDataSorter2.setTriggers(trigArr /*, trigArr*/);

      if(mMatches.size() != mLabelHMP.size()) {
        LOGP(info, "mMatchessize {} mLabelHMP size {} ", mMatches.size(), mLabelHMP.size());
        LOGP(error, "Matches tracks and  MC have different sizes");
        throw std::runtime_error("Matches tracks and  MC have different sizes");
        return;
      }
      
      LOGP(info, "HmpidDataReader2 : useMC {} ; call setClusterMcTruth", mUseMC);
      if(cluArr.size() != cluLblArr.getIndexedSize()) {
        LOGP(info, "cluArr {} cluLblArr size {} ", mMatches.size(), cluLblArr.getIndexedSize());
        LOGP(error, "Clusters and  MC have different sizes");
        throw std::runtime_error("Clusters  and  MC have different sizes");
        return;
      }


      hmpidDataSorter2.organizeAndSortClusters(cluArr);

      hmpidDataSorter2.organizeAndSortMatchInfo(mMatches, mLabelHMP);

      hmpidDataSorter2.setClusterMcTruth(cluLblArr);

      hmpidDataSorter2.iterateOverMatchedTracks();
    }
  }

  ~HmpidDataReader2() {}

  std::vector<Cluster> getClusInEvent(int event) const {
    auto pTgr = &trigArr->at(event);
    std::vector<Cluster> oneEventClusters;
    for (int j = pTgr->getFirstEntry(); j <= pTgr->getLastEntry(); j++) {
      const auto &clu = static_cast<o2::hmpid::Cluster>(cluArrPtr->at(j));
      oneEventClusters.push_back(clu);
    }

    return oneEventClusters;
  }

  std::vector<Cluster> *getcluArrPtr() const { return cluArrPtr; }

  std::vector<Trigger> *getTrigArr() const { return trigArr; }

  void initializeClusterTree();

  void initializeMatchTree();
  void readMatch(int eventID, int &startIndex,
                 std::vector<o2::dataformats::MatchInfoHMP> &filteredMatches,
                 std::vector<o2::MCCompLabel> &filteredLblMatches);

  // std::vector<o2::dataformats::MatchInfoHMP> readMatch(int eventID, int
  // &startIndex);
  /*static TTree *initializeClusterTree(
      std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr);*/

  TTree *initializeMCTree();

  std::vector<o2::MCTrack> *readMC(std::vector<o2::MCTrack> *&mcTrackArrPtr,
                                   TTree *treeKine, int eventId);

  static const o2::MCTrack *getMCEntry(std::vector<o2::MCTrack> *mcTrackArrPtr,
                                       int trackID);
  static void readTreeEntries();
};

void HmpidDataReader2::initializeMatchTree() {

  // std::unique_ptr<TFile> fileMatch(TFile::Open("o2match_hmp.root", "READ"));

  treeMatch = dynamic_cast<TTree *>(fileMatch.Get("matchHMP"));

  if (!treeMatch)
    throw std::runtime_error("treeMatch nullptr");

  // std::vector<o2::dataformats::MatchInfoHMP>* matchArrPtr = nullptr;
  treeMatch->SetBranchAddress("HMPMatchInfo", &matchArrPtr);
  bool mUseMC = true;
  if (mUseMC) {
    treeMatch->SetBranchAddress("MatchHMPMCTruth", &mLabelHMPPtr);
  }

  treeMatch->GetEntry(0);
  treeMatch->Print("toponly");

  Printf("treeMatch entries %lld", treeMatch->GetEntries());

  if (matchArrPtr == nullptr) {
    Printf("HmpidDataReader2::initializeMatchTree matchArrPtr== nullptr");
    // fileMatch.Close();
    throw std::runtime_error("matchArrPtr nullptr");
  }

  if (mUseMC && !mLabelHMPPtr) {
    Printf("HmpidDataReader2::initializeMatchTree mLabelHMPPtr== nullptr");
    // fileMatch.Close();
    throw std::runtime_error("mLabelHMPPtr nullptr");
  } else if (mUseMC && mLabelHMPPtr) {
    LOGP(info, "Match-tracks        : entries {}", matchArrPtr->size());
    LOGP(info, "Got MC-Match-labels : entries {}", mLabelHMPPtr->size());
  }
}

// eventId = eventID to be searched for;
// startIndex : index of where matchArrPtr is to be searched
// newStartIndex startIndex for next event

void HmpidDataReader2::readMatch(
    int eventID, int &startIndex,
    std::vector<o2::dataformats::MatchInfoHMP> &filteredMatches,
    std::vector<o2::MCCompLabel> &filteredLblMatches) {

  Printf("Call HmpidDataReader2::readMatch");
  // std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new
  // std::vector<o2::dataformats::MatchInfoHMP>;

  if (!treeMatch) {
    Printf("TTree not initialized");
    return; // filteredMatches;
  } else {
    Printf("HmpidDataReader2::readMatch : TTree  initialized");
  }

  // Prepare to store filtered matches
  // std::vector<o2::dataformats::MatchInfoHMP> filteredMatches;// = new
  // std::vector<o2::dataformats::MatchInfoHMP>; tracks should be stored in
  // "time" --> when we find our event we can then switch this condition "of"
  // when the event changes:

  bool found = false;

  if (matchArrPtr == nullptr) {
    Printf("HmpidDataReader2::readMatch :: matchArrPtr== nullptr");
    return; // filteredMatches;
  } else {
    Printf("HmpidDataReader2::readMatch : matchArrPtr ok");
  }

  Printf("readMatch : (*matchArrPtr) size : %zu ", (*matchArrPtr).size());
  Printf("readMatch : startIndex %d", startIndex);

  if ((*matchArrPtr).size() < 1) {
    LOGP(info, "matchArrPtr was 0");
    return; // filteredMatches;
  }

  Printf("readMatch : (*matchArrPtr)[startIndex].getEvent() %d eventID %d",
         (*matchArrPtr)[startIndex].getEvent(), eventID);

  if ((*matchArrPtr)[startIndex].getEvent() != eventID) {
    Printf("This shouldnt happen");
    LOGP(info, "(*matchArrPtr)[startIndex].getEvent() {} eventID {}",
         (*matchArrPtr)[startIndex].getEvent(), eventID);
  } else {
    found = true;
  }

  Printf("readMatch : eventID %d startIndex %d", eventID, startIndex);
  LOGP(info, "matchArrPtr->size() {}", matchArrPtr->size());
  LOGP(info, "mLabelHMPPtr->size() {}", mLabelHMPPtr->size());

  // int eventID, int &startIndex
  for (int i = startIndex; i < matchArrPtr->size(); i++) {

    // hvorfor [1] her ?
    const auto &track = (*matchArrPtr)[i];
    const auto &lblTrk = (*mLabelHMPPtr)[i];

    // filteredMatches->push_back(track);
    // startIndex = i;
    //  ef : we cant use this atm, since the clusters from the same trigger
    //  sometimes have different eventNumber!

    LOGP(info, "trackEvent {}, i {} | getMipClusEvent {}", track.getEvent(), i,
         track.getMipClusEvent());

    // ef: TODO is this now needed when checking MCtruth from MClabel?
    if (track.getEvent() != eventID) {
      startIndex = i; // new startIndex for next event
      Printf("readMatch : eventID changed - %d; end of loop ",
             track.getEvent());

      break;
    } else {

      filteredMatches.push_back(track);
      filteredLblMatches.push_back(lblTrk);
    }
  }

  Printf("readMatch : new startIndex %d", startIndex);
}

void HmpidDataReader2::initializeClusterTree(
    /*std::vector<Cluster> *&cluArrPtr, std::vector<Trigger> *&trigArr*/) {
  // TTree *treeClu = (TTree *)fileClu.Get("o2sim");
  treeClu = dynamic_cast<TTree *>(fileClu.Get("o2hmp"));
  if (!treeClu)
    treeClu = dynamic_cast<TTree *>(fileClu.Get("o2sim"));

  if (!treeClu) {
    Printf("Error accessing TTree");
  }
  treeClu->Print("toponly");

  treeClu->SetBranchAddress("HMPIDclusters", &cluArrPtr);
  treeClu->SetBranchAddress("InteractionRecords", &trigArr);

  // ef : add MC clus information
  bool mUseMC = true;
  if (treeClu->GetBranchStatus("HMPIDClusterLabels") && mUseMC) {
    treeClu->SetBranchAddress("HMPIDClusterLabels", &cluLblArrPtr);
    LOGP(info, "HMPIDClusterLabels adress");
  }

  treeClu->GetEntry(0);

  if (!cluLblArrPtr && mUseMC) {
    Printf("HmpidDataReader2::initializeClusterTree cluLblArrPtr== nullptr");
    throw std::runtime_error("cluLblArrPtr nullptr");
  } else if (cluLblArrPtr && mUseMC) {
    LOGP(info, "triggers          : entries {}", trigArr->size());
    LOGP(info, "Got MC-clu-labels : entries {}",
         cluLblArrPtr->getIndexedSize());
    LOGP(info, "clusters          :  entries {}", cluArrPtr->size());
  }

  Printf("treeClu entries %lld", treeClu->GetEntries());
}

TTree *HmpidDataReader2::initializeMCTree() {

  treeKine = dynamic_cast<TTree *>(fileKine.Get("o2sim"));
  // if (!treeKine)
  //   treeKine = dynamic_cast<TTree*>(fileKine.Get("o2sim"));

  if (!treeKine) {
    Printf("Error accessing TTree");
    /*fileKine->Close();
    delete fileKine;*/
    return nullptr;
  }

  std::vector<o2::MCTrack> *mcTrackArrPtr;
  treeKine->SetBranchAddress("MCTrack", &mcTrackArrPtr);
  // treeKine->SetBranchAddress("TrackRefs", &mcTrackArrPtr);

  bool mUseMC = true;
  if (/*!mcTrackArrPtr &&*/ mUseMC) {
    Printf("HmpidDataReader2::initializeMCTree mcTrackArrPtr== nullptr");
    throw std::runtime_error("mcTrackArrPtr nullptr");
  }

  treeKine->GetEntry(0);
  treeKine->Print("toponly");
  Printf("treeKine entries %lld", treeKine->GetEntries());

  return treeKine;
}

std::vector<o2::MCTrack> *
HmpidDataReader2::readMC(std::vector<o2::MCTrack> *&mcTrackArrPtr,
                         TTree *treeKine, int eventID) {
  if (treeKine == nullptr) {
    Printf("Error : treeKine == nullptr");
    return nullptr;
  }

  if (eventID < treeKine->GetEntries()) {
    treeKine->GetEntry(eventID);
    Printf("readMC at entry %d", eventID);
    return mcTrackArrPtr;
  } else {
    Printf("eventId > treeKine->GetEntries()");
    return nullptr;
  }
}

// for the given eventID; read trackID
const o2::MCTrack *
HmpidDataReader2::getMCEntry(std::vector<o2::MCTrack> *mcTrackArrPtr,
                             int trackID) {

  if (trackID < 0 || trackID >= mcTrackArrPtr->size()) {
    return nullptr;
  } else {
    // const auto& track = mcTrackArrPtr->at(trackID);
    for (int i = 0; i < mcTrackArrPtr->size(); ++i) {
      const auto &mcTrack = (*mcTrackArrPtr)[i];
      if (i == trackID) {

        /*TParticlePDG* pPDG =
        TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdgCode()); auto pT =
        TMath::Abs(mcTrack->GetStartVertexMomentumX() *
                              mcTrack->GetStartVertexMomentumX() +
                          mcTrack->GetStartVertexMomentumY() *
                              mcTrack->GetStartVertexMomentumY());
        auto pdgCode = mcTrack->GetPdgCode();
        LOGP(info, "Particle {}: pdg = {}, pT = {}", i, pdgCode, pT);
        */
        /*
        Printf("Particle %d: pdg = %d, pT = %f, px = %f, py = %f, pz = %f, vx "
               "= %f, vy = %f, vz = %f",
               i, mcTrack->GetPdgCode(),
               TMath::Abs(mcTrack->GetStartVertexMomentumX() *
                              mcTrack->GetStartVertexMomentumX() +
                          mcTrack->GetStartVertexMomentumY() *
                              mcTrack->GetStartVertexMomentumY()),
               mcTrack->GetStartVertexMomentumX(),
               mcTrack->GetStartVertexMomentumY(),
               mcTrack->GetStartVertexMomentumZ(),
               mcTrack->GetStartVertexCoordinatesX(),
               mcTrack->GetStartVertexCoordinatesY(),
               mcTrack->GetStartVertexCoordinatesZ());*/
        return &mcTrack;
      }
    }
  }
  return nullptr;
}

void HmpidDataReader2::readTreeEntries() {
  // Open the ROOT file
  /*
  int i;
  auto trackVector = readTrack(0,0,i);
  auto clusterVector = readClu(0,0,i);
  std::vector<o2::MCTrack>* mcVector = readMC(0,0,i);
  for() {

  }

        Printf("numTracks %d numClusters %d numMC %d", trackVector->size(),
  clusterVector->size(), mcVector->size()); */
}

#endif // HMPID_DATA_READER_H
