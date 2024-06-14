
#include "MCRecon.cpp"
#include <vector>

#ifndef Tracks_H
#define Tracks_H

namespace Tracks // fix namespace etc..
{

struct Tracks {
  std::vector<float> xMips, yMips, xRads, yRads, xPCs, yPCs, thetaPs, phiPs,
      momentums, mipPcDists, ckovThPion, ckovThKaon, ckovThProton, refIndexes;
  std::vector<int> sizeMips, qMips;

  std::vector<int> numCkovHough, numCkovHoughMH;

  std::vector<float> ckovRecon, ckovReconMassHyp;
  void addTrackCands(MCRecon *mcReconObj) {
    if (mcReconObj) {

      LOGP(info, " xRads.empl getRadX xRad {}", mcReconObj->getRadX());
      xMips.emplace_back(mcReconObj->getMipX());
      yMips.emplace_back(mcReconObj->getMipY());
      xRads.emplace_back(mcReconObj->getRadX());
      yRads.emplace_back(mcReconObj->getRadY());
      xPCs.emplace_back(mcReconObj->getXpc());
      yPCs.emplace_back(mcReconObj->getYpc());
      thetaPs.emplace_back(mcReconObj->getThetaP());
      phiPs.emplace_back(mcReconObj->getPhiP());
      momentums.emplace_back(mcReconObj->getMomentum());
      qMips.emplace_back(mcReconObj->getQMip());
      sizeMips.emplace_back(mcReconObj->getSizeMip());
      mipPcDists.emplace_back(mcReconObj->getMipPcDist());
      ckovThPion.emplace_back(mcReconObj->getCkovThPion());
      ckovThKaon.emplace_back(mcReconObj->getCkovThKaon());
      ckovThProton.emplace_back(mcReconObj->getCkovThProton());
      refIndexes.emplace_back(mcReconObj->getRefIndex());

      ckovRecon.emplace_back(mcReconObj->getCkovRecon());
      
      // only for expansive ML datataking, not part of official MC HMPID
      //ckovReconMassHyp.emplace_back(mcReconObj->getCkovReconMassHyp());

      int ckovHough = mcReconObj->getNumCkovHough();

      // only for expansive ML datataking, not part of official MC HMPID
      //int ckovHoughMH = mcReconObj->getNumCkovHoughMH();

      numCkovHough.push_back(ckovHough);
      
      // only for expansive ML datataking, not part of official MC HMPID      
      // numCkovHoughMH.push_back(ckovHoughMH);

      // Print the values just added
      std::cout << "Tracks :: Number of Cerenkov photons detected: "
                << ckovHough << std::endl;

    }
  }
};

struct TrackAttributes {
  // Attributes for "this track"
  float xMipThisTrack;
  float yMipThisTrack;
  float xRadThisTrack;
  float yRadThisTrack;
  float xPCThisTrack;
  float yPCThisTrack;
  float thetaPThisTrack;
  float phiPThisTrack;
  float momentumThisTrack;
  int qMipThisTrack;
  int sizeMipThisTrack;
  float mipPcDistThisTrack;
  float ckovThPionThisTrack;
  float ckovThKaonThisTrack;
  float ckovThProtonThisTrack;
  float refIndexThisTrack;

  float ckovReconThisTrack;
  int numCkovHough; // number of Hough selected photons

  // only for expansive ML datataking, not part of official MC HMPID
  //int numCkovHoughMH; // number of Hough selected photons
  //float ckovReconMassHypThisTrack;

  



  // Vectors for "other tracks"
  std::vector<float> xMipsOtherTracks;
  std::vector<float> yMipsOtherTracks;
  std::vector<float> xRadsOtherTracks;
  std::vector<float> yRadsOtherTracks;
  std::vector<float> xPCsOtherTracks;
  std::vector<float> yPCsOtherTracks;
  std::vector<float> thetaPsOtherTracks;
  std::vector<float> phiPsOtherTracks;
  std::vector<float> momentumsOtherTracks;
  std::vector<int> qMipsOtherTracks;
  std::vector<int> sizeMipsOtherTracks;
  std::vector<float> mipPcDistsOtherTracks;
  std::vector<float> ckovThPionOtherTracks;
  std::vector<float> ckovThKaonOtherTracks;
  std::vector<float> ckovThProtonOtherTracks;
  std::vector<float> refIndexesOtherTracks;

  std::vector<float> ckovReconOtherTracks;
  
  // number of Hough selected photons  
  //std::vector<float> ckovReconMassHypOtherTracks;

  void clear() {
    xMipsOtherTracks.clear();
    yMipsOtherTracks.clear();
    xRadsOtherTracks.clear();
    yRadsOtherTracks.clear();
    xPCsOtherTracks.clear();
    yPCsOtherTracks.clear();
    thetaPsOtherTracks.clear();
    phiPsOtherTracks.clear();
    momentumsOtherTracks.clear();
    mipPcDistsOtherTracks.clear();
    ckovThPionOtherTracks.clear();
    ckovThKaonOtherTracks.clear();
    ckovThProtonOtherTracks.clear();
    refIndexesOtherTracks.clear();
    qMipsOtherTracks.clear();
    sizeMipsOtherTracks.clear();
    ckovReconOtherTracks.clear();
    //ckovReconMassHypOtherTracks.clear();
  }

  void setNumberOfTracks(int nTracks) {

    // nTracks - 1 bc we take only the other tracks,
    //

    if (nTracks > 1) {
      nTracks = nTracks - 1;
    }

    xMipsOtherTracks.clear();
    yMipsOtherTracks.clear();
    xRadsOtherTracks.clear();
    yRadsOtherTracks.clear();
    xPCsOtherTracks.clear();
    yPCsOtherTracks.clear();
    thetaPsOtherTracks.clear();
    phiPsOtherTracks.clear();
    momentumsOtherTracks.clear();
    mipPcDistsOtherTracks.clear();
    ckovThPionOtherTracks.clear();
    ckovThKaonOtherTracks.clear();
    ckovThProtonOtherTracks.clear();
    refIndexesOtherTracks.clear();
    qMipsOtherTracks.clear();
    sizeMipsOtherTracks.clear();
    ckovReconOtherTracks.clear();
    //ckovReconMassHypOtherTracks.clear();

    xMipsOtherTracks.reserve(nTracks);
    yMipsOtherTracks.reserve(nTracks);
    xRadsOtherTracks.reserve(nTracks);
    yRadsOtherTracks.reserve(nTracks);
    xPCsOtherTracks.reserve(nTracks);
    yPCsOtherTracks.reserve(nTracks);
    thetaPsOtherTracks.reserve(nTracks);
    phiPsOtherTracks.reserve(nTracks);
    momentumsOtherTracks.reserve(nTracks);
    mipPcDistsOtherTracks.reserve(nTracks);
    ckovThPionOtherTracks.reserve(nTracks);
    ckovThKaonOtherTracks.reserve(nTracks);
    ckovThProtonOtherTracks.reserve(nTracks);
    refIndexesOtherTracks.reserve(nTracks);
    qMipsOtherTracks.reserve(nTracks);
    sizeMipsOtherTracks.reserve(nTracks);
    ckovReconOtherTracks.reserve(nTracks);
    //ckovReconMassHypOtherTracks.reserve(nTracks);
  }
};

struct ClusterData {
  int numEntries; // Store the number of entries for resetting

  std::vector<float> xValues;
  std::vector<float> yValues;
  std::vector<float> qValues;
  std::vector<float> thetaCerValues;
  std::vector<float> phiCerValues;
  std::vector<float> sigmaRingValues;
  std::vector<float> pionProbs;
  std::vector<float> kaonProbs;
  std::vector<float> protonProbs;
  std::vector<float> protonProbsNorm;
  std::vector<float> kaonProbsNorm;
  std::vector<float> pionProbsNorm;
  std::vector<float> sumProbabilityTrack;
  std::vector<int> sizeValues;


  // not part of official MC HMPID
  //std::vector<int> rawSizeValues;
  //std::vector<int> numRawClustersValues;

  ClusterData() {}

  // Method to clear and reset all vectors
  void clear() {
    xValues.assign(numEntries, 0.0f);
    yValues.assign(numEntries, 0.0f);
    qValues.assign(numEntries, 0.0f);
    thetaCerValues.assign(numEntries, 0.0f);
    phiCerValues.assign(numEntries, 0.0f);
    sigmaRingValues.assign(numEntries, 0.0f);
    pionProbs.assign(numEntries, 0.0f);
    kaonProbs.assign(numEntries, 0.0f);
    protonProbs.assign(numEntries, 0.0f);
    protonProbsNorm.assign(numEntries, 0.0f);
    kaonProbsNorm.assign(numEntries, 0.0f);
    pionProbsNorm.assign(numEntries, 0.0f);
    sumProbabilityTrack.assign(numEntries, 0.0f);
    sizeValues.assign(numEntries, 0);
    //rawSizeValues.assign(numEntries, 0);
    //numRawClustersValues.assign(numEntries, 0);
  }

  void setNumberOfEntries(int nEntries) {
    xValues.resize(nEntries, 0.0f);
    yValues.resize(nEntries, 0.0f);
    qValues.resize(nEntries, 0.0f);
    thetaCerValues.resize(nEntries, 0.0f);
    phiCerValues.resize(nEntries, 0.0f);
    sigmaRingValues.resize(nEntries, 0.0f);
    pionProbs.resize(nEntries, 0.0f);
    kaonProbs.resize(nEntries, 0.0f);
    protonProbs.resize(nEntries, 0.0f);
    protonProbsNorm.resize(nEntries, 0.0f);
    kaonProbsNorm.resize(nEntries, 0.0f);
    pionProbsNorm.resize(nEntries, 0.0f);
    sumProbabilityTrack.resize(nEntries, 0.0f);
    sizeValues.resize(nEntries, 0);
    //rawSizeValues.resize(nEntries, 0);
    //numRawClustersValues.resize(nEntries, 0);
  }
};

// struct  holding clusters that exceeds a certain charge-trheshold
// for evaluating FB photons
struct HighChargeClu {

  void addToFields(float _x, float _y, float _q, int _size) {
    x.push_back(_x);
    y.push_back(_y);
    q.push_back(_q);
    size.push_back(_size);
  }

  std::vector<float> x, y, q;
  std::vector<int> size;

  void setNumberOfEntries(int numberOfEntries) { // not using resize here bc we
                                                 // dont excactly know size
    x.reserve(numberOfEntries);
    y.reserve(numberOfEntries);
    q.reserve(numberOfEntries);
    size.reserve(numberOfEntries);
  }

  void clearDataFields() {
    x.clear();
    y.clear();
    q.clear();
    size.clear();
  }
};

} // namespace Tracks
#endif
