//#if !defined(__CLING__) || defined(__ROOTCLING__)

#pragma once
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/RangeReference.h"

#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "HMPIDBase/Param.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

// #include "CkovToolsSingle.cpp"
//#include "ParticleUtils2.cpp"
#include "HmpidDataReader2.cpp"

#include <utility>
#include <vector>

#include <TMath.h>
#include <cmath>


// #endif


float calcCkovFromMass(float p, float n, int pdg);


void evaluateClusterTrack(
    std::vector<o2::hmpid::ClusterCandidate> &clusterPerChamber,
    const o2::dataformats::MatchInfoHMP &track,
    const std::vector<float> &mipCharges, int mcTrackPdg, int trackNumber,
    int &plotNumber);

std::array<float, 3> calcCherenkovHyp(float p, float n);


void SegmentationCkov2(double _sigmaSep = 1.5) {


  /*
  treeOut->Branch("reconCkov", &reconCkovBranch);
  treeOut->Branch("cluCharge", &cluChargeBranch);
  treeOut->Branch("cluSize", &cluSizeBranch);
  treeOut->Branch("refIndex", &refIndexBranch);
  treeOut->Branch("xRad", &xRadBranch);
  treeOut->Branch("yRad", &yRadBranch);
  treeOut->Branch("xMip", &xMipBranch);
  treeOut->Branch("yMip", &yMipBranch);
  treeOut->Branch("th", &thBranch);
  treeOut->Branch("ph", &phBranch);
  treeOut->Branch("p", &pBranch);
  */ 

  // clusters and triggers
  std::vector<o2::hmpid::Cluster> *clusterArr = nullptr;
  /*std::vector<o2::hmpid::Topology> mTopologyFromFile,
      *mTopologyFromFilePtr = &mTopologyFromFile;*/ 


  auto matchFileName = "o2match_hmp.root";
  auto cluFileName = "hmpidclusters.root";
  auto mcFileName = "o2sim_Kine.root";
	
  
  HmpidDataReader2 hmpidDataReader(matchFileName, cluFileName, mcFileName);


  /*LOGP(info, "fileOut->Close();"); * / 
  mParticleEvents.writeH5();
  LOGP(info, " mParticleEvents.writeH5();");
  */
}

float calcCkovFromMass(float p, float n, int pdg) {
  // Constants for particle masses (in GeV/c^2)
  const float mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938;

  float m; // variable to hold the mass

  // Switch based on absolute value of PDG code
  switch (std::abs(pdg)) {
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
    return 0; // return 0 if PDG code doesn't match any known codes
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

void evaluateClusterTrack(
    std::vector<o2::hmpid::ClusterCandidate> &clusterPerChamber,
    const o2::dataformats::MatchInfoHMP &track,
    const std::vector<float> &mipCharges, int mcTrackPdg, int trackNumber,
    int &plotNumber) {

  const auto eventCnt =
      track.getEvent(); // check it corresponds to entry in loop of events?
  const auto momentum = track.getHmpMom();

  const auto nF =
      1.2929 - 0.0025; // track.getRefIndex(); // ef: aon it only gets the mean
                       // value ; TODO: in MatchHMP.cxx get calibration value
  // https://github.com/eflatlan/AliceO2/blob/811fcc6d00b363b1e96e0aa8269d46eed95d879b/Detectors/GlobalTracking/src/MatchHMP.cxx#L432C28-L432C38
  // https://github.com/AliceO2Group/AliceO2/blob/e6b603e4c92f98733ff9f7954100140e72bd99f6/Detectors/HMPID/base/include/HMPIDBase/Param.h#L159C37-L159C64

  const auto nQ = 1.583;
  const auto nG = 1.0005;

  // make static?

  // make account for varying nF ? +- 2 std ?
  const auto &ckovHypsMin = calcCherenkovHyp(momentum, nF);

  const auto &ckovHypsMax = calcCherenkovHyp(momentum, nF);

  float xRad, yRad, xPc, yPc, thetaP, phiP;
  track.getHMPIDtrk(xRad, yRad, xPc, yPc, thetaP, phiP);

  float xMip, yMip;
  int nph, q;
  track.getHMPIDmip(xMip, yMip, q, nph);

  const auto &L = 0.5; //
  double radParams[7] = {
      xRad,
      yRad,
      L,
      thetaP,
      phiP,
      momentum,
      static_cast<double>(mcTrackPdg * 1.0)}; // ef : TODO add PID to MLinfoHMP?

  double refIndexes[3] = {nF, nQ, nG};

}

const float mass_Pion = 0.1396, mass_Kaon = 0.4937,
            mass_Proton = 0.938; // masses in
std::array<float, 3> masses = {mass_Pion, mass_Kaon, mass_Proton};
const float mass_Pion_sq = mass_Pion * mass_Pion,
            mass_Kaon_sq = mass_Kaon * mass_Kaon,
            mass_Proton_sq = mass_Proton * mass_Proton;


std::array<float, 3> calcCherenkovHyp(float p, float n) {
  const float p_sq = p * p;
  const float cos_ckov_denom = p * n;
  const auto cos_ckov_Pion = static_cast<float>(
      TMath::Sqrt(p_sq + mass_Pion_sq) /
      (cos_ckov_denom)); // n = refIndexFreon 1.289 later make it random?

  const float cos_ckov_Kaon =
      static_cast<float>(TMath::Sqrt(p_sq + mass_Kaon_sq) / (cos_ckov_denom));
  const float cos_ckov_Proton =
      static_cast<float>(TMath::Sqrt(p_sq + mass_Proton_sq) / (cos_ckov_denom));

  const float ckovAnglePion = static_cast<float>(TMath::ACos(cos_ckov_Pion));
  const float ckovAngleKaon = static_cast<float>(TMath::ACos(cos_ckov_Kaon));
  const float ckovAngleProton =
      static_cast<float>(TMath::ACos(cos_ckov_Proton));

  Printf("Pion %.3f Kaon %.3f Proton %.3f", ckovAnglePion, ckovAngleKaon,
         ckovAngleProton);

  return {ckovAnglePion, ckovAngleKaon, ckovAngleProton};
}
