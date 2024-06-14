#include "DataFormatsHMP/Cluster.h"
#include "FilterPhotons.cpp"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include <boost/math/distributions/normal.hpp>

#include "HMPIDBase/Param.h"

#ifndef MCRecon_H
#define MCRecon_H

struct McTruth {
  int pdgCodeTrack = -1;
  int pdgCodeClu = -1;
  bool isMipMatchedCorrectly = false;
  bool isTrackToReconKnownPdg =
      false; // if we take pdg code being known when making the requirements

  int numCountedPhotonsPC =
      -1; // numbef of photon-clusters on PC having the track as MID

  McTruth(bool isMipMatched, int pdgCode, int mipPDG, bool isTrackToRecon,
          int nPc)
      : isMipMatchedCorrectly(isMipMatched), pdgCodeClu(mipPDG),
        pdgCodeTrack(pdgCode), isTrackToReconKnownPdg(isTrackToRecon),
        numCountedPhotonsPC(nPc) {}

  McTruth() {}
};

class MCRecon {

public:

  // cuts for accepting clusters as single-photon candidates
  static constexpr double mThreshMipQ = 200;
  static constexpr int mThreshMipSize = 3;

  // ef > if we take into account PDG code is known
  bool isTrackToReconKnownPdg() { return mIsTrackToReconstructWithKnownPdg; }

  // ef > if we take into account PDG code is known
  bool isTrackToReconKnown() { return mIsTrackToReconstruct; }

  void setTrackReconKnownPdg(bool isTrackToReconstructWithKnownPdg) {
    mIsTrackToReconstructWithKnownPdg = isTrackToReconstructWithKnownPdg;
  }


  // calculate the normalized probability across all species and tracks
  // passed by reference: each time calling this function the sum acroos all species for the given track is added
  void calculateNormalizedProb(std::vector<float> &sumProbAllTracks) {

    /*
    LOGP(info, "sizes of vectors >> pPion {} pKaon {} pProton {} sumProbTrack
    {}", pPion.size(), pKaon.size(), pProton.size(), sumProbTrack.size() );


    std::cout << "\n sumProbTrackIn : ";
    for(size_t i = 0; i < sumProbTrack.size(); ++i) {
        if(sumProbTrack[i] != 0.0) {
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2)
    << sumProbTrack[i] << " ";
        }
    }*/

    // add all fields
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pProton.begin(),
                   sumProbTrack.begin(), std::plus<>());
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pKaon.begin(),
                   sumProbTrack.begin(), std::plus<>());
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pPion.begin(),
                   sumProbTrack.begin(), std::plus<>());

    // Normalize pion probabilities
    std::transform(pPion.begin(), pPion.end(), sumProbTrack.begin(),
                   pNormPion.begin(), [](float pionProb, float sumProb) {
                     return sumProb > 0 ? pionProb / sumProb
                                        : 0.0f; // Avoid division by zero
                   });

    std::transform(pKaon.begin(), pKaon.end(), sumProbTrack.begin(),
                   pNormKaon.begin(), [](float kaonProb, float sumProb) {
                     return sumProb > 0 ? kaonProb / sumProb
                                        : 0.0f; // Avoid division by zero
                   });

    std::transform(pProton.begin(), pProton.end(), sumProbTrack.begin(),
                   pNormProton.begin(), [](float protonProb, float sumProb) {
                     return sumProb > 0 ? protonProb / sumProb
                                        : 0.0f; // Avoid division by zero
                   });

    std::transform(sumProbAllTracks.begin(), sumProbAllTracks.end(),
                   sumProbTrack.begin(), sumProbAllTracks.begin(),
                   std::plus<>());


  }

  void setMcTruth(const McTruth &mcTruth) { mMcTruth = mcTruth; }
  void getMcTruth(McTruth &mcTruth) const { mcTruth = mMcTruth; }


  // track is matched, and has ok number of photons, but not necessarily meeting the tight distance cut
  bool isTrackMatched() const { return mIsTrackMatched; }

  // track is matched,  has ok number of photons, and meets the tight distance cut
  // this track will be reconstructed
  bool isTrackAppropriate() const { return mIsTrackToReconstruct; }


  // position, charge and size for the clusters
  std::vector<float> getXValues() const { return xValues; }
  std::vector<float> getYValues() const { return yValues; }
  std::vector<float> getQValues() const { return qValues; }
  std::vector<int> getSizeValues() const { return sizeValues; }

  // reconstructed Ckov photon attrs for clusters
  std::vector<float> getThetaCerValues() const { return thetaCerValues; }
  std::vector<float> getPhiCerValues() const { return phiCerValues; }
  std::vector<float> getSigmaRingValues() const { return sigmaRingValues; }

  // holding information about the raw clusters
  // deconvoluted clusters stems from raw clusters
  std::vector<int> getNumClustersValues() const { return numRawClustersValues; }  // number of clusters forming raw clusters
  std::vector<int> getRawSizeValues() const { return rawSizeValues; }             // total size of raw cluster



  // Specie probability for the photons
  std::vector<float> getProtonProb() const { return pProton; }
  std::vector<float> getKaonProb() const { return pKaon; }
  std::vector<float> getPionProb() const { return pPion; }


  // normalized porbability per photon: dividing by the the sum of all species for hte photon
  std::vector<float> getProtonProbNorm() const { return pNormProton; }
  std::vector<float> getKaonProbNorm() const { return pNormKaon; }
  std::vector<float> getPionProbNorm() const { return pNormPion; }
  
  // sum across all species and tracks in chamber and event
  std::vector<float> getSumProbTrack() const { return sumProbTrack; }

  double zScoreToProb(double z) {
    boost::math::normal_distribution<> dist(0.0, 1.0); // Mean = 0, Standard Deviation = 1 for a standard normal distribution
    return boost::math::cdf(dist, z); // Cumulative Distribution Function value for z
  }


  // calcualte the mass from hte predicted ckov 
  // ckov will deviate from th ckov, so predicted mass also deviates
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


  // for a given mass, calculate the corresponding th ckov value
  // 
  float calcCkovFromMass(float p, float n, int pdg) {
    // Constants for particle masses (in GeV/c^2)
    const float mass_Muon = 0.10566, mass_Pion = 0.1396, mass_Kaon = 0.4937,
                mass_Proton = 0.938;
    auto particlePDG = TDatabasePDG::Instance()->GetParticle(pdg);
    double mass = particlePDG ? particlePDG->Mass() : 0.;

    float m;          // mass
    p = std::abs(p);  // momentum


    // if passed pdg is neither of 13, 211, 321 or 2212, the function returns 0
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

  bool findPhotCkov(double cluX, double cluY, double &thetaCer, double &phiCer) {
    // Finds Cerenkov angle  for this photon candidate
    // Arguments: cluX,cluY - position of cadidate's cluster
    // Returns: Cerenkov angle

    TVector3 dirCkov;

    double zRad =
        -0.5 * radThick() - 0.5 * winThick(); // z position of middle of RAD
    TVector3 rad(fTrkPos.X(), fTrkPos.Y(),
                 zRad); // impact point at middle of RAD
    TVector3 pc(cluX, cluY, 0.5 * winThick() + gapThick()); // mip at PC
    double cluR = TMath::Sqrt(
        (cluX - fMipPos.X()) * (cluX - fMipPos.X()) +
        (cluY - fMipPos.Y()) *
            (cluY - fMipPos.Y())); // ref. distance impact RAD-CLUSTER

    double phi = (pc - rad).Phi(); // phi of photon

    double ckov1 = 0;
    double ckov2 =
        0.75 + fTrkDir.Theta(); // start to find theta cerenkov in DRS
    const double kTol = 0.01;
    Int_t iIterCnt = 0;

    while (1) {

      if (iIterCnt >= 50) {

        return kFALSE;
      }

      double ckov = 0.5 * (ckov1 + ckov2);
      dirCkov.SetMagThetaPhi(1, ckov, phi);
      auto posC = traceForward(dirCkov); // trace photon with actual angles

      double dist =
          cluR -
          (posC - fMipPos)
              .Mod(); // get distance between trial point and cluster position

      if (posC.X() == -999) {
        dist = -999; // total reflection problem
      }

      iIterCnt++; // counter step

      if (dist > kTol) {
        ckov1 = ckov; // cluster @ larger ckov
      } else if (dist < -kTol) {
        ckov2 = ckov; // cluster @ smaller ckov
      } else { // precision achived: ckov in DRS found
        dirCkov.SetMagThetaPhi(1, ckov, phi);
        lors2Trs(dirCkov, thetaCer,
                 phiCer); // find ckov (in TRS:the effective Cherenkov angle!)
        return kTRUE;
      }

      // std::cout << "ckovTh , phi" << ckovTh << " "<< phi << "  dist s "<<
      // dist << std::endl;
    }
  } // FindPhotTheta()

  MCRecon(double xRad, double yRad,
        const std::vector<o2::hmpid::Cluster> &clustersIn,
        const o2::dataformats::MatchInfoHMP &matchInfo) {

    numCkovHough = matchInfo.getNPhots();
    numCkovHoughMH = matchInfo.getMassHypNumPhot();

    std::cout << "Recon :: Number of Cerenkov photons detected: "
              << numCkovHough << std::endl;
    std::cout << "Recon :: Number of Cerenkov photons under mass hypothesis: "
              << numCkovHoughMH << std::endl;

    mSizeMip = matchInfo.getMipClusSize();

    mNphots, mQMip;
    matchInfo.getHMPIDmip(mxMip, myMip, mQMip, mNphots);
    matchInfo.getHMPIDtrk(mXpc, mYpc, mThetaP, mPhiP);

    fTrkDir.SetMagThetaPhi(1., mThetaP, mPhiP);
    fTrkPos.Set(xRad, yRad);
    fMipPos.Set(mxMip, myMip);

    fPc.Set(mXpc, mYpc);
    mClusters = clustersIn;
    mMatchInfo = matchInfo;

    mxRa = xRad;
    myRa = yRad;


    // momentum
    mP = std::abs(mMatchInfo.getHmpMom());


    // theoretical ckov angles for the track
    pionCkovTh = calcCkovFromMass(mP, refIndexTmp, 211);
    kaonCkovTh = calcCkovFromMass(mP, refIndexTmp, 321);
    protonCkovTh = calcCkovFromMass(mP, refIndexTmp, 2212);


    // reconstructed ckov values using standard HTM and HTM with masshyp
    ckovRecon = matchInfo.getHMPsignal();
    ckovReconMassHyp = matchInfo.getHMPsignalMassHyp();

    auto nClusters = clustersIn.size();



    // reserving space for the vectors 
    pNormPion.resize(nClusters, 0.0f), pNormKaon.resize(nClusters, 0.0f),
        pNormProton.resize(nClusters, 0.0f);

    sumProbTrack.resize(nClusters, 0.0f);

    pProton.resize(nClusters, 0.0f);
    pKaon.resize(nClusters, 0.0f);
    pPion.resize(nClusters, 0.0f);

    xValues.resize(nClusters, 0.0f);
    yValues.resize(nClusters, 0.0f);
    qValues.resize(nClusters, 0.0f);
    thetaCerValues.resize(nClusters, 0.0f);
    phiCerValues.resize(nClusters, 0.0f);
    sigmaRingValues.resize(nClusters, 0.0f);
    sizeValues.resize(nClusters, 0);

    rawSizeValues.resize(nClusters, 0.0f);
    numRawClustersValues.resize(nClusters, 0.0f);

    zProton.resize(nClusters, 0.0f);
    zKaon.resize(nClusters, 0.0f);
    zPion.resize(nClusters, 0.0f);
  }

  // double thetaP() const { return mThetaP; }
  // double phiP() const { return mPhiP; }

  // ef : TODO : why do we not take into account the spread in ref-index here?
  double getThetaCerPhotThresh() 
  {
    return 1/getRefIndex();
  }


  void process(int pdgCodeTrack,
               const std::vector<int> &indexOfCkovClusters) {

    // to get thetaCer, phiCer...

    // phi_cer_values[i] = -13.;
    // theta_cer_values[i] = -13.;
    // sigma_ring_values[i] = -13.;

    // bool findPhotCkov(double cluX, double cluY, double& thetaCer, double&
    // phiCer)k
    const float ckovRecon = mMatchInfo.getHMPsignal();
    const int vectorSize = mClusters.size();

    int cnt = 0;

    std::vector<float> phiProton, phiKaon, phiPion;

    auto diffPion = std::abs(getCkovThPion() - ckovRecon);
    auto diffKaon = std::abs(getCkovThKaon() - ckovRecon);
    auto diffProton = std::abs(getCkovThProton() - ckovRecon);



    /* ef > verbose prining
    Printf("process : mP %.2f getThetaP %.2f getPhiP %.2f ", mP, getThetaP(),
           getPhiP());

    Printf("process : ckovRecon %.4f theoPi %.4f theoKa %.4f theoPro %.4f ",
           ckovRecon, getCkovThPion(), getCkovThKaon(), getCkovThProton());

    Printf("process : DIFFS theoPi %.4f theoKa %.4f theoPro %.4f ", diffPion,
           diffKaon, diffProton);

    if (diffPion <= diffKaon && diffPion <= diffProton) {
      std::cout << "Least difference is for Pion: " << std::endl;
    } else if (diffKaon <= diffPion && diffKaon <= diffProton) {
      std::cout << "Least difference is for Kaon: " << std::endl;
    } else {
      std::cout << "Least difference is for Proton: " << std::endl;
    }
    */

    // all values initialzied to zero > only set them if they are good
    // candidates of one of the species

    double sumCkovWindow = 0;
    int nCkovPhots = 0;

    int indexOfclu = 0;

    for (const auto &cluster : mClusters) {
      double x = cluster.x();
      double y = cluster.y();
      double q = cluster.q();
      double size = cluster.size();
      double thetaCer, phiCer, sigmaRing;

      Printf("\n");
      for (const auto &idx : indexOfCkovClusters) {
        if (idx == indexOfclu) {
          LOGP(info, "this is a ckov cluster!");
        }
      }
      indexOfclu++;

      // only consider clusters that show traits
      // of being single-photon clusters
      if (q < mThreshMipQ || size < mThreshMipSize) { // attributes of MIP is the opposite

        // thetaCer, phiCer passed by reference
        if (this->findPhotCkov(x, y, thetaCer, phiCer)) {


          // if thetaCer is over threshold, do not consider photon
          if (thetaCer > 1 / getThetaCerPhotThresh())
            continue;

          // naive average ckov, only used for debugging
          sumCkovWindow += thetaCer;
          
          
          nCkovPhots++;

          o2::hmpid::Param *p = o2::hmpid::Param::instanceNoGeo();


          // angular resolution for the cluster
          sigmaRing = std::sqrt(p->sigma2(getThetaP(), getPhiP(), thetaCer, phiCer));

          auto dist2mip = std::sqrt((getMipX() - x) * (getMipX() - x) +
                                    (getMipY() - y) * (getMipY() - y));


          double sigmaRingCap = sigmaRing;
          // saturate the angular resolution >
          // only for calculating z-score and prob
          if (sigmaRing > 0.035)
            sigmaRingCap = 0.035;


          // calculate teh z score per specie 
          zProton[cnt] = (thetaCer - protonCkovTh) / sigmaRingCap;
          zKaon[cnt] = (thetaCer - kaonCkovTh) / sigmaRingCap;
          zPion[cnt] = (thetaCer - pionCkovTh) / sigmaRingCap;

          Printf(" ZPi %.3f ZKa %.3f ZPi %.3f", zPion[cnt], zKaon[cnt],
                 zProton[cnt]);

          // to evaluate the number of photon-candidates
          // to filter out tracks with too few candidates

          // phiProton > used to evaluate if the number of photons is within defined expected range
          // check for all the phi region, and also in the half-plane

          // pProton > will remain zero for this index if not within 2 std-devs

          if (std::abs(zProton[cnt]) < 2) {
            phiProton.push_back(phiCer);
            pProton[cnt] = zScoreToProb(zProton[cnt]);
          }

          if (std::abs(zKaon[cnt]) < 2) {
            phiKaon.push_back(phiCer);
            pKaon[cnt] = zScoreToProb(zKaon[cnt]);
          }

          if (std::abs(zPion[cnt]) < 2) {
            phiPion.push_back(phiCer);
            pPion[cnt] = zScoreToProb(zPion[cnt]);
          }

          bool isCandOfAnyHadron =
              ((std::abs(zProton[cnt]) < 2) || (std::abs(zKaon[cnt]) < 2) ||
               (std::abs(zPion[cnt]) < 2));

          // if (isCandOfAnyHadron) {
          xValues[cnt] = x;
          yValues[cnt] = y;
          qValues[cnt] = q;
          thetaCerValues[cnt] = thetaCer;
          phiCerValues[cnt] = phiCer;
          sigmaRingValues[cnt] = sigmaRing;
          sizeValues[cnt] = size;

          numRawClustersValues[cnt] = cluster.numLocMax();
          rawSizeValues[cnt] = cluster.sizeRaw();
          //}
        }
      }
      cnt++;
    }

    std::cout << std::endl;
    if (nCkovPhots > 1) {
      LOGP(info, "naive average ckov Th = {}", sumCkovWindow / nCkovPhots);
    }

    FilterPhotons filterPhotons;

    LOGP(info, "Thresholds:");

    for (const auto &pdg : pdgs) {
      std::pair<int, int> thresholds =
          filterPhotons.expectedNumPhotons(refIndexTmp, mP, pdg);
      LOGP(info, "PDG {} >>  full {} half {}", pdg, thresholds.first,
           thresholds.second);
    }

    mipPcDists = std::sqrt((mxMip - mXpc) * (mxMip - mXpc) +
                           (myMip - mYpc) * (myMip - mYpc));

    // only consider track if it has the following
    // 1) From MatchHMP > distCut 6, mipCharge and mipSize
    // 2) from FilterPhotons : one of the ckov zones of the species fulffils the
    // requierement

    bool fixedPdgIsNumPhotsOk = false;
    if (filterPhotons.evaluateNumPhotThresh(phiProton, phiKaon, phiPion,
                                            refIndexTmp, mP, pdgCodeTrack,
                                            fixedPdgIsNumPhotsOk)) {
      mIsTrackMatched = true;

      // tighter cut for tracks to actually be reconstructed
      if (mipPcDists < distThresh) {
        mIsTrackToReconstruct = true;

        // ef > TODO add the ckov ring average ang resolution to the
        // momentum-threshold_?


        //  use this as a requirenemnt on the numPhots requirement?
        // might biase the data, as we here imply the PDG reading from the MC truth
        if (fixedPdgIsNumPhotsOk) {
          setTrackReconKnownPdg(fixedPdgIsNumPhotsOk);
        }
      }

      LOGP(info, "had enough photons > pr {} ka {} pi {}", phiProton.size(),
           phiKaon.size(), phiPion.size());
    } else {
      LOGP(info, "not enough photons > pr {} ka {} pi {}", phiProton.size(),
           phiKaon.size(), phiPion.size());
    }

  } // end process

  float getMipX() const { return mxMip; }

  float getMipY() const { return myMip; }

  float getRadX() const { return mxRa; }

  float getRadY() const { return myRa; }

  float getXpc() const { return mXpc; }

  float getYpc() const { return mYpc; }

  float getThetaP() const { return mThetaP; }

  float getPhiP() const { return mPhiP; }

  float getMomentum() const { return mP; }

  int getQMip() const { return mQMip; }

  int getNumPhots() const { return mNphots; }

  int getSizeMip() const { return mSizeMip; }

  float getMipPcDist() const { return mipPcDists; }

  // Theoretical Cherenkov values for species for the given track
  float getCkovThPion() const { return pionCkovTh; }
  float getCkovThKaon() const { return kaonCkovTh; }
  float getCkovThProton() const { return protonCkovTh; }

  float getRefIndex() const { return refIndexTmp; }


  // the ckov reconstructd using the classic method
  float getCkovRecon() const { return ckovRecon; }

  // the ckov reconstructd using the Mass hypothesis
  float getCkovReconMassHyp() const { return ckovReconMassHyp; }

  // the number of selected Hough photons using the classic method
  int getNumCkovHough() const { return numCkovHough; }

  // the number of selected Hough photons using the Mass hypothesis
  int getNumCkovHoughMH() const { return numCkovHoughMH; }

private:
  static constexpr auto pdgs = {211, 321, 2212};
  std::vector<float> mCkovTheoretical = {-1., -1., -1.};

  // cluster attributes
  std::vector<int> sizeValues;
  std::vector<float> xValues, yValues, qValues, thetaCerValues, phiCerValues,
      sigmaRingValues;


  std::vector<int> numRawClustersValues;// number of raw clusters that forms the convoluted clusters
  std::vector<int> rawSizeValues;       // size of the covoluted cluster
  

  // z-szore and prob for all teh clusters belonging to a specific track
  std::vector<float> zProton, zKaon, zPion;
  std::vector<float> pProton, pKaon, pPion;


  // flag which states wether we should reconstruct the given track
  // needs to be matched with a MIP (only meeting distance and Charge requirement, not being matched to the correct MIP)
  bool mIsTrackToReconstruct = false;

  // ef > if we take into account that particle specie is known
  bool mIsTrackToReconstructWithKnownPdg = false;


  // track attributes
  bool mIsTrackMatched = false;
  float mXpc, mYpc;
  const float refIndexTmp = 1.2904;
  float distThresh = 2.0f;
  float mxMip, myMip;
  float mxRa, myRa;

  float protonCkovTh, kaonCkovTh, pionCkovTh;

  float ckovReconMassHyp, ckovRecon;

  int mQMip;    // MIP charge ADC
  int mSizeMip; // MIP size
  int mNphots;  // num photons in matched track-MIP

  float mP;     // track momentum

  float mipPcDists;

  std::vector<o2::hmpid::Cluster> mClusters = {};

  o2::dataformats::MatchInfoHMP mMatchInfo;

  TVector3 fTrkDir; // track direction in LORS at RAD

  TVector2 fTrkPos; // track positon in LORS at RAD   // XY mag
  TVector2 fMipPos; // mip positon for a given trackf // XY
  TVector2 fPc;     // track position at PC           // XY
  float mThetaP, mPhiP;

  std::vector<float> sumProbTrack;
  std::vector<float> pNormPion, pNormKaon, pNormProton;

  int numCkovHoughMH, numCkovHough;

  McTruth mMcTruth;

  double radThick() const {
    return 1.5;
  } //<--TEMPORAR--> to be removed in future. Radiator thickness
  double winThick() const {
    return 0.5;
  } //<--TEMPORAR--> to be removed in future. Window thickness
  double gapThick() const {
    return 8.0;
  } //<--TEMPORAR--> to be removed in future. Proximity gap thickness
  double winIdx() const {
    return 1.583;
  } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN
    // material (SiO2)
  double gapIdx() const {
    return 1.0005;
  } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP
    // material (CH4)
  double getRefIdx() const {
    return 1.2904;
  } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP
    // material (CH4)

  /*
  def winIdx(self):
      return np.sqrt(1 + 46.411 / (10.666 * 10.666 - self.eV * self.eV) + 228.71
  / (18.125 * 18.125 - self.eV * self.eV))

  def gapIdx(self):
      return 1 + 0.12489e-6 / (2.62e-4 - self.eV * self.eV / 1239.84 / 1239.84)

  def nIdxRad(self, eV, temp=25):  # Default temperature is 20 unless provided
  otherwise eV_term = 1239.84 / self.eV return np.sqrt(1 + 0.554 * eV_term *
  eV_term / (eV_term * eV_term - 5769)) - 0.0005 * (temp - 20) */

  void propagate(const TVector3 dir, TVector3 &pos, double z) const {
    // Finds an intersection point between a line and XY plane shifted along Z.
    // Arguments:  dir,pos   - vector along the line and any point of the line
    //             z         - z coordinate of plain
    //   Returns:  none
    //   On exit:  pos is the position if this intesection if any
    static TVector3 nrm(0, 0, 1);
    TVector3 pnt(0, 0, z);

    TVector3 diff = pnt - pos;
    double sint = (nrm * diff) / (nrm * dir);
    pos += sint * dir;
  } // Propagate()
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void refract(TVector3 &dir, double n1, double n2) const {
    // Refract direction vector according to Snell law
    // Arguments:
    //            n1 - ref idx of first substance
    //            n2 - ref idx of second substance
    //   Returns: none
    //   On exit: dir is new direction

    double sinref = (n1 / n2) * TMath::Sin(dir.Theta());
    if (TMath::Abs(sinref) > 1.) {
      dir.SetXYZ(-999, -999, -999);
    } else {
      dir.SetTheta(TMath::ASin(sinref));
    }
  }

  TVector2 traceForward(TVector3 dirCkov) const {
    // Trace forward a photon from (x,y) up to PC
    //  Arguments: dirCkov photon vector in LORS
    //    Returns: pos of traced photon at PC

    TVector2 pos(-999, -999);
    double thetaCer = dirCkov.Theta();
    if (thetaCer > TMath::ASin(1. / getRefIdx())) {
      return pos;
    } // total refraction on WIN-GAP boundary
    double zRad =
        -0.5 * radThick() - 0.5 * winThick(); // z position of middle of RAD
    TVector3 posCkov(
        fTrkPos.X(), fTrkPos.Y(),
        zRad); // RAD: photon position is track position @ middle of RAD

    propagate(dirCkov, posCkov, -0.5 * winThick()); // go to RAD-WIN boundary

    refract(dirCkov, getRefIdx(), winIdx()); // RAD-WIN refraction

    propagate(dirCkov, posCkov, 0.5 * winThick()); // go to WIN-GAP boundary

    refract(dirCkov, winIdx(), gapIdx()); // WIN-GAP refraction

    propagate(dirCkov, posCkov, 0.5 * winThick() + gapThick()); // go to PC

    pos.Set(posCkov.X(), posCkov.Y());
    return pos;

  } // TraceForward()

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void lors2Trs(TVector3 dirCkov, double &thetaCer, double &phiCer) const {
    // Theta Cerenkov reconstruction
    //  Arguments: dirCkov photon vector in LORS
    //    Returns: thetaCer of photon in TRS
    //               phiCer of photon in TRS
    //  TVector3 dirTrk;
    //  dirTrk.SetMagThetaPhi(1,fTrkDir.Theta(),fTrkDir.Phi()); ->
    //  dirTrk.SetCoordinates(1,fTrkDir.Theta(),fTrkDir.Phi()) double  thetaCer
    //  = TMath::ACos(dirCkov*dirTrk);

    TRotation mtheta;
    mtheta.RotateY(-fTrkDir.Theta());

    TRotation mphi;
    mphi.RotateZ(-fTrkDir.Phi());

    TRotation mrot = mtheta * mphi;

    TVector3 dirCkovTRS;
    dirCkovTRS = mrot * dirCkov;
    phiCer = dirCkovTRS.Phi(); // actual value of the phi of the photon
    thetaCer =
        dirCkovTRS.Theta(); // actual value of thetaCerenkov of the photon
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  void trs2Lors(TVector3 dirCkov, double &thetaCer, double &phiCer) const {
    // Theta Cerenkov reconstruction
    //  Arguments: dirCkov photon vector in TRS
    //    Returns: thetaCer of photon in LORS
    //               phiCer of photon in LORS

    // TRotation mtheta;
    // mtheta.RotateY(fTrkDir.Theta()); ef : changed to :

    TRotation mtheta;
    mtheta.RotateY(fTrkDir.Theta());

    TRotation mphi;
    mphi.RotateZ(fTrkDir.Phi());

    TRotation mrot = mphi * mtheta;

    TVector3 dirCkovLORS;
    dirCkovLORS = mrot * dirCkov;

    phiCer = dirCkovLORS.Phi(); // actual value of the phi of the photon
    thetaCer =
        dirCkovLORS.Theta(); // actual value of thetaCerenkov of the photon
  }


};

#endif
