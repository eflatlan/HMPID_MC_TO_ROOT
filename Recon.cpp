#include "SimpleCluster.cpp"
#include "Alisigma2_.cpp"
#include "FilterPhotons.cpp"
#include "ReconstructionDataFormats/MatchInfoHMP.h"
#include <boost/math/distributions/normal.hpp>
#include "DataFormatsHMP/Cluster.h"

#include "HMPIDBase/Param.h"


#ifndef Recon_H
#define Recon_H


struct McTruth {
  int pdgCodeTrack = -1;
  int pdgCodeClu = -1;
  bool isMipMatchedCorrectly = false;  
  bool isTrackToReconKnownPdg = false; // if we take pdg code being known when making the requirements
 
  int numCountedPhotonsPC = -1; // numbef of photon-clusters on PC having the track as MID 

  McTruth(bool isMipMatched, int pdgCode, int mipPDG, bool isTrackToRecon, int nPc) :  isMipMatchedCorrectly(isMipMatched), pdgCodeClu(mipPDG), pdgCodeTrack(pdgCode), isTrackToReconKnownPdg(isTrackToRecon), numCountedPhotonsPC(nPc)
  {}

  McTruth() {}

};


class Recon {

public:


  // ef > if we take into account PDG code is known
  bool isTrackToReconKnownPdg()
  {
    return mIsTrackToReconstructWithKnownPdg;
  }

    // ef > if we take into account PDG code is known
  bool isTrackToReconKnown()
  {
    return mIsTrackToReconstruct;
  }


  void setTrackReconKnownPdg(bool isTrackToReconstructWithKnownPdg)
  {
    mIsTrackToReconstructWithKnownPdg = isTrackToReconstructWithKnownPdg;
  }


  void calculateNormalizedProb(std::vector<float>& sumProbAllTracks)  
  {


    /*
    LOGP(info, "sizes of vectors >> pPion {} pKaon {} pProton {} sumProbTrack {}", pPion.size(), pKaon.size(), pProton.size(), sumProbTrack.size() );
 

    std::cout << "\n sumProbTrackIn : ";
    for(size_t i = 0; i < sumProbTrack.size(); ++i) {
        if(sumProbTrack[i] != 0.0) { 
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2) << sumProbTrack[i] << " "; 
        }
    }*/

    // add all fields 
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pProton.begin(), sumProbTrack.begin(), std::plus<>());
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pKaon.begin(), sumProbTrack.begin(), std::plus<>());
    std::transform(sumProbTrack.begin(), sumProbTrack.end(), pPion.begin(), sumProbTrack.begin(), std::plus<>());

    // Normalize pion probabilities
    std::transform(pPion.begin(), pPion.end(), sumProbTrack.begin(), pNormPion.begin(),
                   [](float pionProb, float sumProb) {
                       return sumProb > 0 ? pionProb / sumProb : 0.0f; // Avoid division by zero
                   });

    std::transform(pKaon.begin(), pKaon.end(), sumProbTrack.begin(), pNormKaon.begin(),
                   [](float kaonProb, float sumProb) {
                       return sumProb > 0 ? kaonProb / sumProb : 0.0f; // Avoid division by zero
                   });

    std::transform(pProton.begin(), pProton.end(), sumProbTrack.begin(), pNormProton.begin(),
                   [](float protonProb, float sumProb) {
                       return sumProb > 0 ? protonProb / sumProb : 0.0f; // Avoid division by zero
                   });


    std::transform(sumProbAllTracks.begin(), sumProbAllTracks.end(), sumProbTrack.begin(), sumProbAllTracks.begin(), std::plus<>());


    /*
    std::cout << "\npPion  >>>";
    for(size_t i = 0; i < pPion.size(); ++i) {
        if(pPion[i] != 0.0) { 
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2) << pPion[i] << " "; 
        }
    }
    std::cout << "\npKaon  >>>";
    for(size_t i = 0; i < pKaon.size(); ++i) {
        if(pKaon[i] != 0.0) { 
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2) << pKaon[i] << " "; 
        }
    }


    std::cout << "\npProt  >>>";
    for(size_t i = 0; i < pProton.size(); ++i) {
        if(pProton[i] != 0.0) { 
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2) << pProton[i] << " "; 
        }
    }


    std::cout << "\nsumProbTrackOut  >>>";

    for(size_t i = 0; i < sumProbTrack.size(); ++i) {
        if(sumProbTrack[i] != 0.0) { 
            std::cout << "i: " << i << " " << std::fixed << std::setprecision(2) << sumProbTrack[i] << " "; 
        }
    }
    std::cout<<" \n";*/

  }


  void setMcTruth(const McTruth& mcTruth)  
  {
    mMcTruth = mcTruth;
  }

  void getMcTruth(McTruth& mcTruth) const 
  {
    mcTruth = mMcTruth;
  }



  bool isTrackMatched() const 
  {
    return mIsTrackMatched;
  }


  bool isTrackAppropriate() const 
  {
    return mIsTrackToReconstruct;
  }

  std::vector<float> getSumProbTrack() const 
  {
    return sumProbTrack;
  }

  std::vector<float> getProtonProbNorm() const 
  {
    return pNormProton;
  }

  std::vector<float> getKaonProbNorm() const 
  {
    return pNormKaon;
  }

  std::vector<float> getPionProbNorm() const 
  {
    return pNormPion;
  }

  std::vector<float> getXValues() const {
      return xValues;
  }

  std::vector<float> getYValues() const {
      return yValues;
  }

  std::vector<float> getQValues() const {
      return qValues;
  }

  std::vector<float> getThetaCerValues() const {
      return thetaCerValues;
  }

  std::vector<float> getPhiCerValues() const {
      return phiCerValues;
  }

  std::vector<float> getSigmaRingValues() const {
      return sigmaRingValues;
  }

  std::vector<int> getSizeValues() const {
      return sizeValues;
  }

  std::vector<int> getNumClustersValues() const {
    return numRawClustersValues;
  }

  std::vector<int> getRawSizeValues() const {
    return rawSizeValues;
  }


  std::vector<float> getProtonProb() const 
  {
    return pProton;
  }

  std::vector<float> getKaonProb() const 
  {
    return pKaon;
  }

  std::vector<float> getPionProb() const 
  {
    return pPion;
  }

  double zScoreToProb(double z) {
      boost::math::normal_distribution<> dist(0.0, 1.0); // Mean = 0, Standard Deviation = 1 for a standard normal distribution
      return boost::math::cdf(dist, z);                  // Cumulative Distribution Function value for z
  }


  // cyclic dependance, can be resolved ofc
  /*void addTracks(Tracks::Tracks& tracks) 
  {
    tracks.addTrackCands(this);
  }*/ 


  float calcMassFromCkov(float p, float n, float ckov)
  {
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

  float calcCkovFromMass(float p, float n, int pdg)
  {
    // Constants for particle masses (in GeV/c^2)
    const float mass_Muon = 0.10566, mass_Pion = 0.1396, mass_Kaon = 0.4937, mass_Proton = 0.938;
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



  bool findPhotCkov(double cluX, double cluY, double &thetaCer,
                    double &phiCer) {
    // Finds Cerenkov angle  for this photon candidate
    // Arguments: cluX,cluY - position of cadidate's cluster
    // Returns: Cerenkov angle

    TVector3 dirCkov;

    double zRad =
        -0.5 * radThick() - 0.5 * winThick(); // z position of middle of RAD
    TVector3 rad(fTrkPos.X(), fTrkPos.Y(),
                 zRad); // impact point at middle of RAD
    TVector3 pc(cluX, cluY, 0.5 * winThick() + gapThick()); // mip at PC
    double cluR =
        TMath::Sqrt((cluX - fMipPos.X()) * (cluX - fMipPos.X()) +
                    (cluY - fMipPos.Y()) *
                        (cluY - fMipPos.Y())); // ref. distance impact RAD-CLUSTER



    double phi = (pc - rad).Phi();         // phi of photon

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

  Recon(double xRad, double yRad, const std::vector<o2::hmpid::Cluster> &clustersIn, const o2::dataformats::MatchInfoHMP &matchInfo) {


    numCkovHough = matchInfo.getNPhots();
    numCkovHoughMH = matchInfo.getMassHypNumPhot();

    std::cout << "Recon :: Number of Cerenkov photons detected: " << numCkovHough << std::endl;
    std::cout << "Recon :: Number of Cerenkov photons under mass hypothesis: " << numCkovHoughMH << std::endl;


    // ef : this was added now
    mSizeMip = matchInfo.getMipClusSize();



    mNphots, mQMip;
    matchInfo.getHMPIDmip(mxMip, myMip, mQMip, mNphots);
    matchInfo.getHMPIDtrk(mXpc, mYpc, mThetaP, mPhiP);


    fTrkDir.SetMagThetaPhi(1., mThetaP, mPhiP);
    fTrkPos.Set(xRad, yRad);
    fMipPos.Set(mxMip, myMip);

    // should this be xPc, yPc ? 
    fPc.Set(mXpc, mYpc);
    mClusters = clustersIn;
    mMatchInfo = matchInfo;

    mxRa = xRad;
    myRa = yRad;

    mP =  std::abs(mMatchInfo.getHmpMom());


    pionCkovTh = calcCkovFromMass(mP, refIndexTmp, 211);                        
    kaonCkovTh = calcCkovFromMass(mP, refIndexTmp, 321);                        
    protonCkovTh = calcCkovFromMass(mP, refIndexTmp, 2212);                        

    ckovRecon = matchInfo.getHMPsignal();
    ckovReconMassHyp = matchInfo.getHMPsignalMassHyp();

    auto nClusters = clustersIn.size();

    pNormPion.resize(nClusters, 0.0f), pNormKaon.resize(nClusters, 0.0f), pNormProton.resize(nClusters, 0.0f);
    
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

  double thetaP() const {
    return mThetaP;
  } 
  double phiP() const {
    return mPhiP;
  } 




  void process(const Alisigma2_& aliSigma2, int pdgCodeTrack, const std::vector<int>& indexOfCkovClusters) {

    // to get thetaCer, phiCer...

    //phi_cer_values[i] = -13.;
    //theta_cer_values[i] = -13.;
    //sigma_ring_values[i] = -13.;

    // bool findPhotCkov(double cluX, double cluY, double& thetaCer, double&
    // phiCer)k
    const float ckovRecon = mMatchInfo.getHMPsignal();
    const int vectorSize = mClusters.size();


    int cnt = 0;

    std::vector<float> phiProton, phiKaon, phiPion;

    auto diffPion = std::abs(getCkovThPion()-ckovRecon);
    auto diffKaon = std::abs(getCkovThKaon()-ckovRecon);
    auto diffProton = std::abs(getCkovThProton()-ckovRecon);
    
    Printf("process : mP %.2f getThetaP %.2f getPhiP %.2f ",  mP, getThetaP(), getPhiP());

    Printf("process : ckovRecon %.4f theoPi %.4f theoKa %.4f theoPro %.4f ",  ckovRecon, getCkovThPion(), getCkovThKaon(), getCkovThProton());

    Printf("process : DIFFS theoPi %.4f theoKa %.4f theoPro %.4f ",  diffPion, diffKaon, diffProton);

    if (diffPion <= diffKaon && diffPion <= diffProton) {
        std::cout << "Least difference is for Pion: " << std::endl;
    } else if (diffKaon <= diffPion && diffKaon <= diffProton) {
        std::cout << "Least difference is for Kaon: " << std::endl;
    } else {
        std::cout << "Least difference is for Proton: " << std::endl;
    }



    // ef > can we do this when we know the specie-type?
    




    // all values initialzied to zero > only set them if they are good candidates 
    // of one of the species



    double sumCkovWindow = 0;
    int nCkovPhots = 0;

    int indexOfclu = 0;

    for(const auto& cluster : mClusters) {
      double x = cluster.x();
      double y = cluster.y();
      double q = cluster.q();
      double size = cluster.size();
      double thetaCer, phiCer, sigmaRing;

      Printf("\n");
      for(const auto& idx : indexOfCkovClusters) {
        if(idx == indexOfclu) {
          LOGP(info, "this is a ckov cluster!");
        }
      }
      indexOfclu++;

      // only consider clusters that show traits 
      // of being single-photon clusters
      if(q < 200 || size < 3) { // attributes of MIP is the opposite

        if (this->findPhotCkov(x, y, thetaCer, phiCer)) {
          
          if(thetaCer > 1/1.2)
            continue;


          sumCkovWindow+=thetaCer;
          nCkovPhots++;
          sigmaRing = aliSigma2.sigma2(thetaP(), phiP(), thetaCer, phiCer);



          //Try using fParam here >
          //
          o2::hmpid::Param* p = o2::hmpid::Param::instanceNoGeo();
          auto sigmaRing2 = aliSigma2.sigma2(thetaP(), phiP(), thetaCer, phiCer);
          sigmaRing = std::sqrt(p->sigma2(thetaP(), phiP(), thetaCer, phiCer));

          auto dist2mip = std::sqrt((getMipX()-x)*(getMipX()-x) +(getMipY()-y)*(getMipY()-y));
          Printf("dist2mip %.2f",  dist2mip);

          Printf(" sigmaRing %.4f  sigmaRing2 %.4f", sigmaRing2, sigmaRing);


          
          
          Printf("xMip %.4f x %.4f  yMip %.4f y %.4f",  getMipX(), x, getMipY(), y);

          Printf("findPhotCkov : mP %.2f ckovRecon %.2f thetaCer %.2f phiCer %.2f  sigmaRing %.4f",  mP, ckovRecon, thetaCer, phiCer, sigmaRing);


          double sigmaRingCap = sigmaRing;
          // saturate the angular resolution > 
          // only for calculating z-score and prob
          if(sigmaRing > 0.035)
            sigmaRingCap = 0.035;

          zProton[cnt] = (thetaCer-protonCkovTh)/sigmaRingCap;
          zKaon[cnt] = (thetaCer-kaonCkovTh)/sigmaRingCap;
          zPion[cnt] = (thetaCer-pionCkovTh)/sigmaRingCap;


          Printf(" ZPi %.3f ZKa %.3f ZPi %.3f", zPion[cnt], zKaon[cnt], zProton[cnt]);


          // to evaluate the number of photon-candidates
          // to filter out tracks with too few candidates
          if(std::abs(zProton[cnt]) < 2) {
            phiProton.push_back(phiCer);
            pProton[cnt] = zScoreToProb(zProton[cnt]);
          }

          if(std::abs(zKaon[cnt]) < 2) {
            phiKaon.push_back(phiCer);
            pKaon[cnt] = zScoreToProb(zKaon[cnt]);
          }
          
          if(std::abs(zPion[cnt]) < 2) {
            phiPion.push_back(phiCer);
            pPion[cnt] = zScoreToProb(zPion[cnt]);            
          }

          bool isCandOfAnyHadron = ((std::abs(zProton[cnt]) < 2) || (std::abs(zKaon[cnt]) < 2) || (std::abs(zPion[cnt]) < 2));

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

    std::cout<<std::endl;
    if (nCkovPhots > 1) {
      LOGP(info, "naive average ckov Th = {}", sumCkovWindow/nCkovPhots);
    }

    FilterPhotons filterPhotons;


    LOGP(info, "Thresholds:");

    for(const auto& pdg : pdgs)
    {
      std::pair<int, int> thresholds = filterPhotons.expectedNumPhotons(refIndexTmp, mP, pdg);
      LOGP(info, "PDG {} >>  full {} half {}", pdg, thresholds.first, thresholds.second);
    }


    mipPcDists = std::sqrt((mxMip-mXpc)*(mxMip-mXpc)+(myMip-mYpc)*(myMip-mYpc));

    // only consider track if it has the following 
    // 1) From MatchHMP > distCut 6, mipCharge and mipSize
    // 2) from FilterPhotons : one of the ckov zones of the species fulffils the requierement

    
    bool fixedPdgIsNumPhotsOk = false;
    if(filterPhotons.evaluateNumPhotThresh(phiProton, phiKaon, phiPion, refIndexTmp, mP, pdgCodeTrack, fixedPdgIsNumPhotsOk))
    {
      mIsTrackMatched = true;

      // tighter cut for tracks to actually be reconstructed
      if(mipPcDists < distThresh) {
        mIsTrackToReconstruct = true;

        // and use this as a requirenemnt on the numPhots requirement?
        // ef > TODO add the ckov ring average ang resolution to the momentum-threshold_?
        if(fixedPdgIsNumPhotsOk) {
          setTrackReconKnownPdg(fixedPdgIsNumPhotsOk);
        } 

      }

      LOGP(info, "had enough photons > pr {} ka {} pi {}", phiProton.size(),  phiKaon.size(),  phiPion.size());
    } else {
      LOGP(info, "not enough photons > pr {} ka {} pi {}", phiProton.size(),  phiKaon.size(),  phiPion.size());
    }

  } // end process


  float getMipX() const {
      return mxMip;
  }

  float getMipY() const {
      return myMip;
  }

  float getRadX() const {

      return mxRa;
  }

  float getRadY() const {
      return myRa;
  }

  float getXpc() const {
      return mXpc;
  }

  float getYpc() const {
      return mYpc;
  }

  float getThetaP() const {
      return mThetaP;
  }

  float getPhiP() const {
      return mPhiP;
  }

  float getMomentum() const {
      return mP;
  }

  int getQMip() const {
      return mQMip;
  }

  int getNumPhots() const {
      return mNphots;
  }

  int getSizeMip() const {
      return mSizeMip;
  }  

  float getMipPcDist() const {
      return mipPcDists;
  }

  float getCkovThPion() const {
      return pionCkovTh;
  }

  float getCkovThKaon() const {
      return kaonCkovTh;
  }

  float getCkovThProton() const {
      return protonCkovTh;
  }

  float getRefIndex() const {
    return refIndexTmp;
  }


  float getCkovRecon() const {
      return ckovRecon;
  }

  float getCkovReconMassHyp() const {
    return ckovReconMassHyp;
  }


  int getNumCkovHough() const {
    return numCkovHough;
  }

  int getNumCkovHoughMH() const {
    return numCkovHoughMH;
  }

private:

  static constexpr auto pdgs = {211, 321, 2212};
  std::vector<float> mCkovTheoretical = {-1., -1., -1.};
  std::vector<int> sizeValues;
  std::vector<float> xValues, yValues, qValues,  thetaCerValues, phiCerValues, sigmaRingValues;


  std::vector<int> numRawClustersValues, rawSizeValues;

 
  std::vector<float> zProton, zKaon, zPion;
  std::vector<float> pProton, pKaon, pPion;
  bool mIsTrackToReconstruct = false;


  // ef > if we take into account that particle specie is known
  bool mIsTrackToReconstructWithKnownPdg = false;


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

  float mP;

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



  /*
  void ckovAngle(o2::dataformats::MatchInfoHMP* match, const
  std::vector<o2::hmpid::Cluster> clusters, int index, double nmean, float xRa,
  float yRa)
  {
    // Pattern recognition method based on Hough transform
    // Arguments:   pTrk     - track for which Ckov angle is to be found
    //              pCluLst  - list of clusters for this chamber
    //   Returns:            - track ckov angle, [rad],

    const int nMinPhotAcc = 3; // Minimum number of photons required to perform
  the pattern recognition

    int nClusTot = clusters.size();

    initVars(nClusTot);

    float xPc, yPc, th, ph;

    match->getHMPIDtrk(xPc, yPc, th, ph); // initialize this track: th and ph
  angles at middle of RAD

    setTrack(xRa, yRa, th, ph);

    fParam->setRefIdx(nmean);

    float mipQ = -1, mipX = -1, mipY = -1;
    int chId = -1, sizeClu = -1;

    fPhotCnt = 0;

    int nPads = 0;

    for (int iClu = 0; iClu < clusters.size(); iClu++) { // clusters loop

      o2::hmpid::Cluster cluster = clusters.at(iClu);
      nPads += cluster.size();
      if (iClu == index) { // this is the MIP! not a photon candidate: just
  store mip info mipX = cluster.x(); mipY = cluster.y(); mipQ = cluster.q();
        sizeClu = cluster.size();
        continue;
      }

      chId = cluster.ch();
      if (cluster.q() > 2 * fParam->qCut() || cluster.size() > 4) {
        continue;
      }
      double thetaCer, phiCer;
      if (findPhotCkov(cluster.x(), cluster.y(), thetaCer, phiCer)) { // find
  ckov angle for this  photon candidate fPhotCkov[fPhotCnt] = thetaCer; //
  actual theta Cerenkov (in TRS) fPhotPhi[fPhotCnt] = phiCer;
        fPhotClusIndex[fPhotCnt] = iClu; // actual phi   Cerenkov (in TRS): -pi
  to come back to "unusual" ref system (X,Y,-Z) fPhotCnt++; // increment counter
  of photon candidates
      }
    } // clusters loop

    match->setHMPIDmip(mipX, mipY, mipQ, fPhotCnt);     // store mip info in any
  case match->setIdxHMPClus(chId, index + 1000 * sizeClu); // set index of
  cluster match->setMipClusSize(sizeClu);

    if (fPhotCnt < nMinPhotAcc) {         // no reconstruction with <=3 photon
  candidates match->setHMPsignal(kNoPhotAccept); // set the appropriate flag
      return;
    }

    fMipPos.Set(mipX, mipY);

    // PATTERN RECOGNITION STARTED:
    if (fPhotCnt > fParam->multCut()) {
      fIsWEIGHT = kTRUE;
    } // offset to take into account bkg in reconstruction
    else {
      fIsWEIGHT = kFALSE;
    }

    float photCharge[10] = {0x0};

    int iNrec = flagPhot(houghResponse(), clusters, photCharge); // flag photons
  according to individual theta ckov with respect to most probable
    // int iNrec = flagPhot(houghResponse(), clusters); // flag photons
  according to individual theta ckov with respect to most probable

    match->setPhotCharge(photCharge);
    match->setHMPIDmip(mipX, mipY, mipQ, iNrec); // store mip info

    if (iNrec < nMinPhotAcc) {
      match->setHMPsignal(kNoPhotAccept); // no photon candidates are accepted
      return;
    }

    int occupancy = (int)(1000 * (nPads / (6. * 80. * 48.)));

    double thetaC = findRingCkov(clusters.size()); // find the best
  reconstructed theta Cherenkov findRingGeom(thetaC, 2);

    match->setHMPsignal(thetaC + occupancy); // store theta Cherenkov and
  chmaber occupancy
    // match->SetHMPIDchi2(fCkovSigma2); //store experimental ring angular
  resolution squared

    // deleteVars(); ef : in case of smart-pointers, should not be necessary?
  } // CkovAngle()
  */
};

#endif
