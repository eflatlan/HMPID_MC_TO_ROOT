#include <cmath>
#include <iostream>
#include <utility> // For std::pair

#include <vector>
#include <cmath>
#include <algorithm>

#pragma once

class FilterPhotons
{
  public:

    // Function to filter out tracks that do not have enough photons in the relevant Ckov zone from pdg
    std::pair<int, int> expectedNumPhotons(double n, double p, int pdg) {
      double MASS;

      // Set MASS based on the pdg-code
      if (std::abs(pdg) == 211) {
          MASS = 0.1396;
      } else if (std::abs(pdg) == 321) {
          MASS = 0.4937;
      } else if (std::abs(pdg) == 2212) {
          MASS = 0.9383;
      } else {
          std::cerr << "Invalid PDG code" << std::endl;
          return std::make_pair(-1, -1); // Error case
      }

      double pLim = MASS / std::sqrt(n * n - 1);
      double pSquared = p * p;

      double pSquaredPlusMSquared = pSquared + MASS * MASS;

      bool pMask = p > pLim;
      if (!pMask) {
          std::cerr << "Momentum below threshold" << std::endl;
          return std::make_pair(999, 99); // Error case
      }

      double cosThetaC = std::sqrt(pSquaredPlusMSquared) / (p * n);
      double sinThetaC = std::sin(std::acos(cosThetaC));

      // Saturation level for expected number of photons
      const int NUM_EXPECTED_PHOTONS_SATURATION = 13;
      double numExpectedPhotons = sinThetaC * sinThetaC * NUM_EXPECTED_PHOTONS_SATURATION / 0.4;

      Printf("pdg %d num exp phot %.3f", pdg, numExpectedPhotons);

      // Calculate limits for evaluation using Poisson distribution mean - std
      int limitEvaluate = std::floor(numExpectedPhotons - std::sqrt(numExpectedPhotons));
      int limitEvaluateHalf = std::floor(numExpectedPhotons / 2 - std::sqrt(numExpectedPhotons / 2));

      if (limitEvaluateHalf < 1) {
        limitEvaluateHalf = 1;
      }

      if (limitEvaluate < 3) {
        limitEvaluate = 3;
      }      

      return std::make_pair(limitEvaluate, limitEvaluateHalf);
    }



    // ef > TODO moce all these to central  HMP libs


    // dn/dE is 0.0172 eV−1
    // σE = (dn/dE) * σ_det^E
    // σE = 6.33 × 10−4

    // EF > TODO make the sigmaRefIndex on sigmaRefIndex = sigmaRefIndex/sqrt(expectedPhotonsPerSpecie) ??
    double momThresh(double mass, double refIndex){

      // refIndex is mean value of refractive index;
      // we should adapt teh mom-threshold to be 
      // n = meanRefIndex - 2*sigmaRefIndex 
      auto meanRefIndex = refIndex; // ef > todo change notation to avoid disambiguity >> IP to func meanRefIndex
      auto sigmaRefIndex = 0.0055; // ef > TODO make this more accurate and not hardcoded
      auto limitRefIndex = meanRefIndex - 2*sigmaRefIndex;

      return mass/std::sqrt(limitRefIndex*limitRefIndex-1);
    }





    // ef > TODO take into account ring average ang res?
    bool  evaluateMomentumThreshold(bool& momThreshPi, bool& momThreshKa, bool& momThreshPr, double p, double n){

      const float PION_MASS = 0.1396;
      const float KAON_MASS = 0.4937;
      const float PROTON_MASS = 0.9383;

      momThreshPi = p > momThresh(PION_MASS, n);
      momThreshKa = p > momThresh(KAON_MASS, n);
      momThreshPr = p > momThresh(PROTON_MASS, n);

      return (momThreshPi || momThreshKa || momThreshPr);      
    }



    bool evaluateNumPhotThresh(const std::vector<float>& pr_phi, const std::vector<float>& k_phi, const std::vector<float>& pi_phi, double n, double p, int pdgCodeTrack, bool fixedPdgIsNumPhotsOk) {

      

      // ef > can we do this?
      // only reconstructing when we also now the specie-type?



      // Lambda function to count non-zero elements
      auto countNonZero = [](const std::vector<float>& vec) {
          return std::count_if(vec.begin(), vec.end(), [](float val) { return val != 0; });
      };

      // Lambda function to count photons in the half plane
      auto countInHalfPlane = [](const std::vector<float>& vec) {
          return std::count_if(vec.begin(), vec.end(), [](float phi) {
              // return std::abs(phi) > M_PI / 2 && std::abs(phi) < M_PI;
              return phi < 0 / 2 && phi > - M_PI;
          });
      };

      // Lambda function to count photons in the half plane
      auto countInHalfPlane2 = [](const std::vector<float>& vec) {
          return std::count_if(vec.begin(), vec.end(), [](float phi) {
              return std::abs(phi) > M_PI / 2 && std::abs(phi) < M_PI;
          });
      };




      // Calculate the number of photons for each particle type
      int numProtons = countNonZero(pr_phi);
      int numKaons = countNonZero(k_phi);
      int numPions = countNonZero(pi_phi);

      // Calculate photons in the half plane for protons, kaons, and pions
      int numPhotonInHalfProton = countInHalfPlane(pr_phi);
      int numPhotonInHalfKaon = countInHalfPlane(k_phi);
      int numPhotonInHalfPion = countInHalfPlane(pi_phi);


      Printf("numPhotonInHalfPion %d", numPhotonInHalfPion);
      Printf("numPhotonInHalfKaon %d", numPhotonInHalfKaon);
      Printf("numPhotonInHalfProton %d", numPhotonInHalfProton);


      int numPhotonInHalfProton2 = countInHalfPlane2(pr_phi);
      int numPhotonInHalfKaon2 = countInHalfPlane2(k_phi);
      int numPhotonInHalfPion2 = countInHalfPlane2(pi_phi);

      Printf("numPhotonInHalfPion2 %d", numPhotonInHalfPion2);
      Printf("numPhotonInHalfKaon2 %d", numPhotonInHalfKaon2);
      Printf("numPhotonInHalfProton2 %d", numPhotonInHalfProton2);


      // Calculate thresholds for each particle type
      auto [numPhotThreshProton, numPhotHalfThreshProton] = expectedNumPhotons(n, p, 2212);
      auto [numPhotThreshKaon, numPhotHalfThreshKaon] = expectedNumPhotons(n, p, 321);
      auto [numPhotThreshPion, numPhotHalfThreshPion] = expectedNumPhotons(n, p, 211);

      bool momThreshPi = false, momThreshKa = false, momThreshPr = false;

      // ef > TODO take into account ring average ang res?
      evaluateMomentumThreshold(momThreshPi, momThreshKa, momThreshPr, p, n);

      Printf("momThreshPi %d Ka %d Pr %d", momThreshPi, momThreshKa, momThreshPr);


      Printf("Pr %d hlf %d mom %d", (numProtons >= numPhotThreshProton) , (numPhotonInHalfProton2 >= numPhotHalfThreshProton), momThreshPr);
      Printf("Ka %d hlf %d mom %d", (numKaons >= numPhotThreshKaon) , (numPhotonInHalfKaon2 >= numPhotHalfThreshKaon), momThreshKa);
      Printf("Pi %d hlf %d mom %d", (numPions >= numPhotThreshPion) , (numPhotonInHalfPion2 >= numPhotHalfThreshPion), momThreshPi);


      // Check if thresholds are met for each particle type
      bool condPr = (numProtons >= numPhotThreshProton) && (numPhotonInHalfProton2 >= numPhotHalfThreshProton) && (momThreshPr);
      bool condKa = (numKaons >= numPhotThreshKaon) && (numPhotonInHalfKaon2 >= numPhotHalfThreshKaon) && (momThreshKa);
      bool condPi = (numPions >= numPhotThreshPion) && (numPhotonInHalfPion2 >= numPhotHalfThreshPion) && (momThreshPi);


      Printf("condPr %d condKa %d condPi %d", condPr, condKa, condPi);


      // Combine conditions for all particle types
      bool combinedCondition = condPr || condKa || condPi;


      // ef> can we do this ? taking into account that we know hte specie-type
      switch (std::abs(pdgCodeTrack)) {
        case 211:
          fixedPdgIsNumPhotsOk = condPi && momThreshPi;
        break;

        case 321:
          fixedPdgIsNumPhotsOk = condKa && momThreshKa;
        break;
          
        case 2212:
          fixedPdgIsNumPhotsOk = condPr && momThreshPr;
        break;
      }


      return combinedCondition;
    }

};