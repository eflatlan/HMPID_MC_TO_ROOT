#include "SimpleCluster.cpp"

#ifndef Alisigma2_H
#define Alisigma2_H

class Alisigma2_ {

public:
  Alisigma2_(double refIdx) { mRefIdx = refIdx; }

  /* Member function declarations

  #double sigma2(double trkTheta, double trkPhi, double ckovTh, double ckovPh);

  double sigLoc(double trkTheta, double trkPhi, double thetaC, double phiC,
  double betaM);

  double sigCrom(double trkTheta, double trkPhi, double thetaC, double phiC,
  double betaM);

  double sigGeom(double trkTheta, double trkPhi, double thetaC, double phiC,
  double betaM);

  double sigmaCorrFact(int iPart, double occupancy); */

  double getRefIdx() const { return mRefIdx; } // running refractive index

  double mRefIdx = 1.2904; // running refractive index

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
    return 1.5833;
  } //<--TEMPORAR--> to be removed in future. Mean refractive index of WIN
    // material (SiO2)

  // double winIdx() const { return 1.5787; } //<--TEMPORAR--> to be removed in
  // future. Mean refractive index of WIN material (SiO2)

  double gapIdx() const {
    return 1.0005;
  } //<--TEMPORAR--> to be removed in future. Mean refractive index of GAP
    // material (CH4)

  double sigma2(double trkTheta, double trkPhi, double ckovTh, double ckovPh) const
  {

    // Analithical calculation of total error (as a sum of localization,
    // geometrical and chromatic errors)

    // on Cerenkov angle for a given Cerenkov photon

    // created by a given MIP. Fromulae according to CERN-EP-2000-058

    // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]

    //            dip and azimuthal angles for MIP taken at the entrance to
    //            radiator, [radians]

    //            MIP beta

    //   Returns: absolute error on Cerenkov angle, [radians]

    double trkBeta = 1. / (std::cos(ckovTh) * getRefIdx());

    if (trkBeta > 1) {

      trkBeta = 1; // protection against bad measured thetaCer
    }

    if (trkBeta < 0) {

      trkBeta = 0.0001; //
    }

    double sigLocVar = sigLoc(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);

    double sigGeomVar = sigGeom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);

    double sigCromVar = sigCrom(trkTheta, trkPhi, ckovTh, ckovPh, trkBeta);

    const double trackErrorSigma = 0.002; // 2mRad sigma for track errors

    double sigmaRes =
        std::sqrt(sigLocVar * sigLocVar + sigGeomVar * sigGeomVar +
                  sigCromVar * sigCromVar + trackErrorSigma * trackErrorSigma);

    return sigmaRes;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  double sigLoc(double trkTheta, double trkPhi, double thetaC, double phiC,
                double betaM) const

  {

    // Analitical calculation of localization error (due to finite segmentation
    // of PC) on Cerenkov angle for a given

    // Cerenkov photon

    // created by a given MIP. Fromulae according to CERN-EP-2000-058

    // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]

    //            dip and azimuthal angles for MIP taken at the entrance to
    //            radiator, [radians]

    //            MIP beta

    //   Returns: absolute error on Cerenkov angle, [radians]

    double phiDelta = phiC - trkPhi;

    double sint = std::sin(trkTheta);

    double cost = std::cos(trkTheta);

    double sinf = std::sin(trkPhi);

    double cosf = std::cos(trkPhi);

    double sinfd = std::sin(phiDelta);

    double cosfd = std::cos(phiDelta);

    double tantheta = std::tan(thetaC);

    double alpha = cost - tantheta * cosfd * sint; // formula (11)

    double k = 1. - getRefIdx() * getRefIdx() +
               alpha * alpha / (betaM * betaM); // formula (after 8 in the text)

    if (k < 0) {

      return 1e10;
    }

    double mu = sint * sinf +
                tantheta * (cost * cosfd * sinf + sinfd * cosf); // formula (10)

    double e = sint * cosf +
               tantheta * (cost * cosfd * cosf - sinfd * sinf); // formula (9)

    double kk =
        betaM * std::sqrt(k) / (gapThick() * alpha); // formula (6) and (7)

    // formula (6)

    double dtdxc = kk * (k * (cosfd * cosf - cost * sinfd * sinf) -
                         (alpha * mu / (betaM * betaM)) * sint * sinfd);

    // formula (7)            pag.4

    double dtdyc = kk * (k * (cosfd * sinf + cost * sinfd * cosf) +
                         (alpha * e / (betaM * betaM)) * sint * sinfd);

    double errX = 0.2, errY = 0.25; // end of page 7

    return std::sqrt(errX * errX * dtdxc * dtdxc + errY * errY * dtdyc * dtdyc);
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  double sigCrom(double trkTheta, double trkPhi, double thetaC, double phiC,
                 double betaM) const

  {

    // Analitical calculation of chromatic error (due to lack of knowledge of
    // Cerenkov photon energy)

    // on Cerenkov angle for a given Cerenkov photon

    // created by a given MIP. Fromulae according to CERN-EP-2000-058

    // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]

    //            dip and azimuthal angles for MIP taken at the entrance to
    //            radiator, [radians]

    //            MIP beta

    //   Returns: absolute error on Cerenkov angle, [radians]

    double phiDelta = phiC - trkPhi;

    double sint = std::sin(trkTheta);

    double cost = std::cos(trkTheta);

    double cosfd = std::cos(phiDelta);

    double tantheta = std::tan(thetaC);

    double alpha = cost - tantheta * cosfd * sint; // formula (11)

    double dtdn =
        cost * getRefIdx() * betaM * betaM / (alpha * tantheta); // formula (12)

    //  double f = 0.00928*(7.75-5.635)/std::sqrt(12.);

    double f = 0.0172 * (7.75 - 5.635) / std::sqrt(24.);

    return f * dtdn;

  } // SigCrom()

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  double sigGeom(double trkTheta, double trkPhi, double thetaC, double phiC,
                 double betaM) const

  {

    // Analitical calculation of geometric error (due to lack of knowledge of
    // creation point in radiator)

    // on Cerenkov angle for a given Cerenkov photon

    // created by a given MIP. Formulae according to CERN-EP-2000-058

    // Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]

    //            dip and azimuthal angles for MIP taken at the entrance to
    //            radiator, [radians]

    //            MIP beta

    //   Returns: absolute error on Cerenkov angle, [radians]

    double phiDelta = phiC - trkPhi;

    double sint = std::sin(trkTheta);

    double cost = std::cos(trkTheta);

    double sinf = std::sin(trkPhi);

    double cosfd = std::cos(phiDelta);

    double costheta = std::cos(thetaC);

    double tantheta = std::tan(thetaC);

    double alpha = cost - tantheta * cosfd * sint; // formula (11)

    double k = 1. - getRefIdx() * getRefIdx() +
               alpha * alpha / (betaM * betaM); // formula (after 8 in the text)

    if (k < 0) {

      return 1e10;
    }

    double eTr = 0.5 * radThick() * betaM * std::sqrt(k) /
                 (gapThick() * alpha); // formula (14)

    double lambda = (1. - sint * sinf) * (1. + sint * sinf); // formula (15)

    double c1 =
        1. / (1. + eTr * k /
                       (alpha * alpha * costheta * costheta)); // formula (13.a)

    double c2 = betaM * std::pow(k, 1.5) * tantheta * lambda /
                (gapThick() * alpha * alpha); // formula (13.b)

    double c3 = (1. + eTr * k * betaM * betaM) /
                ((1 + eTr) * alpha * alpha); // formula (13.c)

    double c4 = std::sqrt(k) * tantheta * (1 - lambda) /
                (gapThick() * betaM); // formula (13.d)

    double dtdT = c1 * (c2 + c3 * c4);

    double trErr = radThick() / (std::sqrt(12.) * cost);

    return trErr * dtdT;

  } // SigGeom()
};
#endif
