#include "vertex.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace diag;
using namespace std;

extern parameter Para;

// double sum2(const momentum &Mom) {
//   double Sum2 = 0.0;
//   for (int i = 0; i < D; i++)
//     Sum2 += Mom[i] * Mom[i];
//   return Sum2;
// }

// double norm2(const momentum &Mom) { return sqrt(sum2(Mom)); }

double bose::Interaction(double Tau, const momentum &Mom, int VerType) {
  if (VerType >= 0) {
    double interaction =
        8.0 * PI / (Mom.squaredNorm() + Para.Mass2 + Para.Lambda);
    if (VerType > 0) {
      // the interaction contains counter-terms
      interaction *=
          pow(Para.Lambda / (Mom.squaredNorm() + Para.Mass2 + Para.Lambda),
              VerType);
      // interaction *= pow(-1, VerType);
    }
    // cout << "Interaction: " << interaction << ", Mom: " << Mom.squaredNorm()
    //      << ", VerType: " << VerType << ", Lambda: " << Para.Lambda << endl;
    return interaction;
  } else if (VerType == -1) {
    return 1.0;
  } else if (VerType == -2) {
    return 0.0;
  } else {
    ABORT("VerType can not be " << VerType);
  }
}

fermi::fermi() {
  // UpperBound = 5.0 * Para.Ef;
  // LowerBound = 0.0;
  // DeltaK = UpperBound / MAXSIGMABIN;
  // UpperBound2 = 1.2 * Para.Ef;
  // LowerBound2 = 0.8 * Para.Ef;
  // DeltaK2 = UpperBound2 / MAXSIGMABIN;
  if (Para.SelfEnergyType == READ)
    BuildFockSigma();
  cout << "fermi initialized" << endl;
}

double fermi::Fock(double k) {
  // warning: this function only works for T=0!!!!
  double kF = Para.Kf;
  if (D == 3) {
    double l = sqrt(Para.Mass2 + Para.Lambda);
    double fock = 1.0 + l / kF * atan((k - kF) / l);
    fock -= l / kF * atan((k + kF) / l);
    fock -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
            log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
    fock *= (-2.0 * kF) / PI;

    double shift = 1.0 - l / kF * atan(2.0 * kF / l);
    shift -= l * l / 4.0 / kF / kF * log(l * l / (l * l + 4.0 * kF * kF));
    shift *= (-2.0 * kF) / PI;
    return fock - shift;
  } else if (D == 2) {
    double l2 = Para.Mass2 + Para.Lambda;
    double x = Para.Kf * Para.Kf + l2 - k * k;
    double c = 4.0 * k * k * l2;
    double fock = -2.0 * log((sqrt(x * x + c) + x) / 2.0 / l2);

    double shift = -2.0 * log((Para.Kf * Para.Kf + l2) / l2);
    return fock - shift;
  }
}

double fermi::BuildFockSigma() {
  double k, s;
  if (Para.SelfEnergyType == READ) {
    ifstream sigmafile("sigma.data");
    ASSERT_ALLWAYS(sigmafile.is_open(), "failed to load sigma.data!");
    // int size = -1;
    // cout << "reading ..." << endl;
    int size = 1024;
    sigmafile >> size;
    ASSERT_ALLWAYS(size >= 0, "error in reading sigma size!");
    kgrid.resize(size);
    Sigma.resize(size);
    for (int i = 0; i < size; ++i) {
      sigmafile >> k >> s >> dMu;
      //   // cout << "got: " << k << ", " << s << ", " << dMu << endl;
      //   // kgrid[i] = k * Para.Kf;
      // k = 12.0 * Para.Kf / 1023 * i;
      kgrid[i] = k * Para.Kf;
      // Sigma[i] = Fock(k + 1.0e-4);
      Sigma[i] = s;
    }
    LowerBound = kgrid[0];
    UpperBound = kgrid[kgrid.size() - 1];
    // cout << kgrid.size() << endl;
    // cout << Sigma.size() << endl;
    // cout << LowerBound / Para.Kf << " -> " << UpperBound / Para.Kf << endl;
    // cout << Sigma[0] << " -> " << Sigma[Sigma.size() - 1] << endl;
    // cout << kgrid[1022] << ", " << kgrid[1023] << endl;
    // cout << Sigma[1022] << ", " << Sigma[1023] << endl;
    sigmafile.close();
  }
  cout << "sigma initialized" << endl;
  return 0.0;
}

double fermi::FockSigma(const momentum &Mom) {
  double k = Mom.norm(); // bare propagator
  double fock = 0.0;
  // double fock;
  // if (k >= LowerBound2 && k < UpperBound2) {
  //   int i = (k - LowerBound2) / DeltaK2;
  //   fock = Sigma2[i];
  // } else if ((k >= LowerBound && k < LowerBound2) ||
  //            (k >= UpperBound2 && k < UpperBound)) {
  //   int i = (k - LowerBound) / DeltaK;
  //   fock = Sigma[i];
  // } else {
  //   fock = Fock(k);
  // }
  // ASSERT_ALLWAYS(
  //     Equal(fock, Fock(k), 5.0e-5),
  //     fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n",
  //     k,
  //                 fock, Fock(k)));
  if (Para.SelfEnergyType == FOCK) {
    fock = Fock(k);
  } else if (Para.SelfEnergyType == READ) {
    fock = LinearInterpolate(k);
    // cout << fock << endl;
  }
  return fock;
}

double fermi::LinearInterpolate(double k) {
  // if (k >= UpperBound || k < LowerBound)
  //   return 0.0;
  if (k >= UpperBound)
    return dMu;
  double dk = UpperBound / (kgrid.size() - 1);
  int xi0 = int(k / dk);
  // xi0 = floor(xgrid, x)
  int xi1 = xi0 + 1;
  if (xi1 == Sigma.size())
    return Sigma[Sigma.size() - 1];

  double dx0 = k - kgrid[xi0];
  double dx1 = kgrid[xi1] - k;

  double d0 = Sigma[xi0];
  double d1 = Sigma[xi1];
  double sigma = (Sigma[xi0] * dx1 + Sigma[xi1] * dx0) / (dx0 + dx1);
  ASSERT_ALLWAYS(sigma >= Sigma[xi0] && sigma <= Sigma[xi1],
                 "interpolation error!");
  return sigma + dMu;
}

double fermi::PhyGreen(double Tau, const momentum &Mom, bool IsFock) {
  // if tau is exactly zero, set tau=0^-
  double green, Ek;
  if (Tau == 0.0) {
    Tau = -1.0e-10;
  }

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  Ek = Mom.squaredNorm(); // bare propagator
  if (IsFock)
    Ek += FockSigma(Mom); // Fock diagram dressed propagator
  Ek -= Para.Mu;

  //// enforce an UV cutoff for the Green's function ////////
  // if(Ek>8.0*EF) then
  //   PhyGreen=0.0
  //   return
  // endif

  double x = Para.Beta * Ek / 2.0;
  double y = 2.0 * Tau / Para.Beta - 1.0;
  if (x > 100.0)
    green = exp(-x * (y + 1.0));
  else if (x < -100.0)
    green = exp(x * (1.0 - y));
  else
    green = exp(-x * y) / (2.0 * cosh(x));

  green *= s;

  // cout << "x: " << x << ", y: " << y << ", G: " << green << endl;
  // cout << "G: " << green << endl;

  if (std::isnan(green))
    ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
                  << ", Ek=" << Ek << ", Green=" << green << ", Mom"
                  << ToString(Mom));
  return green;
}

double fermi::TwoPhyGreen(double Tau, const momentum &Mom, bool IsFock) {
  // if tau is exactly zero, set tau=0^-
  // cout << Tau << endl;

  double green, Ek;
  if (Tau == 0.0) {
    Tau = -1.0e-10;
  }

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  Ek = Mom.squaredNorm(); // bare propagator
  if (IsFock)
    Ek += FockSigma(Mom); // Fock diagram dressed propagator
  Ek -= Para.Mu;

  // double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  // double y = 2.0 * Tau / Para.Beta - 1.0;
  // if (x > 100.0)
  //   green = exp(-x * (y + 1.0));
  // else if (x < -100.0)
  //   green = exp(x * (1.0 - y));
  // else
  //   green = exp(-x * y) / (2.0 * cosh(x));

  green = exp(-Ek * Tau) / pow((1.0 + exp(-Para.Beta * Ek)), 2.0) *
          (Tau - (Para.Beta - Tau) * exp(-Para.Beta * Ek));

  // green = 0.5 / (1.0 + cosh(Para.Beta * Ek)) * (-Para.Beta);

  green *= s;

  if (std::isfinite(green) == false)
    ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
                  << ", Ek=" << Ek << ", Green=" << green << ", Mom"
                  << ToString(Mom));
  // if (std::isnan(green))
  //   ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
  //                 << ", Ek=" << Ek << ", Green=" << green << ", Mom"
  //                 << ToString(Mom));
  return green;
}

double fermi::ThreePhyGreen(double Tau, const momentum &Mom, bool IsFock) {
  // if tau is exactly zero, set tau=0^-
  // cout << Tau << endl;

  double green, Ek;
  if (Tau == 0.0) {
    Tau = -1.0e-10;
  }

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  Ek = Mom.squaredNorm(); // bare propagator
  if (IsFock)
    Ek += FockSigma(Mom); // Fock diagram dressed propagator
  Ek -= Para.Mu;

  // double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  // double y = 2.0 * Tau / Para.Beta - 1.0;
  // if (x > 100.0)
  //   green = exp(-x * (y + 1.0));
  // else if (x < -100.0)
  //   green = exp(x * (1.0 - y));
  // else
  //   green = exp(-x * y) / (2.0 * cosh(x));
  if (Ek > 0.0) {
    double Factor = exp(-Para.Beta * Ek);
    green =
        exp(-Ek * Tau) / pow(1.0 + Factor, 3.0) *
        (Tau * Tau / 2.0 -
         (Para.Beta * Para.Beta / 2.0 + Para.Beta * Tau - Tau * Tau) * Factor +
         pow(Para.Beta - Tau, 2.0) * Factor * Factor / 2.0);
  } else {
    double Factor = exp(Para.Beta * Ek);
    green =
        exp(Ek * (Para.Beta - Tau)) / pow(1.0 + Factor, 3.0) *
        (Tau * Tau / 2.0 * Factor * Factor -
         (Para.Beta * Para.Beta / 2.0 + Para.Beta * Tau - Tau * Tau) * Factor +
         pow(Para.Beta - Tau, 2.0) / 2.0);
  }
  // cout << green << " vs " << green1 << endl;
  // if (abs(green - green1) > 1.0e-10) {
  //   cout << Para.Beta * Ek << endl;
  //   ABORT("wrong! " << green << ", " << green1 << endl);
  // }

  green *= s;

  if (std::isfinite(green) == false)
    ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
                  << ", Ek=" << Ek << ", Green=" << green << ", Mom"
                  << ToString(Mom));
  // if (std::isnan(green))
  //   ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
  //                 << ", Ek=" << Ek << ", Green=" << green << ", Mom"
  //                 << ToString(Mom));
  return green;
}

double fermi::Green(double Tau, const momentum &Mom, spin Spin, int GType) {
  double green;
  bool IsFock = false;
  if (Para.SelfEnergyType == FOCK || Para.SelfEnergyType == READ)
    IsFock = true;
  if (GType == 0) {
    green = PhyGreen(Tau, Mom, IsFock);
  } else if (GType == 1) {
    green = TwoPhyGreen(Tau, Mom, IsFock);
    // equal time green's function
    // green = PhyGreen(-1.0e-10, Mom);
  } else if (GType == 2) {
    green = ThreePhyGreen(Tau, Mom, IsFock);
  } else if (GType == -1) {
    green = PhyGreen(Tau, Mom, false);
    // green = 1.0;
  } else if (GType == -2) {
    // green = PhyGreen(Tau, Mom, IsFock);
    green = 1.0;
  } else {
    ABORT("GType " << GType << " has not yet been implemented!");
    // return FakeGreen(Tau, Mom);
  }
  return green;
}

verfunc::verfunc() {
  // test angle utility

  // TODO: implement D=3
  if (D == 3)
    return;

  _TestAngle2D();
  _TestAngleIndex();

  // initialize UV ver4 table
  momentum KInL = {1.0, 0.0};
  for (int inin = 0; inin < InInAngBinSize; ++inin)
    for (int inout = 0; inout < InOutAngBinSize; ++inout) {
      double AngleInIn = Index2Angle(inin, InInAngBinSize);
      double AngleInOut = Index2Angle(inout, InOutAngBinSize);
      momentum KInR = {cos(AngleInIn), sin(AngleInIn)};
      momentum KOutL = {cos(AngleInOut), sin(AngleInOut)};
      momentum KOutR = KInL + KInR - KOutL;
      Ver4AtUV[inin][inout] = (KInL - KInR).dot(KOutL - KOutR);
    }
}

void verfunc::Vertex4(const momentum &InL, const momentum &InR,
                      const momentum &OutL, const momentum &OutR,
                      int Ver4TypeDirect, int Ver4TypeExchange, double &Direct,
                      double &Exchange) {
  if (Ver4TypeDirect != 0 || Ver4TypeExchange != 0)
    ABORT("Ver4Type is only implemented for 0!");

  /**************   Yokawar Interaction ************************/
  Direct = 8.0 * PI / ((OutL - InL).squaredNorm() + Para.Mass2 + Para.Lambda);
  Exchange = 8.0 * PI / ((OutR - InL).squaredNorm() + Para.Mass2 + Para.Lambda);

  /**************   Generic Interaction ************************/
}

double verfunc::Angle2D(const momentum &K1, const momentum &K2) {
  // Returns the angle in radians between vectors 'K1' and 'K2'
  double dotp = K1.dot(K2);
  double det = K1[0] * K2[1] - K1[1] * K2[0];
  double Angle2D = atan2(det, dotp);
  if (Angle2D < 0)
    Angle2D += 2.0 * PI;
  return Angle2D;
}

double verfunc::Index2Angle(const int &Index, const int &AngleNum) {
  // Map index [0...AngleNum-1] to the theta range [0.0, 2*pi)
  return Index * 2.0 * PI / AngleNum;
}

int verfunc::Angle2Index(const double &Angle, const int &AngleNum) {
  // Map theta range  [0.0, 2*pi) to index [0...AngleNum-1]
  double dAngle = 2.0 * PI / AngleNum;
  if (Angle >= 2.0 * PI - dAngle / 2.0 || Angle < dAngle / 2.0)
    return 0;
  else
    return int(Angle / dAngle + 0.5);
}

void verfunc::_TestAngle2D() {
  // Test Angle functions
  momentum K1 = {1.0, 0.0};
  momentum K2 = {1.0, 0.0};

  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2)) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not zero! It is {:.13f}",
                  Angle2D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {-1.0, 0.0};
  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2) - PI) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not Pi! Instead, it is {:.13f}",
                  Angle2D(K1, K2)));

  K1 = {1.0, 0.0};
  K2 = {1.0, -EPS};
  ASSERT_ALLWAYS(
      abs(Angle2D(K1, K2) - 2.0 * PI) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not 2.0*Pi! It is {:.13f}",
                  Angle2D(K1, K2)));
}

void verfunc::_TestAngleIndex() {
  // Test Angle functions
  int AngleNum = 64;
  ASSERT_ALLWAYS(abs(Index2Angle(0, AngleNum) - 0.0) < 1.0e-10,
                 "Angle for index 0 should be zero!");

  ASSERT_ALLWAYS(abs(Index2Angle(AngleNum - 1, AngleNum) -
                     (2.0 * PI * (1.0 - 1.0 / AngleNum))) < 1.0e-10,
                 "Angle for index AngleNum should be 2.0*pi-0^+!");

  ASSERT_ALLWAYS(Angle2Index(0.0, AngleNum) == 0,
                 "Angle zero should have index 0!");
  ASSERT_ALLWAYS(
      Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) + EPS, AngleNum) == 0,
      "Angle 2*pi-pi/AngleNum should have index 1!");

  ASSERT_ALLWAYS(Angle2Index(2.0 * PI * (1.0 - 0.5 / AngleNum) - EPS,
                             AngleNum) == AngleNum - 1,
                 "Angle 2*pi-pi/AngleNum-0^+ should have index AngleNum!");
}