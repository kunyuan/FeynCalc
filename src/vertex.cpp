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
<<<<<<< HEAD
  UpperBound = 12.0 * Para.Kf;  //6.0
=======
  UpperBound = 24.0 * Para.Kf;  //6.0
>>>>>>> dev_WDM
  LowerBound = 0.0;
  DeltaK = UpperBound / (MAXSIGMABIN-1);
  // LowerBound = DeltaK;
  // Mu_shift = Para.Ef;

  UpperBound2 = 1.2 * Para.Ef;
  LowerBound2 = 0.8 * Para.Ef;
  DeltaK2 = UpperBound2 / MAXSIGMABIN;
  if (Para.SelfEnergyType == FOCK)
     BuildFockSigma();
}

double fermi::Fock(double k) {
  // warning: this function only works for T=0!!!!
  double kF = Para.Kf;

  // if(k>UpperBound)  cout << 'k=' << k<< '>'<< UpperBound<< endl;
  if (D == 3) {
    double l = sqrt(Para.Mass2 + Para.Lambda);
    double fock = 1.0 + l / kF * atan((k - kF) / l);
    fock -= l / kF * atan((k + kF) / l);
    fock -= (l * l - k * k + kF * kF) / 4.0 / k / kF *
            log((l * l + (k - kF) * (k - kF)) / (l * l + (k + kF) * (k + kF)));
    fock *= (-2.0 * kF) / PI;

    // double shift = 1.0 - l / kF * atan(2.0 * kF / l);
    // shift -= l * l / 4.0 / kF / kF * log(l * l / (l * l + 4.0 * kF * kF));
    // shift *= (-2.0 * kF) / PI;
    // return fock - shift;
    return fock;
  } else if (D == 2) {
    double l2 = Para.Mass2 + Para.Lambda;
    double x = Para.Kf * Para.Kf + l2 - k * k;
    double c = 4.0 * k * k  * l2;
    double fock = -2.0 * log((sqrt(x * x + c) + x) / 2.0 / l2);

    // double shift = -2.0 * log((Para.Kf * Para.Kf + l2) / l2);
    // return fock - shift;
    return fock;
  }
}

void fermi::BuildFockSigma() {
  ASSERT_ALLWAYS(D == 3, "The Fock self energy is for 3D!");
  double fockE, k, beta, rs, lambda, kmax;
  int num_k;
  double dk;
  
  string FileName = "sigma3D.txt";
  ifstream FockFile(FileName);
  ASSERT_ALLWAYS(FockFile.is_open(),
                  "Unable to find the file " << FileName << endl);
  LOG_INFO("Find " << FileName << "\n");
  
  FockFile >> beta >> rs >> lambda;
  FockFile >> kmax >> num_k;
  FockFile >> Mu_shift >> Mu_ideal;
  dk = kmax/num_k;
  Mu_ideal = Mu_ideal* Para.Ef;
  if(!(Equal(beta/Para.Ef, Para.Beta, 1.0e-6)&&Equal(rs, Para.Rs, 1.0e-6)&&
     Equal(lambda, Para.Lambda, 1.0e-6)&&Equal(dk, DeltaK, 1.0e-6)&&
     (num_k>=MAXSIGMABIN))){
      ASSERT_ALLWAYS(false,
                  "The parameters in the file"<< FileName<<" unmatch." << endl);
      FockFile.close();
  }
  LOG_INFO("Read " << FileName << "with right parameters\n");
  for (int i=0; i < MAXSIGMABIN; i++){
    FockFile >> fockE;
    Sigma[i]=fockE-Mu_shift;
    k = i * DeltaK + LowerBound;
<<<<<<< HEAD
    cout << k << " : " << Sigma[i] << " vs " << Fock(k)-Mu_shift << endl;
=======
    // cout << k << " : " << Sigma[i] << " vs " << Fock(k)-Mu_shift << endl;
>>>>>>> dev_WDM
  }
  FockFile.close();

  //   for (int i = 0; i < MAXSIGMABIN; ++i) {
  //   // k: (0^+, UpBound^-)
  //   // i=0 ==> k==0.5*DeltaK
  //   // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
  //   k = (i + 0.5) * DeltaK + LowerBound;
  //   Sigma[i] = Fock(k);
  //   if (i > 0 && k <= LowerBound2 && k >= UpperBound2) {
  //     ASSERT_ALLWAYS(
  //         Equal(Sigma[i - 1], Sigma[i], 5.0e-5),
  //         fmt::format("Fock are not accurate enough! At k={0}: {1} vs {2}\n", k,
  //                     Sigma[i - 1], Sigma[i]));
  //   }
  //   cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  // }
  // for (int i = 0; i < MAXSIGMABIN; ++i) {
  //   // k: (0^+, UpBound^-)
  //   // i=0 ==> k==0.5*DeltaK
  //   // i=MAXSIGMABIN-1 ==> k==(MAXSIGMABIN-0.5)*DeltaK
  //   k = (i + 0.5) * DeltaK2 + LowerBound2;
  //   Sigma2[i] = Fock(k);
  //   if (i > 0) {
  //     ASSERT_ALLWAYS(Equal(Sigma2[i - 1], Sigma2[i], 5.0e-5),
  //                    fmt::format("The 2rd level Fock are not accurate enough!"
  //                                "level! At k={0}: {1} vs {2}\n",
  //                                k, Sigma2[i - 1], Sigma2[i]));
  //   }
  //   // cout << k << " : " << Sigma[i] << " vs " << Fock(k) << endl;
  // }
};

double fermi::FockSigma(const momentum &Mom) {
  double k = Mom.norm(); // bare propagator
  double fock;
  if (k >= LowerBound && k < UpperBound) {
    int i = (k - LowerBound) / DeltaK;
    fock = Sigma[i]+(k-DeltaK*i-LowerBound)*(Sigma[i+1]-Sigma[i])/DeltaK;
  } 
  // else if(k<LowerBound){
    // fock = Sigma[0]+(k-LowerBound)*(Sigma[1]-Sigma[0])/DeltaK;}
  else {
    fock = Fock(k) - Mu_shift;
  }
  // cout << k <<" : "<< fock <<" vs "<< Fock(k)-Mu_shift << endl;

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
  // double fock = Fock(k);
  return fock;
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
  else 
    Ek -= Mu_ideal;

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
  else 
    Ek -= Mu_ideal;

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
  else
    Ek -= Mu_ideal;

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
  if (Para.SelfEnergyType == FOCK)
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
