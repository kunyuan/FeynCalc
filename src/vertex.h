#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"

double sum2(const momentum &);
double norm2(const momentum &);

namespace diag {

<<<<<<< HEAD
const int MAXSIGMABIN = 16385; //8193; //2**13+1
=======
const int MAXSIGMABIN = 32769; 
>>>>>>> dev_WDM
class fermi {
public:
  fermi();
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType);

private:
  // beyond which the expense sigma function will be called
  double UpperBound, LowerBound;
  double DeltaK;
  double Mu_ideal;
  double Mu_shift;
  double UpperBound2, LowerBound2; // lower upbound for better sigma
  double DeltaK2;
  double PhyGreen(double Tau, const momentum &Mom, bool IsFock);
  double TwoPhyGreen(double Tau, const momentum &Mom, bool IsFock);
  double ThreePhyGreen(double Tau, const momentum &Mom, bool IsFock);
  double FockSigma(const momentum &Mom);
  void BuildFockSigma();
  double Fock(double k);
  // warning: this function only works for T=0 and 3D!!!!
  double GetSigma(double k);
  double Sigma[MAXSIGMABIN];
  double Sigma2[MAXSIGMABIN];
};

class bose {
public:
  double Interaction(double Tau, const momentum &Momentum, int VerType);
};

class verfunc {
public:
  verfunc();
  void Vertex4(const momentum &InL, const momentum &InR, const momentum &OutL,
               const momentum &OutR, int Ver4TypeDirect, int Ver4TypeExchange,
               double &Direct, double &Exchange);

private:
  double Ver4AtUV[InInAngBinSize][InOutAngBinSize];
  double Angle2D(const momentum &K1, const momentum &K2);
  double Index2Angle(const int &Index, const int &AngleNum);
  int Angle2Index(const double &Angle, const int &AngleNum);
  void _TestAngleIndex();
  void _TestAngle2D();
};

}; // namespace diag
#endif
