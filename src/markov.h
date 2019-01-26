#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/estimator.h"
#include "utility/rng.h"
#include "weight.h"
#include <string>
#include <vector>

namespace mc {
using namespace std;
const int MCUpdates = 5;

typedef array<double, ExtMomBinSize> polar;

class markov {
public:
  markov() : Var(Weight.Var), Groups(Weight.Groups){};
  long long Counter;

  void Initialization(std::string FilePrefix);
  void Hop(const int);
  void PrintMCInfo();
  void AdjustGroupReWeight();

  void Measure();
  void SaveToFile(std::string FilePrefix);

  int DynamicTest();

private:
  diag::weight Weight;
  diag::variable &Var;
  vector<diag::group> &Groups;
  // double SumofProbofCall[MCUpdates];
  std::string UpdateName[MCUpdates];
  double Accepted[MCUpdates][MaxGroupNum];
  double Proposed[MCUpdates][MaxGroupNum];

  void ChangeTau();
  void ChangeMomentum();
  void IncreaseOrder();
  void DecreaseOrder();
  void ChangeGroup();

  double ShiftExtK(const int &, int &);
  double ShiftK(const momentum &, momentum &);
  double ShiftTau(const double &, double &);

  double GetNewTau(double &);
  double GetNewK(momentum &);
  double RemoveOldTau(double &);
  double RemoveOldK(momentum &);

  enum Updates {
    INCREASE_ORDER = 0,
    DECREASE_ORDER,
    CHANGE_GROUP,
    CHANGE_TAU,
    CHANGE_MOM,
    END
  };
  std::string _DetailBalanceStr(Updates op);

  // polarizatoin for each group
  vector<polar> Polar;
  // polarizatoin for each group at the zero momentumr;
  vector<double> PolarStatic;
};
}; // namespace mc

#endif
