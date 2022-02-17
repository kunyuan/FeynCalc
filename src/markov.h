#ifndef markov_H
#define markov_H

#include "global.h"
#include "utility/rng.h"
#include "weight.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace mc {
using namespace std;
const int MCUpdates = 5;

typedef array<double, ExtMomBinSize> polar_q;
typedef array<array<double, ExtTauBinSize>, ExtMomBinSize> polar;

class markov {
public:
  markov();
  long long Counter;

  void PrintMCInfo();
  void PrintDeBugMCInfo();
  void AdjustGroupReWeight();

  // MC updates
  void ChangeTau();
  void ChangeMomentum();
  void IncreaseOrder();
  void DecreaseOrder();
  void ChangeGroup();

  void Measure();
  void SaveToFile();

  int DynamicTest();

  // MC variables
  diag::weight Weight;
  diag::variable &Var;
  vector<diag::group> &Groups;

private:
  // polarization for each group at external Momentum
  unordered_map<int, polar_q> PolarStatic;

  // polarization for each group at the external Tau;
  unordered_map<int, polar> Polar;

  // polarization for each diagrams
  unordered_map<string, polar_q> Polar_Diag;

  // kinetic energy for each group
  unordered_map<int, polar_q> E_k;

  // MC updates

  double ShiftExtK(const int &, int &);
  double ShiftK(const momentum &, momentum &);
  double ShiftExtTau(const int &, int &);
  double ShiftTau(const double &, double &);

  double GetNewTau(double &);
  double GetNewK(momentum &);
  double RemoveOldTau(double &);
  double RemoveOldK(momentum &);

  // MC updates information
  std::string UpdatesName[MCUpdates];
  double Accepted[MCUpdates][MaxGroupNum];
  double Proposed[MCUpdates][MaxGroupNum];
  // temp of the last Proposed[3] for AdjustGroupReWeight
  double Temp_Proposed[MaxGroupNum]; 

  enum Updates {
    INCREASE_ORDER = 0,
    DECREASE_ORDER,
    CHANGE_GROUP,
    CHANGE_TAU,
    CHANGE_MOM,
    END
  };
  std::string _DetailBalanceStr(Updates op);
};
}; // namespace mc

#endif
