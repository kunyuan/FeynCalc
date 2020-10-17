#ifndef diagram_H
#define diagram_H

#include "global.h"
#include "utility/utility.h"
#include <array>
#include <string>
#include <vector>

extern parameter Para;

namespace diag {
using namespace std;
const size_t MaxBranchNum = 1 << (MaxOrder - 1); // 2**(MaxOrder-1)

// column-major two dimensional array
template <class T, size_t ROW, size_t COL>
using matrix = std::array<std::array<T, COL>, ROW>;

typedef std::array<double, MaxLoopNum>
    loop; // array to store the loop basis for a propagator or interaction line
typedef std::array<int, 2> tau; // array to store the tau basis (In and Out)
                                // for a propagator or interaction line

// G pool to store all basis for Green's functions
struct green {
  long long int Version; // keep track of the version
  int Type;              // type of each green's function
  loop LoopBasis;        // loop basis for momentum
  tau TauBasis;          // tau basis
  double Weight;         // weight of each green's function
  double NewWeight;      // weight of each green's function
  bool Excited;
};

// Ver pool to store all basis for 4-vertex functions
struct vertex4 {
  long long int Version; // keep track of the version
  array<int, 2> Type;    // type of each vertex function
  array<green *, 4>
      Ver4Legs; // the GIndex of four legs of every indepdent 4-vertex
  array<loop, 4> LoopBasis;    // loop basis for momentum transfer
  array<loop, 2> IntLoopBasis; // interaction loop basis for momentum transfer
  // tau TauBasis;               //tau basis, Left and Right
  array<double, 2> Weight; // direct/exchange weight of each 4-vertex function
  array<double, 2> NewWeight;
  // direct/exchange weight of each 4-vertex function
  array<bool, 2> Excited;
  // status of direct/exchange weight of each 4-vertex function
};

struct pool {
  std::array<green, MaxGPoolSize> GPool; // array to store indepdent G
  std::array<vertex4, MaxVerPoolSize>
      Ver4Pool; // array to store indepdent vertex4
  int GPoolSize;
  int Ver4PoolSize;
};

struct diagram {
  int ID;
  array<int, MaxGNum> Permutation;
  array<int, 4> LegTau;                   // tau index for the 4 legs
  double SymFactor;                       // the symmetry factor of a diagram
  array<double, MaxBranchNum> SpinFactor; // the spin factor of a diagram
  array<green *, 2 * MaxOrder> G;         // the index of all indepdent G
  array<vertex4 *, 2 * MaxOrder> Ver4;    // the index of all indepdent 4-vertex
  double Weight;
  double NewWeight;
};

// A group can be diagrams with different orders,
// or diagrams with same order but have little sign cancellation
// store G and Ver indexes pointing to the corresponding pool
struct group {
  std::string Name;
  int ID;
  int HugenNum;        // Number of Hugenholz diagrams in each group
  int Order;           // diagram order of the group
  int Ver4Num;         // number of 4-vertex
  int GNum;            // number of G
  int LoopNum;         // dimension of loop basis
  int InternalLoopNum; // dimension of internal loop basis
  int ExtLoopNum;      // dimension of external loop basis
  int TauNum;          // dimension of tau basis
  int ExtTauNum;       // dimension of external tau basis
  int InternalTauNum;  // dimension of internal tau basis
  double ReWeight;
  double Weight;
  double NewWeight;
  array<bool, MaxLoopNum> IsExtLoop;
  array<bool, MaxLoopNum> IsLockedLoop;
  array<bool, MaxTauNum> IsExtTau;
  array<bool, MaxTauNum> IsLockedTau;
  vector<diagram> Diag; // diagrams
};

// diagram type in the group
// const int HUGEN = 1;
// const int NORMAL = 1;
// const int RG = 2;

group ReadOneGroup(istream &, pool &);

void Test(group &);

}; // namespace diag

std::string ToString(const diag::diagram &);
std::string ToString(const diag::group &);
std::string ToString(const diag::green &);
std::string ToString(const diag::vertex4 &);

#endif