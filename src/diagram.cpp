#include "diagram.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/utility.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

using namespace diag;
using namespace std;

/*Find Case Insensitive Sub String in a given substring */
size_t _findCaseInsensitive(std::string data, std::string toSearch,
                            size_t pos = 0) {
  // Convert complete given String to lower case
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  // Convert complete given Sub String to lower case
  std::transform(toSearch.begin(), toSearch.end(), toSearch.begin(), ::tolower);
  // Find sub string in given string
  return data.find(toSearch, pos);
}

int _findIndex(vector<int> permutation, int key) {
  vector<int>::iterator itr = find(permutation.begin(), permutation.end(), key);
  int index = 0;

  if (itr != permutation.cend())
    index = std::distance(permutation.begin(), itr);
  else
    ABORT("the key " << key << " doesn't exist in the permutation!");
  return index;
}

string _CheckKeyWord(istream &file, string KeyWord) {
  // a line only contains space will be skipped
  string buff;
  do {
    getline(file, buff);
    if (file.bad())
      ABORT("Fail to read the file!");
  } while (buff.find_first_not_of(' ') == buff.npos);

  auto found = _findCaseInsensitive(buff, KeyWord);
  ASSERT_ALLWAYS(found != string::npos,
                 fmt::format("{0} is not in: {1}", KeyWord, buff));

  return buff;
}

template <typename T>
vector<T> _ExtractNumbers(istream &file, string KeyWord = "") {
  // This function extracts the first T type number in one line in the file
  // stream.
  string line;
  getline(file, line);

  if (KeyWord.size() > 0) {
    auto found = _findCaseInsensitive(line, KeyWord);
    ASSERT_ALLWAYS(found != string::npos,
                   fmt::format("{0} is not in: {1}", KeyWord, line));
  }

  stringstream ss;
  vector<T> IntList;
  /* Storing the whole string into string stream */
  ss << line;
  /* Running loop till the end of the stream */
  string temp;
  T found;
  while (!ss.eof()) {
    /* extracting word by word from stream */
    ss >> temp;
    /* Checking the given word is integer or not */
    if (stringstream(temp) >> found)
      IntList.push_back(found);
    /* To save from space at the end of string */
    temp = "";
  }
  return IntList;
}

vector<loop> _Transpose(vector<vector<double>> &Basis) {
  vector<loop> NewBasis;
  if (Basis.size() == 0)
    return NewBasis;

  for (int i = 0; i < Basis[0].size(); i++) {
    loop TempLoopBasis;
    TempLoopBasis.fill(0);

    for (int j = 0; j < Basis.size(); j++)
      TempLoopBasis[j] = Basis[j][i];

    NewBasis.push_back(TempLoopBasis);
  }
  return NewBasis;
}

green *_AddOneGToPool(pool &Pool, green &Green) {
  vector<green *> GPool_Type;
  // select all green in GPooll have the type as GType[i]
  for (int i = 0; i < Pool.GPoolSize; i++) {
    green *g = &(Pool.GPool[i]);
    if (g->Type == Green.Type)
      GPool_Type.push_back(g);
  }
  // select all green's function with same TauIn and TauOut
  vector<green *> GPool_Tau;
  for (auto g : GPool_Type) {
    if (g->TauBasis[IN] == Green.TauBasis[IN] &&
        g->TauBasis[OUT] == Green.TauBasis[OUT]) {
      GPool_Tau.push_back(g);
    }
  }
  // select all green's function with same LoopBasis
  vector<green *> GPool_Basis;
  for (auto g : GPool_Tau)
    if (Equal<double>(g->LoopBasis.data(), Green.LoopBasis.data(), MaxLoopNum))
      GPool_Basis.push_back(g);

  //   cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum <<
  //   endl;

  // if GPool_Filter is not empty, means the green's function already exists
  if (GPool_Basis.size() > 0) {
    for (auto &g : GPool_Basis) {
      //   cout << ToString(Green) << endl;
      //   cout << ToString(*g) << endl;
    }
    ASSERT_ALLWAYS(GPool_Basis.size() == 1,
                   "There are more than two same G in the GPool!");
    return GPool_Basis[0];
  } else {
    ASSERT_ALLWAYS(Pool.GPoolSize < MaxGPoolSize, "MaxGPoolSize is too small!");
    Pool.GPoolSize++;
    Pool.GPool[Pool.GPoolSize - 1] = Green;
    return &Pool.GPool[Pool.GPoolSize - 1];
  }
}

vertex4 *_AddOneInteractionToPool(pool &Pool, vertex4 &Vertex4) {
  vector<vertex4 *> Ver4Pool_Type;
  // select all vertex in VerPooll have the type as Vertex

  for (int i = 0; i < Pool.Ver4PoolSize; i++) {
    vertex4 *v = &(Pool.Ver4Pool[i]);
    if (v->Type[0] == Vertex4.Type[0])
      if (v->Type[1] == Vertex4.Type[1])
        Ver4Pool_Type.push_back(v);
  }

  vector<vertex4 *> Ver4Pool_Basis;
  for (auto v : Ver4Pool_Type)
    if (Equal<double>(v->IntLoopBasis[0].data(), Vertex4.IntLoopBasis[0].data(),
                      MaxLoopNum) &&
        Equal<double>(v->IntLoopBasis[1].data(), Vertex4.IntLoopBasis[1].data(),
                      MaxLoopNum))
      Ver4Pool_Basis.push_back(v);
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (Ver4Pool_Basis.size() > 0) {
    ASSERT_ALLWAYS(Ver4Pool_Basis.size() == 1,
                   "There are more than two same Vertex in the VerPool!");
    return Ver4Pool_Basis[0];
  } else {
    // add new vertex4 to pool
    ASSERT_ALLWAYS(Pool.Ver4PoolSize < MaxVerPoolSize,
                   "MaxVerPoolSize is too small!");
    Pool.Ver4PoolSize++;
    Pool.Ver4Pool[Pool.Ver4PoolSize - 1] = Vertex4;
    return &Pool.Ver4Pool[Pool.Ver4PoolSize - 1];
  }
}

vertex4 *_AddOneVer4ToPool(pool &Pool, vertex4 &Vertex4) {
  vector<vertex4 *> Ver4Pool_Type;
  // select all vertex in VerPooll have the type as Vertex
  for (int i = 0; i < Pool.Ver4PoolSize; i++) {
    vertex4 *v = &(Pool.Ver4Pool[i]);
    if (v->Type[0] == Vertex4.Type[0])
      if (v->Type[1] == Vertex4.Type[1])
        Ver4Pool_Type.push_back(v);
  }

  vector<vertex4 *> Ver4Pool_Basis;
  for (auto v : Ver4Pool_Type) {
    bool Flag = true;
    for (int leg = 0; leg < 4; leg++)
      if (!Equal<double>(v->LoopBasis[leg].data(),
                         Vertex4.LoopBasis[leg].data(), MaxLoopNum))
        Flag = false;
    if (Flag)
      Ver4Pool_Basis.push_back(v);
  }
  // cout << "Basis pool: " << GPool_Basis.size() << "," << MaxLoopNum << endl;
  // if GPool_Filter is not empty, means the green's function already exists
  if (Ver4Pool_Basis.size() > 0) {
    ASSERT_ALLWAYS(Ver4Pool_Basis.size() == 1,
                   "There are more than two same vertex4 in the Ver4Pool!");
    return Ver4Pool_Basis[0];
  } else {
    // add new vertex4 to pool
    ASSERT_ALLWAYS(Pool.Ver4PoolSize < MaxVerPoolSize,
                   "MaxVerPoolSize is too small!");
    Pool.Ver4PoolSize++;
    Pool.Ver4Pool[Pool.Ver4PoolSize - 1] = Vertex4;
    return &Pool.Ver4Pool[Pool.Ver4PoolSize - 1];
  }
}

vector<green *> _AddAllGToPool(pool &Pool, vector<tau> &VerBasis,
                               vector<loop> &LoopBasis, vector<int> &GType,
                               int GNum) {
  vector<green *> GIndex;
  for (int i = 0; i < GNum; i++) {
    // construct a new green's function
    green Green;
    Green.Type = GType[i];
    Green.LoopBasis.fill(0);
    std::copy(LoopBasis[i].begin(), LoopBasis[i].end(),
              Green.LoopBasis.begin());
    Green.TauBasis = VerBasis[i];
    // GList.push_back(Green);
    GIndex.push_back(_AddOneGToPool(Pool, Green));
  }
  //   cout << "All First: " << ToString(GPool[0]) << endl;
  //   cout << "All Last: " << ToString(GPool.back()) << endl;
  //   cout << "All First: " << ToString(*(GIndex[0])) << endl;
  //   cout << "All Last: " << ToString(*(GIndex.back())) << endl;
  //   cout << GIndex[0] << " vs " << &(GPool[0]) << endl;
  return GIndex;
}

vector<vertex4 *> _AddAllVer4ToPool(pool &Pool, vector<tau> &VerBasis,
                                    vector<int> &Ver4Legs,
                                    vector<loop> &LoopBasisG,
                                    vector<int> &VerType, int Ver4Num) {
  vector<vertex4 *> Ver4Index;
  for (int i = 0; i < Ver4Num; i++) {
    int Inidx = 2 * i, Outidx = 2 * i + 1;
    // construct a new vertex 4 function
    vertex4 Vertex4;
    Vertex4.Type = {VerType[Inidx], VerType[Outidx]};

    // build 4-leg loop basis
    for (int leg = 0; leg < 4; leg++) {
      int legidx = 4 * i + leg; // index shift
      Vertex4.LoopBasis[leg].fill(0);
      int gidx = Ver4Legs[legidx];
      std::copy(LoopBasisG[gidx].begin(), LoopBasisG[gidx].end(),
                Vertex4.LoopBasis[leg].begin());
    }

    // build Interaction loop basis
    Vertex4.IntLoopBasis[DIRECT].fill(0);
    Vertex4.IntLoopBasis[EXCHANGE].fill(0);
    for (int i = 0; i < Vertex4.IntLoopBasis[DIRECT].size(); ++i) {
      Vertex4.IntLoopBasis[DIRECT][i] =
          Vertex4.LoopBasis[INL][i] - Vertex4.LoopBasis[OUTL][i];
      Vertex4.IntLoopBasis[EXCHANGE][i] =
          Vertex4.LoopBasis[INL][i] - Vertex4.LoopBasis[OUTR][i];
    }

    if (Para.UseVer4)
      Ver4Index.push_back(_AddOneVer4ToPool(Pool, Vertex4));
    else
      Ver4Index.push_back(_AddOneInteractionToPool(Pool, Vertex4));
  }
  return Ver4Index;
}

diagram ReadOneDiagram(istream &DiagFile, pool &Pool, int Order, int LoopNum,
                       int GNum, int Ver4Num) {
  string buff;
  diagram Diagram;

  Diagram.G.fill(nullptr);
  Diagram.Ver4.fill(nullptr);

  //////// Diagram Topology  ////////////////////////
  _CheckKeyWord(DiagFile, "Permutation"); // title
  vector<int> Permutation = _ExtractNumbers<int>(DiagFile);
  std::copy(Permutation.begin(), Permutation.end(),
            Diagram.Permutation.begin());

  Diagram.LegTau[INL] = Permutation[0];
  Diagram.LegTau[OUTL] = _findIndex(Permutation, 0);
  Diagram.LegTau[INR] = Permutation[1];
  Diagram.LegTau[OUTR] = _findIndex(Permutation, 1);
  // cout << Diagram.LegTau[OUTL] << ", " << Diagram.LegTau[OUTR] << endl;

  //////// symmetry factor //////////////////
  _CheckKeyWord(DiagFile, "SymFactor"); // title
  Diagram.SymFactor = _ExtractNumbers<double>(DiagFile)[0];

  //////// Propagator type //////////////////
  _CheckKeyWord(DiagFile, "GType"); // title
  auto GType = _ExtractNumbers<int>(DiagFile);

  /////// Ver Basis  /////////////////////////
  _CheckKeyWord(DiagFile, "VertexBasis"); // title
  vector<int> StartVer = _ExtractNumbers<int>(DiagFile);
  vector<int> EndVer = _ExtractNumbers<int>(DiagFile);
  vector<tau> VerBasis;
  for (int i = 0; i < GNum; i++)
    VerBasis.push_back(tau({StartVer[i], EndVer[i]}));

  /////// Loop Basis  /////////////////////////
  _CheckKeyWord(DiagFile, "LoopBasis"); // title
  vector<vector<double>> TransposedLoopBasis;
  for (int j = 0; j < LoopNum; j++)
    TransposedLoopBasis.push_back(_ExtractNumbers<double>(DiagFile));

  vector<loop> LoopBasis = _Transpose(TransposedLoopBasis);

  /////// 4 legs of 4-ver  /////////////////////////
  _CheckKeyWord(DiagFile, "Ver4Legs"); // title
  vector<int> Ver4Legs;
  if (Ver4Num > 0)
    Ver4Legs = _ExtractNumbers<int>(DiagFile);

  /////// Interaction Type  /////////////////////////
  _CheckKeyWord(DiagFile, "WType"); // title
  vector<int> VerType;
  if (Ver4Num > 0)
    VerType = _ExtractNumbers<int>(DiagFile);

  /////// Spin Factor  ////////////////////
  _CheckKeyWord(DiagFile, "SpinFactor"); // title
  auto SpinFactor = _ExtractNumbers<double>(DiagFile);
  std::copy(
      SpinFactor.begin(), SpinFactor.end(),
      Diagram.SpinFactor.begin()); // copy spin factor into the group member

  ////////   Add G to GPool ////////////////////////
  vector<green *> GIndex =
      _AddAllGToPool(Pool, VerBasis, LoopBasis, GType, GNum);
  ASSERT_ALLWAYS(GIndex.size() == GNum,
                 "Number of Green's function does not match!");
  std::copy(GIndex.begin(), GIndex.end(), Diagram.G.begin());

  //   cout << "in diagram" << endl;
  //   cout << ToString(*(Diagram.GIndex[0])) << endl;
  //   cout << Diagram.GIndex[0] << endl;

  if (Ver4Num > 0) {
    //////  Add 4-Vertex to Ver4Pool /////////
    vector<vertex4 *> Ver4Index = _AddAllVer4ToPool(
        Pool, VerBasis, Ver4Legs, LoopBasis, VerType, Ver4Num);
    std::copy(Ver4Index.begin(), Ver4Index.end(), Diagram.Ver4.begin());
  }

  return Diagram;
}

group diag::ReadOneGroup(istream &DiagFile, pool &Pool) {
  group Group;
  _CheckKeyWord(DiagFile, "Type"); // group type, simply skip

  Group.HugenNum = _ExtractNumbers<int>(DiagFile, "DiagNum")[0];
  ASSERT_ALLWAYS(Group.HugenNum <= MaxDiagNum,
                 "Diagram Number must be smaller than " << MaxDiagNum);

  Group.Order = _ExtractNumbers<int>(DiagFile, "Order")[0];
  ASSERT_ALLWAYS(Group.Order <= MaxOrder,
                 "Order Number must be smaller than " << MaxOrder);

  Group.GNum = _ExtractNumbers<int>(DiagFile, "GNum")[0];
  ASSERT_ALLWAYS(Group.GNum <= MaxGNum,
                 "G Number must be smaller than " << MaxGNum);

  Group.Ver4Num = _ExtractNumbers<int>(DiagFile, "Ver4Num")[0];
  ASSERT_ALLWAYS(Group.Ver4Num <= MaxVer4Num,
                 "Ver4 Number must be smaller than " << MaxVer4Num);

  Group.LoopNum = _ExtractNumbers<int>(DiagFile, "LoopNum")[0];
  ASSERT_ALLWAYS(Group.LoopNum <= MaxLoopNum,
                 "Loop Number must be smaller than " << MaxLoopNum);

  vector<int> ExtLoop = _ExtractNumbers<int>(DiagFile, "ExtLoopIndex");
  Group.ExtLoopNum = ExtLoop.size();
  Group.IsExtLoop.fill(false);
  for (auto index : ExtLoop)
    Group.IsExtLoop[index] = true;

  vector<int> LockedLoop = _ExtractNumbers<int>(DiagFile, "DummyLoopIndex");
  Group.IsLockedLoop.fill(false);
  for (auto index : LockedLoop)
    Group.IsLockedLoop[index] = true;

  Group.TauNum = _ExtractNumbers<int>(DiagFile, "TauNum")[0];
  ASSERT_ALLWAYS(Group.TauNum <= MaxTauNum,
                 "Tau Number must be smaller than " << MaxTauNum);

  vector<int> ExtTau = _ExtractNumbers<int>(DiagFile, "ExtTauIndex");
  Group.ExtTauNum = ExtTau.size();
  Group.IsExtTau.fill(false);
  for (auto index : ExtTau)
    Group.IsExtTau[index] = true;

  vector<int> LockedTau = _ExtractNumbers<int>(DiagFile, "DummyTauIndex");
  Group.IsLockedTau.fill(false);
  for (auto index : LockedTau)
    Group.IsLockedTau[index] = true;

  Group.InternalLoopNum = Group.LoopNum - Group.ExtLoopNum;
  Group.InternalTauNum = Group.TauNum - Group.ExtTauNum;

  for (int i = 0; i < Group.HugenNum; i++) {
    diagram Diagram = ReadOneDiagram(DiagFile, Pool, Group.Order, Group.LoopNum,
                                     Group.GNum, Group.Ver4Num);
    Diagram.ID = i;
    Group.Diag.push_back(Diagram);
    // for (int i = 0; i < Group.GNum; i++)
    //   cout << ToString(*(Group.DiagList.back().GIndex[i])) << endl;
  }
  diag::Test(Group);
  //   cout << "Group" << endl;
  //   for (int i = 0; i < Group.GNum; i++)
  //   cout << "in Group" << endl;
  //   cout << ToString(Group) << endl;
  //   cout << ToString(*(Group.DiagList[0].GIndex[0])) << endl;

  return Group;
}

std::string ToString(const diagram &Diag) {
  std::ostringstream oss;
  oss << "DiagID: " << Diag.ID << "\n Weight=" << Diag.Weight
      << "\n SymFactor=" << Diag.SymFactor << endl;
  return oss.str();
}
std::string ToString(const group &Group) {
  return fmt::format("GroupID: {0}\n Weight={1}\n NewWeight={2}\n "
                     "HugenNum={3}\n ReWeight={4}\n",
                     Group.ID, Group.Weight, Group.NewWeight, Group.HugenNum,
                     Group.ReWeight);
};
std::string ToString(const green &G) {
  std::ostringstream oss;
  oss << "Version: " << G.Version << "\n Type=" << G.Type
      << "\n Weight=" << G.Weight
      << "\n TauBasis=" << ToString<int>(G.TauBasis.data(), 2)
      << "\n LoopBasis="
      << ToString<double>(G.LoopBasis.data(), (size_t)G.LoopBasis.size())
      << endl;
  return oss.str();
};
std::string ToString(const vertex4 &);

void diag::Test(group &Group) {
  for (auto &d : Group.Diag) {
    for (auto i = 0; i < Group.GNum; i++) {
      CHECKNULL(d.G[i]);
    }
    for (auto i = 0; i < Group.Ver4Num; i++) {
      CHECKNULL(d.Ver4[i]);
    }
  }
}
