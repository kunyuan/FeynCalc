#include "weight.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace diag;
using namespace std;

void weight::ReadDiagrams() {
  Pool.GPoolSize = 0;
  Pool.Ver4PoolSize = 0;

  int ID = 0;
  for (auto &name : Para.GroupName) {
    // construct filename based on format string and group id
    string FileName = fmt::format(Para.DiagFileFormat, name);
    ifstream DiagFile(FileName);
    ASSERT_ALLWAYS(DiagFile.is_open(),
                   "Unable to find the file " << FileName << endl);
    // group Group;
    LOG_INFO("Find " << FileName << "\n");
    // vector<green> GList;
    istream &DiagFileStream = DiagFile;
    Groups.push_back(ReadOneGroup(DiagFileStream, Pool));
    Groups.back().Name = name;
    Groups.back().ID = ID;
    ID++;
  }
  LOG_INFO("Find " << Pool.GPoolSize << " indepdent green's function.");
  LOG_INFO("Find " << Pool.Ver4PoolSize << " indepdent 4-vertex.");

  // cout << "After read" << endl;
  // cout << ToString(*(GroupList[0].DiagList[0].GIndex[0])) << endl;
  Initialization();
}

void weight::Initialization() {

  LOG_INFO("Initializating diagram states ...")
  for (auto &group : Groups) {
    group.ReWeight = 1.0;
    for (auto &diag : group.Diag) {
      for (int i = 0; i < group.GNum; i++) {
        diag.G[i]->Excited = false;
        diag.G[i]->Version = -1;
        diag.G[i]->Weight = 1.0e-10;
      }
      for (int i = 0; i < group.Ver4Num; i++) {
        diag.Ver4[i]->Excited = {false, false};
        diag.Ver4[i]->Version = -1;
        diag.Ver4[i]->Weight = {1.0e-10, -1.0e-10};
      }
    }

    if (Para.ObsType == EQUALTIME) {
      // for (int i = 0; i < group.Ver4Num; ++i)
      for (int i = 0; i < group.TauNum; ++i)
        if (group.IsExtTau[i])
          // to measure equal-time observable, lock all external tau
          group.IsLockedTau[i] = true;
    }
  }

  LOG_INFO("Initializating MC variables ...")
  // initialize momentum variables
  for (auto &mom : Var.LoopMom)
    for (int i = 0; i < D; i++)
      mom[i] = Random.urn() * Para.Kf / sqrt(D);

  // initialize tau variables
  for (int i = 0; i < MaxTauNum / 2; i++) {
    Var.Tau[2 * i] = Random.urn() * Para.Beta;
    Var.Tau[2 * i + 1] = Var.Tau[2 * i]; // assume even and odd tau are the same
  }

  // initialize spin variables
  for (auto &sp : Var.LoopSpin)
    sp = (spin)(Random.irn(0, 1));

  // initialize external momentum
  for (int i = 0; i < ExtMomBinSize; i++) {
    // the external momentum only has x component
    Var.ExtMomTable[i][0] = Para.MinExtMom + i*(Para.MaxExtMom-Para.MinExtMom)/ExtMomBinSize;
    for (int j = 1; j < D; j++)
      Var.ExtMomTable[i][j] = 0.0;
  }
  Var.CurrExtMomBin = 0;
  // Var.LoopMom[0].fill(0.0);
  // for (int i = 0; i < D; i++)
  //   Var.LoopMom[0][i] = Var.ExtMomTable[Var.CurrExtMomBin][i];
  Var.LoopMom[0] = Var.ExtMomTable[Var.CurrExtMomBin];

  // initialize external tau
  Var.Tau[0] = 0.0;
  Var.Tau[1] = 1.0e-10; // do not make Tau[1]==Tau[0], otherwise the Green's
                        // function is not well-defined

  Var.CurrTau = Var.Tau[1] - Var.Tau[0];

  // initialize group

  Var.CurrVersion = 0;
  //   Var.CurrGroup = &Groups[0];

  Var.CurrGroup = &Groups[0];
  LOG_INFO("Calculating the weights of all objects...")

  ChangeGroup(*Var.CurrGroup, true);
  GetNewWeight(*Var.CurrGroup);
  AcceptChange(*Var.CurrGroup);

  LOG_INFO("Initializating variables done.")
}

void weight::ChangeGroup(group &Group, bool Forced) {
  // the objects (G, Ver or Ver4) in the new group will be recalculated if the
  // either of the following conditions is met: 1) Forced=true, then all objects
  // are forced recalculated 2) object.Version<CurrVersion, means the objects
  // are not in the current group, and are already outdated
  for (auto &d : Group.Diag) {
    // cout << "diag ID: " << d.ID << endl;
    for (int i = 0; i < Group.GNum; i++) {
      // cout << "G: " << i << endl;
      green *G = d.G[i];
      if (Forced || G->Version < Var.CurrVersion) {
        double Tau = Var.Tau[G->TauBasis[OUT]] - Var.Tau[G->TauBasis[IN]];
        G->Excited = true;
        GetMom(G->LoopBasis, Group.LoopNum, _Mom);
        G->NewWeight = Fermi.Green(Tau, _Mom, UP, G->Type);
        // if (Group.ID == 2) {
        //   cout << Group.Name << ", " << G->NewWeight << endl;
        //   ;
        // }
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      // cout << "Ver: " << i << endl;
      vertex4 *Ver4 = d.Ver4[i];
      if (Para.UseVer4) {
        if (Forced || Ver4->Version < Var.CurrVersion) {
          Ver4->Excited = {true, true};
          GetMom(Ver4->LoopBasis[INL], Group.LoopNum, _InL);
          GetMom(Ver4->LoopBasis[OUTL], Group.LoopNum, _OutL);
          GetMom(Ver4->LoopBasis[INR], Group.LoopNum, _InR);
          GetMom(Ver4->LoopBasis[OUTR], Group.LoopNum, _OutR);
          VerFunc.Vertex4(_InL, _InR, _OutL, _OutR, 0, 0,
                          Ver4->NewWeight[DIRECT], Ver4->NewWeight[EXCHANGE]);
        }
      } else {
        if (Forced || Ver4->Version < Var.CurrVersion) {
          Ver4->Excited = {true, true};
          GetMom(Ver4->IntLoopBasis[IN], Group.LoopNum, _Mom);
          Ver4->NewWeight[IN] = Bose.Interaction(0.0, _Mom, Ver4->Type[IN]);
          GetMom(Ver4->IntLoopBasis[OUT], Group.LoopNum, _Mom);
          Ver4->NewWeight[OUT] = Bose.Interaction(0.0, _Mom, Ver4->Type[OUT]);
        }
      }
    }
  }
}

void weight::ChangeMom(group &Group, int MomIndex) {
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      if (G->LoopBasis[MomIndex] != 0) {
        double Tau = Var.Tau[G->TauBasis[OUT]] - Var.Tau[G->TauBasis[IN]];
        G->Excited = true;
        GetMom(G->LoopBasis, Group.LoopNum, _Mom);
        G->NewWeight = Fermi.Green(Tau, _Mom, UP, G->Type);
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      vertex4 *Ver4 = d.Ver4[i];

      if (Para.UseVer4) {
        if (Ver4->LoopBasis[INL][MomIndex] != 0 ||
            Ver4->LoopBasis[INR][MomIndex] != 0 ||
            Ver4->LoopBasis[OUTL][MomIndex] != 0 ||
            Ver4->LoopBasis[OUTR][MomIndex] != 0) {
          Ver4->Excited = {true, true};
          GetMom(Ver4->LoopBasis[INL], Group.LoopNum, _InL);
          GetMom(Ver4->LoopBasis[OUTL], Group.LoopNum, _OutL);
          GetMom(Ver4->LoopBasis[INR], Group.LoopNum, _InR);
          GetMom(Ver4->LoopBasis[OUTR], Group.LoopNum, _OutR);
          VerFunc.Vertex4(_InL, _InR, _OutL, _OutR, 0, 0,
                          Ver4->NewWeight[DIRECT], Ver4->NewWeight[EXCHANGE]);
        }
      } else {
        if (Ver4->IntLoopBasis[DIRECT][MomIndex] != 0) {
          Ver4->Excited[DIRECT] = true;
          GetMom(Ver4->IntLoopBasis[IN], Group.LoopNum, _Mom);
          Ver4->NewWeight[DIRECT] =
              Bose.Interaction(0.0, _Mom, Ver4->Type[DIRECT]);
        }
        if (Ver4->IntLoopBasis[EXCHANGE][MomIndex] != 0) {
          Ver4->Excited[EXCHANGE] = true;
          GetMom(Ver4->IntLoopBasis[EXCHANGE], Group.LoopNum, _Mom);
          Ver4->NewWeight[EXCHANGE] =
              Bose.Interaction(0.0, _Mom, Ver4->Type[OUT]);
        }
      }
    }
  }
}

void weight::ChangeTau(group &Group, int TauIndex) {
  // TODO: we assume TauLeft==TauRight for now
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      int TauIn = G->TauBasis[IN];
      int TauOut = G->TauBasis[OUT];
      bool ReCalcFlag = false;

      if (TauIndex == TauIn || TauIndex == TauOut) {
        // trigger recalculation
        double Tau = Var.Tau[TauOut] - Var.Tau[TauIn];
        G->Excited = true;
        GetMom(G->LoopBasis, Group.LoopNum, _Mom);
        G->NewWeight = Fermi.Green(Tau, _Mom, UP, G->Type);
      }
    }
  }
}

double weight::GetNewWeight(group &Group) {
  static double VIn, VOut, TotWeight;
  Group.NewWeight = 0.0;

  for (auto &d : Group.Diag) {
    double GWeight = d.SymFactor;
    for (int i = 0; i < Group.GNum; i++) {
      if (d.G[i]->Excited)
        GWeight *= d.G[i]->NewWeight;
      else
        GWeight *= d.G[i]->Weight;
    }

    double VerWeight;

    if (Group.Ver4Num == 0) {
      VerWeight = d.SpinFactor[0];
      // cout << "spin factor: " << d.SpinFactor[0] << endl;
    } else {
      vertex4 *Ver4 = d.Ver4[0];

      _Tree[0][0] = Ver4->Excited[DIRECT] ? Ver4->NewWeight[DIRECT]
                                          : Ver4->Weight[DIRECT];
      _Tree[0][1] = Ver4->Excited[EXCHANGE] ? Ver4->NewWeight[EXCHANGE]
                                            : Ver4->Weight[EXCHANGE];

      // if (Group.ID == 2) {
      //   cout << "0, 0:  " << _Tree[0][0] << ", 0, 1: " << _Tree[0][1] <<
      //   endl;
      //   // exit(0);
      // }

      int BlockNum = 2;
      for (int level = 1; level < Group.Ver4Num; level++) {

        vertex4 *Ver4 = d.Ver4[level];
        VIn = Ver4->Excited[DIRECT] ? Ver4->NewWeight[DIRECT]
                                    : Ver4->Weight[DIRECT];
        VOut = Ver4->Excited[EXCHANGE] ? Ver4->NewWeight[EXCHANGE]
                                       : Ver4->Weight[EXCHANGE];
        // if (Group.ID == 2) {
        //   cout << "Vin " << VIn << ", Vout: " << VOut << endl;
        //   // exit(0);
        // }

        for (int j = 0; j < BlockNum; j++) {
          _Tree[level][2 * j] = _Tree[level - 1][j] * VIn;
          _Tree[level][2 * j + 1] = _Tree[level - 1][j] * VOut;
        }
        BlockNum *= 2;
      }

      VerWeight = 0.0;
      for (int j = 0; j < BlockNum; j++)
        VerWeight += _Tree[Group.Ver4Num - 1][j] * d.SpinFactor[j];

      //============= for spin case ===========================//

      // double TempWeightIn, TempWeightOut;
      // for (int i = 0; i < Group.Ver4Num; i++) {
      //   //=========== for spinless case ===========================//
      //   vertex *Ver4 = d.Ver4[i];
      //   TempWeightIn =
      //       Ver4->Excited[IN] ? Ver4->NewWeight[IN] : Ver4->Weight[IN];
      //   TempWeightOut =
      //       Ver4->Excited[OUT] ? Ver4->NewWeight[OUT] : Ver4->Weight[OUT];
      //   VerWeight *= TempWeightIn - TempWeightOut;

      // if (Group.ID == 2) {
      //   cout << GWeight << ", " << VerWeight << endl;
      //   // exit(0);
      // }
    }

    d.NewWeight = GWeight * VerWeight / pow(2 * PI, D * Group.InternalLoopNum);

    // Group Weight= sum of all diagram weight in the group
    Group.NewWeight += d.NewWeight;
  }
  // if (Group.ID == 2) {
  //   cout << Group.NewWeight << endl;
  //   // exit(0);
  // }
  return Group.NewWeight;
}

void weight::AcceptChange(group &Group) {
  Var.CurrVersion++;
  Var.CurrGroup = &Group;
  Group.Weight = Group.NewWeight; // accept group  newweight

  for (auto &d : Group.Diag) {
    d.Weight = d.NewWeight; // accept diagram newweight
    for (int i = 0; i < Group.GNum; i++) {
      green *G = d.G[i];
      G->Version = Var.CurrVersion;
      if (G->Excited) {
        G->Excited = false;
        G->Weight = G->NewWeight;
      }
    }
    for (int i = 0; i < Group.Ver4Num; i++) {
      vertex4 *Ver4 = d.Ver4[i];
      Ver4->Version = Var.CurrVersion;
      if (Ver4->Excited[DIRECT]) {
        Ver4->Excited[DIRECT] = false;
        Ver4->Weight[DIRECT] = Ver4->NewWeight[DIRECT];
      }
      if (Ver4->Excited[EXCHANGE]) {
        Ver4->Excited[EXCHANGE] = false;
        Ver4->Weight[EXCHANGE] = Ver4->NewWeight[EXCHANGE];
      }
    }
  }
}

void weight::RejectChange(group &Group) {
  for (auto &d : Group.Diag) {
    for (int i = 0; i < Group.GNum; i++) {
      if (d.G[i]->Excited)
        d.G[i]->Excited = false;
      for (int i = 0; i < Group.Ver4Num; i++) {
        if (d.Ver4[i]->Excited[0])
          d.Ver4[i]->Excited[0] = false;
        if (d.Ver4[i]->Excited[1])
          d.Ver4[i]->Excited[1] = false;
      }
    }
  }
}

void weight::GetMom(const loop &LoopBasis, const int &LoopNum, momentum &Mom) {
  // In C++11, because of the move semantics, there is no additional cost by
  // returning an array

  auto &loopmom = Var.LoopMom;
  for (int d = 0; d < D; ++d)
    Mom[d] = loopmom[0][d] * LoopBasis[0];

  for (int i = 1; i < LoopNum; ++i)
    for (int d = 0; d < D; ++d)
      Mom[d] += loopmom[i][d] * LoopBasis[i];
}

bool weight::IsInteractionReducible(loop &LoopBasisVer, int LoopNum) {
  // check if an interaction is reducible
  if ((!Equal(LoopBasisVer[0], 1.0)) && (!Equal(LoopBasisVer[0], -1.0)))
    return false;

  bool Flag = true;
  for (int i = 1; i < LoopNum; i++) {
    if (!Equal(LoopBasisVer[i], 0.0)) {
      Flag = false;
      break;
    }
  }
  return Flag;
};

bool weight::IsInteractionReducible(loop &LoopBasisG1, loop &LoopBasisG2,
                                    int LoopNum) {
  // check if an interaction is reducible
  if ((!Equal(LoopBasisG1[0] - LoopBasisG2[0], 1.0)) &&
      (!Equal(LoopBasisG1[0] - LoopBasisG2[0], -1.0)))
    return false;

  bool Flag = true;
  for (int i = 1; i < LoopNum; i++) {
    if (!Equal(LoopBasisG1[i] - LoopBasisG2[i], 0.0)) {
      Flag = false;
      break;
    }
  }
  return Flag;
};