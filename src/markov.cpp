//
//  markov.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//
#include "markov.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

extern parameter Para;
extern RandomFactory Random;

using namespace mc;
using namespace diag;
using namespace std;

#define NAME(x) #x
// #define COPYFROMTO(x, y)                                                       \
//   for (int i = 0; i < D; i++)                                                  \
//     y[i] = x[i];

// x=2*i, then PAIR(x)=2*i+1
// x=2*i+1, then PAIR(x)=2*i
#define PAIR(x) (int(x / 2) * 4 + 1 - x)

markov::markov() : Var(Weight.Var), Groups(Weight.Groups)
{
  ///==== initialize Weight ============================//
  Weight.ReadDiagrams();

  //===== initialize updates related variable ==========//

  UpdatesName[INCREASE_ORDER] = NAME(INCREASE_ORDER);
  UpdatesName[DECREASE_ORDER] = NAME(DECREASE_ORDER);
  UpdatesName[CHANGE_GROUP] = NAME(CHANGE_GROUP);
  UpdatesName[CHANGE_MOM] = NAME(CHANGE_MOM);
  UpdatesName[CHANGE_TAU] = NAME(CHANGE_TAU);

  // for(int i=0;i<MCUpdates;i++)
  // UpdatesName[(Updates)i]=NAME((Updates))

  InitialArray(&Accepted[0][0], 1.0e-10, MCUpdates * MaxGroupNum);
  InitialArray(&Proposed[0][0], 1.0e-10, MCUpdates * MaxGroupNum);
  InitialArray(&Temp_Proposed[0], 0.0, MaxGroupNum);

  ///==== initialize observable =======================//
  for (auto &g : Groups)
  {
    PolarStatic[g.ID].fill(1.0e-10);
    for (int i = 0; i < ExtMomBinSize; i++)
      Polar[g.ID][i].fill(1.0e-10);
    E_k[g.ID].fill(1.0e-10);
    for (auto &d : g.Diag)
    {
      string CurrdiagName = g.Name + "_" + to_string(d.ID);
      Polar_Diag[CurrdiagName].fill(1.0e-10);
    }
  }
  ///=== Do all kinds of test  =======================//
  Weight.StaticTest();
  Weight.DynamicTest();

  ///==== Set Reweighting factor =====================//
  // AdjustGroupReWeight();
  for (int i = 0; i < Weight.Groups.size(); i++)
  {
    Weight.Groups[i].ReWeight = Para.ReWeight[i];
    for (int j = 0; j < ExtMomBinSize; j++)
    {
      if (Weight.Groups[i].Name.back() == '1' && Para.ObsType == FREQ_q && Para.MinExtMom != 0.0)
        Weight.Groups[i].ReWeight_ExtMom.push_back(pow(Var.ExtMomTable[j][0], 2.0));
      else if (Weight.Groups[i].Name.back() == '2' && Para.ObsType == FREQ_q && Para.MinExtMom != 0.0)
        Weight.Groups[i].ReWeight_ExtMom.push_back(pow(Var.ExtMomTable[j][0], 4.0));
      else
        Weight.Groups[i].ReWeight_ExtMom.push_back(1.0);
    }
  }
};

int markov::DynamicTest() { return Weight.DynamicTest(); }

void markov::PrintDeBugMCInfo()
{
  string msg;
  msg = string(80, '=') + "\n";
  msg += "\nMC Counter: " + to_string(Para.Counter) + ":\n";
  msg += "Current Group Info:\n " + ToString(*Var.CurrGroup);

  msg += string(80, '=') + "\n";
  msg += "GWeight: \n";
  for (int i = 0; i < Var.CurrGroup->GNum; i++)
    msg += ToString(Var.CurrGroup->Diag[0].G[i]->Weight) + "; ";
  msg += "\n";

  msg += "VerWeight: \n";
  for (int i = 0; i < Var.CurrGroup->Ver4Num; i++)
  {
    msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[0]) + ", ";
    msg += ToString(Var.CurrGroup->Diag[0].Ver4[i]->Weight[1]) + "; ";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "LoopMom: \n";
  for (int d = 0; d < D; d++)
  {
    for (int i = 0; i < Var.CurrGroup->LoopNum; i++)
      msg += ToString(Var.LoopMom[i][d]) + ", ";
    msg += "\n";
  }
  msg += "\n";

  msg += string(80, '=') + "\n";
  msg += "Tau: \n";
  for (int i = 0; i < Var.CurrGroup->TauNum; i++)
    msg += ToString(Var.Tau[i]) + ", ";

  msg += "\n";

  LOG_INFO(msg);
}

void markov::AdjustGroupReWeight()
{
  Para.ReWeight[0] = 1.0;
  Weight.Groups[0].ReWeight = Para.ReWeight[0];

  cout << "ReWeighting:" << endl;
  cout << 0 << "  " << Temp_Proposed[0] << "  " << Proposed[3][0] << endl;
  cout << Para.ReWeight[0] << endl;
  for (int i = 1; i < Weight.Groups.size(); i++)
  {
    cout << i << "  " << Temp_Proposed[i] << "  " << Proposed[3][i] << endl;
    Para.ReWeight[i] = (Proposed[3][0] - Temp_Proposed[0]) / (Proposed[3][i] - Temp_Proposed[i]) * Para.ReWeight[i] / 4.0;
    // if (i==2)
    //   Para.ReWeight[i] = Para.ReWeight[i]*4.0;
    Weight.Groups[i].ReWeight = Para.ReWeight[i];
    Temp_Proposed[i] = Proposed[3][i];
    cout << Para.ReWeight[i] << endl;
  }
  Temp_Proposed[0] = Proposed[3][0];
};

void markov::Measure()
{
  double MCWeight = fabs(Var.CurrGroup->Weight) * Var.CurrGroup->ReWeight * Var.CurrGroup->ReWeight_ExtMom[Var.CurrExtMomBin];
  double WeightFactor = Var.CurrGroup->Weight / MCWeight;

  PolarStatic[Var.CurrGroup->ID][Var.CurrExtMomBin] += WeightFactor;
  Polar[Var.CurrGroup->ID][Var.CurrExtMomBin][Var.CurrExtTauBin] += WeightFactor;
  // if (Var.CurrGroup->ID != 0 && Para.ObsType == KINETIC){
  // momentum _KMom;
  // for (auto &d : Var.CurrGroup->Diag){
  //   green *G = d.G[1];
  //   Weight.GetMom(G->LoopBasis, Var.CurrGroup->LoopNum, _KMom);
  // }
  // E_k[Var.CurrGroup->ID][Var.CurrExtMomBin] += WeightFactor* _KMom.squaredNorm();
  // }

  for (auto &d : Var.CurrGroup->Diag)
  {
    double WeightFactor_Diag = d.Weight / fabs(Var.CurrGroup->Weight) / Var.CurrGroup->ReWeight;
    WeightFactor_Diag /= Var.CurrGroup->ReWeight_ExtMom[Var.CurrExtMomBin];
    // if (d.Weight!=0){
    //   cout << Var.CurrGroup->ID <<"  "<< d.ID << " "<<d.Weight<< " "<< Var.CurrGroup->Weight<<endl;
    // }else{
    //   cout <<d.Weight<<"  error!"<<endl;
    //   exit(-1);
    // }
    string CurrdiagName = Var.CurrGroup->Name + "_" + to_string(d.ID);
    Polar_Diag[CurrdiagName][Var.CurrExtMomBin] += WeightFactor_Diag;
  }
};

void markov::SaveToFile()
{
  if (Para.ObsType == FREQ_tau)
  {
    for (int qi = 0; qi < ExtMomBinSize; qi++)
    {
      string FileName = fmt::format("pid{0}_q{1}.dat", Para.PID, qi);
      ofstream PolarFile;
      PolarFile.open(FileName, ios::out | ios::trunc);
      if (PolarFile.is_open())
      {
        PolarFile << "# Step: " << Para.Counter << endl;
        PolarFile << "# Group: ";
        for (auto &group : Groups)
          PolarFile << group.Name << ", ";
        PolarFile << endl;

        PolarFile << "# ReWeight: ";
        for (auto &group : Groups)
          PolarFile << group.ReWeight << ", ";
        PolarFile << endl;

        PolarFile << "# TauGrid: ";
        for (int j = 0; j < Polar[Groups[0].ID][qi].size(); j++)
          PolarFile << Var.ExtTauTable[j] << " ";
        PolarFile << endl;

        for (auto &group : Groups)
        {
          for (int j = 0; j < Polar[Groups[0].ID][qi].size(); j++)
            PolarFile << Polar[group.ID][qi][j] << " ";
          PolarFile << endl;
        }
      }
      else
      {
        LOG_WARNING("Polarization for ExtQ " << qi << " PID " << Para.PID << " fails to save!");
      }
      PolarFile.close();
    }
  }
  else
  {
    string FileName = fmt::format("pid{0}.dat", Para.PID);
    ofstream PolarFile;
    PolarFile.open(FileName, ios::out | ios::trunc);
    if (PolarFile.is_open())
    {
      PolarFile << "# Step: " << Para.Counter << endl;
      PolarFile << "# Group: ";
      for (auto &group : Groups)
        PolarFile << group.Name << ", ";
      PolarFile << endl;

      PolarFile << "# ReWeight: ";
      for (auto &group : Groups)
        PolarFile << group.ReWeight << ", ";
      PolarFile << endl;

      PolarFile << "# KGrid: ";
      for (int j = 0; j < PolarStatic[Groups[0].ID].size(); j++)
        PolarFile << Var.ExtMomTable[j][0] << " ";
      PolarFile << endl;

      for (auto &group : Groups)
      {
        for (int j = 0; j < PolarStatic[group.ID].size(); j++)
        {
          PolarFile << PolarStatic[group.ID][j] << " ";
        }
        PolarFile << endl;
      }
      // PolarFile << fmt::sprintf(
      //     "#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Group: %s, Step: %d\n",
      //     Para.PID, Para.ObsType, Para.Rs, Para.Beta, group.Name,
      //     Para.Counter);

      // for (int j = 0; j < Polar[group.ID].size(); j++)
      //   PolarFile << fmt::sprintf("%13.6f\t%13.6f\n", Var.ExtMomTable[j][0],
      //                             Polar[group.ID][j]);
    }
    else
    {
      LOG_WARNING("Polarization "
                  << " PID " << Para.PID << " fails to save!");
    }
    PolarFile.close();
  }
  // for (auto &group : Groups) {
  //   // ofstream PolarFile;
  //   // string FileName = fmt::format("group{0}_pid{1}.dat", group.Name, Para.PID);
  //   // PolarFile.open(FileName, ios::out | ios::trunc);
  //   // if (PolarFile.is_open()) {
  //   //   PolarFile << fmt::sprintf(
  //   //       "#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Group: %s, Step: %d\n",
  //   //       Para.PID, Para.ObsType, Para.Rs, Para.Beta, group.Name, Para.Counter);

  //   //   for (int j = 0; j < Polar[group.ID].size(); j++)
  //   //     PolarFile << fmt::sprintf("%13.6f\t%13.6f\n", Var.ExtMomTable[j][0],
  //   //                               Polar[group.ID][j]);
  //   //   PolarFile.close();
  //   // } else {
  //   //   LOG_WARNING("Polarization for PID " << Para.PID <<", Group "<<group.Name<< " fails to save!");
  //   // }
  //   for (auto &d : group.Diag){
  //     ofstream PolarFile;
  //     string CurrdiagName = group.Name+"_"+to_string(d.ID);
  //     string FileName = fmt::format("Diag{0}_pid{1}.dat", CurrdiagName, Para.PID);
  //     PolarFile.open(FileName, ios::out | ios::trunc);
  //     if (PolarFile.is_open()) {
  //       PolarFile << fmt::sprintf(
  //           "#PID:%d, Type:%d, rs:%.3f, Beta: %.3f, Diag: %s, Step: %d\n",
  //           Para.PID, Para.ObsType, Para.Rs, Para.Beta, group.Name+"_"+to_string(d.ID), Para.Counter);

  //       for (int j = 0; j < Polar_Diag[CurrdiagName].size(); j++)
  //         PolarFile << fmt::sprintf("%13.6f\t%13.6f\n", Var.ExtMomTable[j][0],
  //                                   Polar_Diag[CurrdiagName][j]);
  //       PolarFile.close();
  //     } else {
  //       LOG_WARNING("Polarization for PID " << Para.PID <<", Diag "<<CurrdiagName <<" fails to save!");
  //     }
  //   }
  // }

  // if (Para.ObsType == KINETIC){
  //   string fname = fmt::format("Ek_pid{0}.dat", Para.PID);
  //   ofstream EkFile;
  //   EkFile.open(fname, ios::out | ios::trunc);
  //   if (EkFile.is_open()) {
  //     EkFile << "# Step: " << Para.Counter << endl;
  //     EkFile << "# Group: ";
  //     for (auto &group : Groups)
  //       EkFile << group.Name << ", ";
  //     EkFile << endl;

  //     EkFile << "# ReWeight: ";
  //     for (auto &group : Groups)
  //       EkFile << group.ReWeight << ", ";
  //     EkFile << endl;

  //     EkFile << "# KGrid: ";
  //     for (int j = 0; j < E_k[Groups[0].ID].size(); j++)
  //       EkFile << Var.ExtMomTable[j][0] << " ";
  //     EkFile << endl;

  //     for (auto &group : Groups) {
  //       for (int j = 0; j < E_k[group.ID].size(); j++) {
  //         if (group.ID == 0)
  //           EkFile << Polar[group.ID][j] << " ";
  //         else
  //           EkFile << E_k[group.ID][j] << " ";
  //       }
  //       EkFile << endl;
  //     }
  //     EkFile.close();
  //   } else {
  //     LOG_WARNING("Kinetic energy for PID " << Para.PID << " fails to save!");
  //   }
  // }
}

void markov::ChangeGroup()
{
  group &NewGroup = Groups[Random.irn(0, Groups.size() - 1)];
  if (NewGroup.ID == Var.CurrGroup->ID)
    return;
  // cout << NewGroup.ID <<'    '<< Groups[NewGroup.ID].Name<<endl;
  Updates Name;
  double Prop = 1.0;

  if (NewGroup.Order == Var.CurrGroup->Order)
  {
    // change group with the same order
    Name = CHANGE_GROUP;
    Prop = 1.0;
  }
  else if (NewGroup.Order == Var.CurrGroup->Order + 1)
  {
    // change to a new group with one higher order
    Name = INCREASE_ORDER;
    static momentum NewMom;
    double NewTau;
    // Generate New Tau
    Prop = GetNewTau(NewTau);
    int NewTauIndex = Var.CurrGroup->TauNum;
    // ASSUME: NewTauIndex will never equal to 0 or 1
    Var.Tau[NewTauIndex] = NewTau;
    // Generate New Mom
    Prop *= GetNewK(NewMom);
    Var.LoopMom[Var.CurrGroup->LoopNum] = NewMom;
  }
  else if (NewGroup.Order == Var.CurrGroup->Order - 1)
  {
    // change to a new group with one lower order
    Name = DECREASE_ORDER;
    // Remove OldTau
    int TauToRemove = Var.CurrGroup->TauNum - 1;
    Prop = RemoveOldTau(Var.Tau[TauToRemove]);
    // Remove OldMom
    int LoopToRemove = Var.CurrGroup->LoopNum - 1;
    Prop *= RemoveOldK(Var.LoopMom[LoopToRemove]);
  }
  else
  {
    return;
  }

  Proposed[Name][Var.CurrGroup->ID] += 1;

  // if (NewGroup.ID == 3) {
  //   cout << fmt::format("group 3\n");
  //   cout << Weight.DebugInfo(NewGroup);
  // }

  Weight.ChangeGroup(NewGroup);
  double NewWeight = Weight.GetNewWeight(NewGroup) * NewGroup.ReWeight * NewGroup.ReWeight_ExtMom[Var.CurrExtMomBin];
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight) /
             Var.CurrGroup->ReWeight / Var.CurrGroup->ReWeight_ExtMom[Var.CurrExtMomBin];

  // if (NewGroup.ID == 2) {
  //   cout << endl << NewGroup.Name << ", TauNum: " << NewGroup.TauNum << endl;
  //   cout << fmt::format("weight: {0}, prop: {1}, R: {2}\n", NewWeight, Prop,
  //   R); cout << Weight.DebugInfo(NewGroup);
  // }

  if (Random.urn() < R)
  {
    Accepted[Name][Var.CurrGroup->ID]++;
    Weight.AcceptChange(NewGroup);
  }
  else
  {
    Weight.RejectChange(NewGroup);
  }
  return;
};

void markov::ChangeTau()
{
  // int TauIndex = Random.irn(0, Var.CurrGroup->TauNum);
  int TauIndex = Random.irn(1, Var.CurrGroup->TauNum);

  // if TauIndex is a locked tau, skip
  if (Var.CurrGroup->IsLockedTau[TauIndex])
    return;

  Proposed[CHANGE_TAU][Var.CurrGroup->ID]++;

  double CurrTau = Var.Tau[TauIndex];
  int NewExtTauBin;
  double Prop;

  if (Var.CurrGroup->IsExtTau[TauIndex] && Para.ObsType == FREQ_tau)
  {
    Prop = ShiftExtTau(Var.CurrExtTauBin, NewExtTauBin);
    Var.Tau[TauIndex] = Var.ExtTauTable[NewExtTauBin];
  }
  else
  {
    Prop = ShiftTau(CurrTau, Var.Tau[TauIndex]);
  }
  // Var.Tau[TauIndex] = NewTau;

  Weight.ChangeTau(*Var.CurrGroup, TauIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  // if (Var.CurrGroup->IsExtTau[TauIndex] && Para.ObsType == FREQ_tau)
  // cout << Var.CurrGroup->ID <<" " << NewWeight << " " << Var.CurrGroup->Weight <<endl;

  if (Random.urn() < R)
  {
    Accepted[CHANGE_TAU][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
    if (Var.CurrGroup->IsExtTau[TauIndex] && Para.ObsType == FREQ_tau)
      Var.CurrExtTauBin = NewExtTauBin;
  }
  else
  {
    // retore the old Tau if the update is rejected
    // if TauIndex is external, then its partner can be different
    Var.Tau[TauIndex] = CurrTau;
    Weight.RejectChange(*Var.CurrGroup);
  }
};

void markov::ChangeMomentum()
{
  int LoopIndex = Random.irn(0, Var.CurrGroup->LoopNum - 1);

  // if loopIndex is a locked loop, skip
  if (Var.CurrGroup->IsLockedLoop[LoopIndex])
    return;

  Proposed[CHANGE_MOM][Var.CurrGroup->ID]++;

  double Prop;
  int NewExtMomBin;
  static momentum CurrMom;

  CurrMom = Var.LoopMom[LoopIndex];

  if (Var.CurrGroup->IsExtLoop[LoopIndex])
  {
    Prop = ShiftExtK(Var.CurrExtMomBin, NewExtMomBin);
    Var.LoopMom[LoopIndex] = Var.ExtMomTable[NewExtMomBin];
    if (Var.LoopMom[LoopIndex].norm() > Para.MaxExtMom)
    {
      Var.LoopMom[LoopIndex] = CurrMom;
      return;
    }
    if (Var.LoopMom[LoopIndex].norm() < Para.MinExtMom)
    {
      Var.LoopMom[LoopIndex] = CurrMom;
      return;
    }
  }
  else
  {
    Prop = ShiftK(CurrMom, Var.LoopMom[LoopIndex]);
  }

  Weight.ChangeMom(*Var.CurrGroup, LoopIndex);
  double NewWeight = Weight.GetNewWeight(*Var.CurrGroup);
  double R = Prop * fabs(NewWeight) / fabs(Var.CurrGroup->Weight);
  if (Var.CurrGroup->IsExtLoop[LoopIndex])
    R = R * Var.CurrGroup->ReWeight_ExtMom[NewExtMomBin] / Var.CurrGroup->ReWeight_ExtMom[Var.CurrExtMomBin];

  if (Random.urn() < R)
  {
    Accepted[CHANGE_MOM][Var.CurrGroup->ID]++;
    Weight.AcceptChange(*Var.CurrGroup);
    if (Var.CurrGroup->IsExtLoop[LoopIndex])
    {
      Var.CurrExtMomBin = NewExtMomBin;
    }
  }
  else
  {
    Var.LoopMom[LoopIndex] = CurrMom;
    Weight.RejectChange(*Var.CurrGroup);
  }
};

double markov::GetNewTau(double &NewTau)
{
  NewTau = Random.urn() * Para.Beta;
  return Para.Beta;
};
double markov::RemoveOldTau(double &OldTau) { return 1.0 / Para.Beta; }

double markov::GetNewK(momentum &NewMom)
{
  //====== The hard Way ======================//
  //  double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  //  if (dK > Para.Kf / 2)
  //    dK = Para.Kf / 2; // to avoid dK>Kf situation
  double dK = Para.Kf / 2;
  double KAmp = Para.Kf + (Random.urn() - 0.5) * 2.0 * dK;
  // Kf-dK<KAmp<Kf+dK
  double Phi = 2.0 * PI * Random.urn();
  if (D == 3)
  {
    double Theta = PI * Random.urn();
    if (Theta == 0.0)
      return 0.0;
    double K_XY = KAmp * sin(Theta);
    NewMom[0] = K_XY * cos(Phi);
    NewMom[1] = K_XY * sin(Phi);
    NewMom[D - 1] = KAmp * cos(Theta);
    return 2.0 * dK                    // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI                  // prop density of Phi
           * PI                        // prop density of Theta
           * sin(Theta) * KAmp * KAmp; // Jacobian
  }
  else if (D == 2)
  {
    NewMom[0] = KAmp * cos(Phi);
    NewMom[1] = KAmp * sin(Phi);
    return 2.0 * dK   // prop density of KAmp in [Kf-dK, Kf+dK)
           * 2.0 * PI // prop density of Phi
           * KAmp;    // Jacobian
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   NewMom[i] = Para.Kf * (Random.urn() - 0.5) * 2.0;
  // return pow(2.0 * Para.Kf, D);

  //============================================//
};

double markov::RemoveOldK(momentum &OldMom)
{
  //====== The hard Way ======================//
  //  double dK = Para.Kf / sqrt(Para.Beta) / 4.0;
  //  if (dK > Para.Kf / 2)
  //    dK = Para.Kf / 2; // to avoid dK>Kf situation
  double dK = Para.Kf / 2;
  double KAmp = OldMom.norm();
  if (KAmp < Para.Kf - dK || KAmp > Para.Kf + dK)
    // Kf-dK<KAmp<Kf+dK
    return 0.0;
  if (D == 3)
  {
    auto SinTheta = sqrt(OldMom[0] * OldMom[0] + OldMom[1] * OldMom[1]) / KAmp;
    if (SinTheta < EPS)
      return 0.0;
    return 1.0 / (2.0 * dK * 2.0 * PI * PI * SinTheta * KAmp * KAmp);
  }
  else if (D == 2)
  {
    return 1.0 / (2.0 * dK * 2.0 * PI * KAmp);
  }

  //===== The simple way  =======================//
  // for (int i = 0; i < D; i++)
  //   if (fabs(OldMom[i] > Para.Kf))
  //     return 0.0;
  // return 1.0 / pow(2.0 * Para.Kf, D);
  //============================================//
}

double markov::ShiftK(const momentum &OldMom, momentum &NewMom)
{
  double x = Random.urn();
  double Prop;
  if (x < 1.0 / 3)
  {
    // COPYFROMTO(OldMom, NewMom);
    NewMom = OldMom;
    int dir = Random.irn(0, D - 1);
    double STEP = Para.Beta > 1.0 ? Para.Kf / Para.Beta * 3.0 : Para.Kf;
    NewMom[dir] += STEP * (Random.urn() - 0.5);
    Prop = 1.0;
  }
  else if (x < 2.0 / 3)
  {
    double k = OldMom.norm();
    if (k < 1.0e-9)
    {
      Prop = 0.0;
      NewMom = OldMom;
    }
    else
    {
      const double Lambda = 1.5; // Lambda can be larger for large-q
      double knew = k / Lambda + Random.urn() * (Lambda - 1.0 / Lambda) * k;
      double Ratio = knew / k;
      for (int i = 0; i < D; i++)
        NewMom[i] = OldMom[i] * Ratio;
      if (D == 2)
        Prop = 1.0;
      else if (D == 3)
        Prop = Ratio;
    }

    // if (isnan(Var.LoopMom[0][0]) || isnan(Var.LoopMom[0][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[1][0]) || isnan(Var.LoopMom[1][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[2][0]) || isnan(Var.LoopMom[2][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[3][0]) || isnan(Var.LoopMom[3][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[4][0]) || isnan(Var.LoopMom[4][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[5][0]) || isnan(Var.LoopMom[5][0])) {
    //   cout << "Na" << endl;
    // }
    // if (isnan(Var.LoopMom[6][0]) || isnan(Var.LoopMom[6][0])) {
    //   cout << "Na" << endl;
    // }
  }
  else
  {
    for (int i = 0; i < D; i++)
      NewMom[i] = -OldMom[i];
    Prop = 1.0;
  }

  return Prop;
};

double markov::ShiftExtK(const int &OldExtMomBin, int &NewExtMomBin)
{
  NewExtMomBin = Random.irn(0, ExtMomBinSize - 1);
  return 1.0;
};

double markov::ShiftTau(const double &OldTau, double &NewTau)
{
  double x = Random.urn();
  if (x < 1.0 / 3)
  {
    double DeltaT = Para.Beta / 3.0;
    NewTau = OldTau + DeltaT * (Random.urn() - 0.5);
  }
  else if (x < 2.0 / 3)
  {
    NewTau = -OldTau;
  }
  else
  {
    NewTau = Random.urn() * Para.Beta;
  }
  if (NewTau < 0.0)
    NewTau += Para.Beta;
  if (NewTau > Para.Beta)
    NewTau -= Para.Beta;
  return 1.0;
}

double markov::ShiftExtTau(const int &OldExtTauBin, int &NewExtTauBin)
{
  NewExtTauBin = Random.irn(0, ExtTauBinSize - 1);
  return 1.0;
};

std::string markov::_DetailBalanceStr(Updates op)
{
  string Output = string(80, '-') + "\n";
  Output += UpdatesName[op] + ":\n";
  double TotalProposed = 0.0, TotalAccepted = 0.0;
  for (int i = 0; i < Groups.size(); i++)
  {
    if (!Equal(Proposed[op][i], 0.0))
    {
      TotalAccepted += Accepted[op][i];
      TotalProposed += Proposed[op][i];
      Output +=
          fmt::sprintf("\t%8s:%15g%15g%15g\n", Groups[i].Name, Proposed[op][i],
                       Accepted[op][i], Accepted[op][i] / Proposed[op][i]);
      // fmt::format("\t%8s%4s:%15g%15g%15g\n", "Group", Groups[i].Name,
      //             Proposed[op][i], Accepted[op][i],
      //             Accepted[op][i] / Proposed[op][i]);
    }
  }
  if (!Equal(TotalProposed, 0.0))
  {
    Output += fmt::sprintf("\t%10s:%15g%15g%15g\n", "Summation", TotalProposed,
                           TotalAccepted, TotalAccepted / TotalProposed);
  }
  else
    Output += "\tNo updates are proposed/accepted!\n";
  return Output;
}

void markov::PrintMCInfo()
{
  string Output = "";
  Output = string(80, '=') + "\n";
  Output += "MC Counter: " + to_string(Para.Counter) + "\n";
  for (int i = 0; i < MCUpdates; i++)
    Output += _DetailBalanceStr((Updates)i);
  Output += string(80, '=') + "\n";
  LOG_INFO(Output);
}
