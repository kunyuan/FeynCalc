//
//  main.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include "global.h"
#include "markov.h"
#include "utility/abort.h"
#include "utility/logger.h"
#include "utility/timer.h"
#include "weight.h"
#include <iostream>
#include <math.h>
#include <unistd.h>

using namespace std;
using namespace mc;
void InitPara();
void MonteCarlo();

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[]) {
  cout << "Order, Beta, Rs, Mass2, MaxExtMom(*kF), TotalStep(*1e6), Seed, "
          "PID\n";
  cin >> Para.Order >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.MaxExtMom >>
      Para.TotalStep >> Para.Seed >> Para.PID;
  InitPara(); // initialize global parameters
  MonteCarlo();
  return 0;
}

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  Para.Type = POLAR;
  Para.ObsType = FREQ;

  Para.UseVer4 = false;
  // Para.UseVer4 = true;

  // diagram file path: groups/DiagPolar1.dat
  Para.DiagFileFormat = "groups/DiagPolar{}.txt";
  Para.GroupName.clear();
  for (int o = 1; o <= Para.Order; o++)
    for (int v = 0; v <= Para.Order - 1; v++) {
      // order 1 do not allow lambda counterterm
      if (o == 1 && v > 0)
        continue;
      for (int g = 0; g <= (Para.Order - 1) / 2; g++) {
        // cout << o << ", " << v << ", " << g << ", " << o + v + 2 * g << endl;
        // interaction+lambda counterterm+2*self-energy counter <=Order
        if (o + v + 2 * g > Para.Order)
          continue;
        auto name = to_string(o) + "_" + to_string(v) + "_" + to_string(g);
        cout << name << ", ";
        Para.GroupName.push_back(name);
      }
    }
  cout << endl;
  // Para.GroupName = {"1", "2", "3"};

  // Para.GroupName = {"1", "2"};
  Para.ReWeight = {1.0, 1.0, 5.0, 10.0, 0.1};
  // Para.SelfEnergyType = FOCK;
  Para.SelfEnergyType = BARE;

  //// initialize the global parameter //////////////////////
  double Kf;
  if (D == 3) {
    Kf = pow(9.0 * PI / 4.0, 1.0 / 3.0) / Para.Rs; // 3D
  } else if (D == 2) {
    Kf = sqrt(2.0) / Para.Rs; // 2D
  } else {
    ABORT("Dimension " << D << " has not yet been implemented!");
  }
  Para.Kf = Kf;
  Para.Ef = Kf * Kf;
  Para.Mu = Para.Ef;
  Para.MaxExtMom *= Kf;

  // scale all energy with E_F
  Para.Beta /= Para.Ef;
  Para.UVScale = 8.0 * Para.Ef;
  Para.UVCoupling = 1.0 * Para.Ef;

  LOG_INFO("Order: " << Para.Order << "\nInverse Temperature: " << Para.Beta
                     << "\nUV Energy Scale: " << Para.UVScale
                     << "\nUV Coupling: " << Para.UVCoupling << "\n"
                     << "r_s: " << Para.Rs << "\n"
                     << "Fermi Mom: " << Para.Kf << "\n"
                     << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 5;
  Para.SaveFileTimer = 30;
  Para.ReweightTimer = 30;
}

void MonteCarlo() {
  LOG_INFO("Initializing Markov!");
  markov Markov;
  InterruptHandler Interrupt;

  Random.Reset(Para.Seed);
  Para.Counter = 0;

  timer ReweightTimer, PrinterTimer, SaveFileTimer, MessageTimer;
  PrinterTimer.start();
  SaveFileTimer.start();
  MessageTimer.start();
  ReweightTimer.start();

  LOG_INFO("Start simulation ...")

  for (int Block = 0; Block < Para.TotalStep; Block++) {
    for (int i = 0; i < 1000000; i++) {
      Para.Counter++;
      // if (Para.Counter == 9) {
      //   cout << "Before: " << Para.Counter << endl;
      //   PrintDeBugMCInfo();
      // }

      double x = Random.urn();
      if (x < 1.0 / 3.0) {
        Markov.ChangeGroup();
        // ;
      } else if (x < 2.0 / 3.0) {
        Markov.ChangeMomentum();
        // ;
      } else if (x < 3.0 / 3.0) {
        Markov.ChangeTau();
        // ;
      }

      // if (Para.Counter == 8831001) {
      //   cout << "After: " << Para.Counter << endl;
      //   PrintDeBugMCInfo();
      // }

      Markov.Measure();
      // Markov.DynamicTest();

      if (i % 1000 == 0) {
        // Markov.PrintDeBugMCInfo();
        if (PrinterTimer.check(Para.PrinterTimer)) {
          Markov.DynamicTest();
          Markov.PrintDeBugMCInfo();
          Markov.PrintMCInfo();
          LOG_INFO(ProgressBar((double)Block / Para.TotalStep));
        }

        if (SaveFileTimer.check(Para.SaveFileTimer)) {
          Interrupt.Delay(); // the process can not be killed in saving
          Markov.SaveToFile();
          Interrupt.Resume(); // after this point, the process can be killed
        }

        if (ReweightTimer.check(Para.ReweightTimer)) {
          Markov.AdjustGroupReWeight();
          Para.ReweightTimer *= 1.5;
        }
      }
    }
  }

  LOG_INFO("Simulation is ended!");
  Markov.PrintMCInfo();
  Interrupt.Delay(); // the process can not be killed in saving
  Markov.SaveToFile();
  Interrupt.Resume(); // after this point, the process can be killed
  LOG_INFO("Quit Markov.");
}