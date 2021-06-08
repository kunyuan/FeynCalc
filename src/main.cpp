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
#include "utility/fmt/format.h"
#include "utility/logger.h"
#include "utility/timer.h"
#include "weight.h"
#include <iostream>
#include <math.h>
#include <random>
#include <unistd.h>

using namespace std;
using namespace mc;
void InitPara();
void MonteCarlo();

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[]) {
  if (argc > 1) {
    Para.PID = atoi(argv[1]);
    Para.Seed = Para.PID;
  } else {
    std::random_device rd;
    Para.PID = rd() % 1000000;
    Para.Seed = Para.PID;
  }
  ifstream File;
  File.open("parameter", ios::in);
  cout << "Order, Beta, Rs, Mass2, Lambda, MaxExtMom(*kF), TotalStep(*1e6), "
          "Seed, "
          "PID\n";

  File >> Para.Order >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.Lambda >>
      Para.MinExtMom >> Para.MaxExtMom >> Para.TotalStep;
  File.close();
  InitPara(); // initialize global parameters
  MonteCarlo();
  return 0;
}

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  // Para.ObsType = FREQ;   //Polariation counterterm
  Para.ObsType = EQUALTIME;   //Self-energy counterterm

  Para.Type = POLAR;
  Para.SelfEnergyType = FOCK;
  // Para.SelfEnergyType = BARE;
  Para.UseVer4 = false;
  // Para.UseVer4 = true;

  if (Para.ObsType == FREQ) {
    Para.DiagFileFormat = "groups_charge/DiagPolar{}.txt";
    // Para.DiagFileFormat = "groups_spin/DiagPolar{}.txt";
    // Para.DiagFileFormat = "groups_spinless/DiagPolar{}.txt";
    Para.GroupName = {"0"}; // initialized with a normalization diagram
    Para.ReWeight = {1.0};
    for (int o = 1; o <= Para.Order; o++)
      for (int v = 0; v <= Para.Order - 1; v++) {
        // order 1 do not allow lambda counterterm
        if (o == 1 && v > 0)
          continue;
        for (int g = 0; g <= (Para.Order - 1) / 2; g++) {
          // if (g != 0)
          //   continue;
          // interaction+lambda counterterm+2*self-energy counter <=Order
          if (o + v + 2 * g > Para.Order)
            continue;
          auto name = to_string(o) + "_" + to_string(v) + "_" + to_string(g);
          cout << name << ", ";
          Para.GroupName.push_back(name);
          Para.ReWeight.push_back(pow(2.0, o));
        }
      }
    Para.GroupName.push_back("1_0_2");
    Para.ReWeight.push_back(10.0);
    Para.ReWeight[0] = Para.ReWeight[1] * 4.0;
  } else if (Para.ObsType == EQUALTIME) {
    Para.DiagFileFormat = "groups_mu/DiagPolar{}.txt";
    Para.GroupName = {
        "0",     "1_0_0", "1_0_1", "1_0_2", "2_1_0", "2_2_0",
        "2_3_0", "2_0_1", "2_1_1", "3_0_0", "3_1_0", "3_2_0",
        "3_0_1", "4_0_0", "4_1_0", "5_0_0"}; // initialized with a
                                             // normalization diagram
    Para.ReWeight = {4.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
                     2.0, 4.0, 4.0, 4.0, 4.0, 8.0, 4.0, 4.0};
  }

  ASSERT_ALLWAYS(Para.GroupName.size() == Para.ReWeight.size(),
                 "group num doesn't match reweight num!");
  ASSERT_ALLWAYS(Para.GroupName.size() < MaxGroupNum,
                 "MaxGroupNum should be larger!");

  // Para.ReWeight.clear();
  // ifstream File;
  // File.open("reweight.data");
  // if (File.is_open()) {
  //   for (int g = 0; g < Para.GroupName.size(); g++) {
  //     double weight;
  //     File >> weight;
  //     Para.ReWeight.push_back(weight);
  //   }
  // } else {
  //   for (int g = 0; g < Para.GroupName.size(); g++)
  //     Para.ReWeight.push_back(1.0);
  // }
  // Para.ReWeight.clear();
  // for (int g = 0; g < Para.GroupName.size(); g++)
  //   Para.ReWeight.push_back(1.0);

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
  Para.MinExtMom *= Kf;
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
                     << "Fermi Energy: " << Para.Ef << "\n"
                     << "Seed: " << Para.Seed << "\n");

  Para.PrinterTimer = 600;
  Para.SaveFileTimer = 600;
  Para.ReweightTimer = 600;
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

  // for (int Block = 0; Block < Para.TotalStep; Block++) {
  int Block = 0;
  while (Block < Para.TotalStep || Para.TotalStep <= 0) {
    Block++;
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
