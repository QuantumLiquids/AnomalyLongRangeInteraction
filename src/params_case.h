/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2023/11/01.
 *
 */


#ifndef ANOMALY_LONG_RANGE_INTERACTION_SRC_PARAMS_CASE_H
#define ANOMALY_LONG_RANGE_INTERACTION_SRC_PARAMS_CASE_H

#include "gqmps2/case_params_parser.h"

using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Geometry = ParseStr("Geometry");
    L = ParseInt("L");
    omega_0 = ParseDouble("omega_0");
    omega_1 = ParseDouble("omega_1");
    J = ParseDouble("J");
    InteractionDecayLevel = ParseInt("InteractionDecayLevel");
    if (InteractionDecayLevel == 1) {
      J_exp_decay_length = ParseDouble("J_exp_decay_length");
    } else {
      J_exp_decay_length = 0;
    }
    if (InteractionDecayLevel == 0) {
      InteractionRange = ParseInt("InteractionRange");
    } else {
      InteractionRange = 0;
    }
    IsAnomaly = ParseBool("IsAnomaly");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    Threads = ParseInt("Threads");
    Perturbation = ParseDouble("Perturbation");
    wavelength = ParseDouble("wavelength");
    noise = ParseDoubleVec("noise");
  }

  std::string Geometry; // Cylinder, Torus, OBC, Rotated, Ladder (useless actually)
  size_t L;
  double omega_0;
  double omega_1;
  double J;
  int InteractionDecayLevel; // 0 don't decay; -1, 1/x decay; -2; 1/x^2 decay; 1 exponetially decay
  double J_exp_decay_length; //only for exponentially decay
  size_t InteractionRange;// only for const interaction.
  bool IsAnomaly;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t Threads;
  double Perturbation;
  double wavelength;
  std::vector<double> noise;
};

struct LongRangInteractionModelParamters {
 public :
  LongRangInteractionModelParamters(const double omega_0,
                                    const double omega_1,
                                    const double J,
                                    const int decay_level,
                                    const bool is_anomaly
  ) : omega_0(omega_0), omega_1(omega_1), J(J), decay_level(decay_level), is_anomaly(is_anomaly) {}

  LongRangInteractionModelParamters(const CaseParams &params
  ) : LongRangInteractionModelParamters(params.omega_0, params.omega_1, params.J, params.InteractionDecayLevel,
                                        params.IsAnomaly) {}

  inline void Print(void) {
    using std::cout;
    cout << " ===================****** Model Parameter List ******================= " << "\n";
    cout << "anomaly: ";
    if (is_anomaly) {
      cout << "OPEN.";
    } else {
      cout << "CLOSED.";
    }
    cout << "\n";
    cout << "omega_0  = " << omega_0 << "\n";
    cout << "omega_1 = " << omega_1 << "\n";
    cout << "J  = " << J << "\n";
    cout << "decay level = ";
    if (decay_level == 0) {
      cout << "const.";
    } else if (decay_level == -1) {
      cout << "1 / x decay.";
    } else if (decay_level == -2) {
      cout << "1 / x^2 decay.";
    } else if (decay_level == 1) {
      cout << "exponentially decay.";
    } else {
      cout << "Do not support.";
      exit(1);
    }
    cout << std::endl;
  }

  double omega_0;
  double omega_1;
  double J;
  int decay_level;
  bool is_anomaly;
};

#endif //ANOMALY_LONG_RANGE_INTERACTION_SRC_PARAMS_CASE_H
