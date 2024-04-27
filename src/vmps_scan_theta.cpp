#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <stdlib.h>
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include <random>
#include <iostream>
#include <sstream>
#include <filesystem>


using namespace qlmps;
using namespace qlten;
using namespace std;

void DMRGForOneTheta(const CaseParams &params,
                     const double theta);

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  if (params.Threads == 0) {
    const size_t max_threads = std::thread::hardware_concurrency();
    params.Threads = max_threads;
  }

  double start_theta, end_theta;
  if (argc == 2) {
    start_theta = 0;
    end_theta = M_PI / 2;
  } else if (argc == 3) {
    start_theta = std::atof(argv[2]);
    end_theta = M_PI / 2;
  } else if (argc == 4) {
    start_theta = std::atof(argv[2]);
    end_theta = std::atof(argv[3]);
  }

  qlten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  std::vector<double> theta_set = linspace(0, M_PI / 2, 30);

  theta_set.erase(std::remove_if(theta_set.begin(), theta_set.end(),
                                 [start_theta, end_theta](const double &theta) {
                                   return theta < start_theta || theta >= end_theta;
                                 }),
                  theta_set.end());

  for (double theta: theta_set) {
    DMRGForOneTheta(params, theta);
  }
  return 0;
}


std::string GenerateEEFilename(const CaseParams &params,
                               const LongRangInteractionModelParamters &model_params) {
  std::ostringstream filename;
  filename << "eeLRIL" << params.L
           << "omega0" << std::fixed << std::setprecision(4) << model_params.omega_0
           << "omega1" << std::fixed << std::setprecision(4) << model_params.omega_1
           << "J" << std::fixed << std::setprecision(1) << model_params.J
           << "Level" << model_params.decay_level
           << "P" << params.InteractionRange
           << "D" << params.Dmax;
  return filename.str();
}

std::string GenerateEEFilenameTheta0(const CaseParams &params,
                                     const LongRangInteractionModelParamters &model_params) {
  std::ostringstream filename;
  filename << "eeLRIL" << params.L
           << "omega0" << std::fixed << std::setprecision(1) << model_params.omega_0
           << "omega1" << std::fixed << std::setprecision(1) << model_params.omega_1
           << "J" << std::fixed << std::setprecision(1) << model_params.J
           << "Level" << model_params.decay_level
           << "P" << params.InteractionRange
           << "D" << params.Dmax;
  return filename.str();
}


std::string GenerateMPSDirectroyName(const CaseParams &params,
                                     const LongRangInteractionModelParamters &model_params) {
  std::ostringstream dirname;
  dirname << "mpsLRIL" << params.L
          << "omega0" << std::fixed << std::setprecision(4) << model_params.omega_0
          << "omega1" << std::fixed << std::setprecision(4) << model_params.omega_1
          << "J" << std::fixed << std::setprecision(1) << model_params.J
          << "Level" << model_params.decay_level
          << "R" << params.InteractionRange
          << "D" << params.Dmax;
  return dirname.str();
}

void DMRGForOneTheta(const CaseParams &params,
                     const double theta) {
  size_t L = params.L;
  LongRangInteractionModelParamters model_params(params);
  model_params.omega_0 = std::cos(theta);
  model_params.omega_1 = std::sin(theta);

//  model_params.omega_0 = std::round(model_params.omega_0 * 10000.0) / 10000.0;
//  model_params.omega_1 = std::round(model_params.omega_1 * 10000.0) / 10000.0;

  model_params.Print();

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  using QNT = U1QN;
  const SiteVec<TenElemT, QNT> sites = SiteVec<TenElemT, QNT>(L, pb_out);


  qlmps::MPO<Tensor> mpo(L);
  qlmps::MPOGenerator<TenElemT, U1QN> mpo_gen(sites, qn0);
  if (model_params.is_anomaly) {
    for (size_t i = 0; i < L - 2; i++) {
      mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_z, i, sigma_x, i + 1);
      mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_x, i + 1, sigma_z, i + 2);
      mpo_gen.AddTerm(0.5 * model_params.omega_1, sigma_x, i + 1);
      mpo_gen.AddTerm(-0.5 * model_params.omega_1, {sigma_z, sigma_x, sigma_z}, {i, i + 1, i + 2});
    }
  } else {
    for (size_t i = 0; i < L - 2; i++) {
      mpo_gen.AddTerm(0.5 * (model_params.omega_1 + model_params.omega_0), sigma_x, i + 1);
      mpo_gen.AddTerm(-0.5 * (model_params.omega_0 - model_params.omega_1), {sigma_z, sigma_x, sigma_z},
                      {i, i + 1, i + 2});
    }
  }

  if (model_params.decay_level == 0) { //const
    for (size_t i = 0; i < L; i++) {
      for (size_t j = i + 2; j < std::min(i + params.InteractionRange + 2, L); j++) {
        mpo_gen.AddTerm(model_params.J / 2, sigma_z, i, sigma_z, j);
      }
    }
  } else if (model_params.decay_level == -1) {
    for (size_t i = 0; i < L; i++) {
      for (size_t j = i + 2; j < L; j++) {
        mpo_gen.AddTerm(model_params.J / 2 / (j - i), sigma_z, i, sigma_z, j);
      }
    }
  } else if (model_params.decay_level == -2) {
    for (size_t i = 0; i < L; i++) {
      for (size_t j = i + 2; j < L; j++) {
        mpo_gen.AddTerm(model_params.J / 2 / (j - i) / (j - i), sigma_z, i, sigma_z, j);
      }
    }
  } else if (model_params.decay_level == 1) {
    for (size_t i = 0; i < L; i++) {
      for (size_t j = i + 2; j < L; j++) {
        mpo_gen.AddTerm(model_params.J * std::exp(-(double) (j - i) / params.J_exp_decay_length),
                        sigma_z, i, sigma_z, j);
      }
    }
  } else {
    std::cout << "Do not support." << std::endl;
    exit(1);
  }

  mpo = mpo_gen.Gen();

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(L);
  for (size_t i = 0; i < L; i++) {
    stat_labs[i] = (i % 2);
  }

  std::srand(std::time(nullptr));
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(stat_labs.begin(), stat_labs.end(), g);

  if (IsPathExist(kMpsPath)) {
    if (L == GetNumofMps()) {
      cout << "The number of mps files is consistent with mps size." << endl;
      cout << "Directly use mps from files." << endl;
    } else {
      qlmps::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(kMpsPath, true);
    }
  } else {
    cout << "No mps directory." << endl;
    qlmps::DirectStateInitMps(mps, stat_labs);
    cout << "Initial mps as direct product state." << endl;
    mps.Dump(kMpsPath, true);
  }
  double e0;

  std::string mps_path;

  size_t sweeps;
  if (theta > 1e-15) {
    sweeps = 3;
  } else {
    sweeps = params.Sweeps;
  }
  qlmps::FiniteVMPSSweepParams sweep_params(
      sweeps,
      params.Dmin, params.Dmax,
      params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );
  e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
  mps_path = sweep_params.mps_path;

  std::cout << "E0/site: " << e0 / L << std::endl;
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  mps.Load(mps_path);
  auto ee_list = mps.GetEntanglementEntropy(1);
  std::copy(ee_list.begin(), ee_list.end(), std::ostream_iterator<double>(std::cout, " "));

  std::string ee_file_name = GenerateEEFilename(params, model_params);
  std::ofstream ofs(ee_file_name, std::ofstream::binary);
  ofs.write((const char *) ee_list.data(), ee_list.size() * sizeof(double));
  ofs << std::endl;
  ofs.close();

  if (std::abs(theta) < 1e-15) {
    std::string ee_file_name = GenerateEEFilenameTheta0(params, model_params);
    std::ofstream ofs(ee_file_name, std::ofstream::binary);
    ofs.write((const char *) ee_list.data(), ee_list.size() * sizeof(double));
    ofs << std::endl;
    ofs.close();
  }


  std::string mps_dir_name = GenerateMPSDirectroyName(params, model_params);
  std::filesystem::copy(sweep_params.mps_path, mps_dir_name);
  std::filesystem::remove_all(sweep_params.temp_path);
}