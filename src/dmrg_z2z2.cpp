#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <stdlib.h>
#include "qldouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include <random>

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  int boundary_condition(std::atoi(argv[2]));
  if (boundary_condition < 0 || boundary_condition > 3) {
    std::cerr << "boundary condition error!" << std::endl;
  }

  size_t L = params.L;
  LongRangInteractionModelParamters model_params(params);
  model_params.Print();

  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  std::vector<size_t> Dmin_set, Dmax_set;
  if (has_bond_dimension_parameter) {
    Dmin_set = input_D_set;
    Dmax_set = input_D_set;
  } else {
    Dmin_set = {params.Dmin};
    Dmax_set = {params.Dmax};
  }
  size_t DMRG_time = Dmax_set.size();
  std::vector<size_t> MaxLanczIterSet(DMRG_time);
  MaxLanczIterSet.back() = params.MaxLanczIter;
  if (DMRG_time > 1) {
    size_t MaxLanczIterSetSpace;
    MaxLanczIterSet[0] = 3;
    MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
    std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
    for (size_t i = 1; i < DMRG_time - 1; i++) {
      MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
      std::cout << MaxLanczIterSet[i] << ", ";
    }
    std::cout << MaxLanczIterSet.back() << "]" << std::endl;
  } else {
    std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
  }

  clock_t startTime, endTime;
  startTime = clock();
  OperatorInitial();
  using QNT = Z2Z2QN;
  IndexVec<QNT> phy_indices(L);
  phy_indices[0] = pb_out_front;
  phy_indices.back() = pb_out_end;
  for (size_t i = 1; i < L - 1; i++) {
    phy_indices[i] = pb_out;
  }
  const SiteVec<TenElemT, QNT> sites = SiteVec<TenElemT, QNT>(phy_indices);

//  qlmps::MPO<Tensor> mpo(L);
//  const std::string kMpoPath = "mpo";
//  const std::string kMpoTenBaseName = "mpo_ten";
//  if (IsPathExist(kMpoPath)) {
//    for (size_t i = 0; i < mpo.size(); i++) {
//      std::string filename = kMpoPath + "/" +
//                             kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
//      mpo.LoadTen(i, filename);
//    }
//    cout << "MPO loaded." << endl;
//  } else {
//    cout << "No mpo directory. start to generate mpo" << std::endl;
  qlmps::MPOGenerator<TenElemT, Z2Z2QN> mpo_gen(sites, qn0);
  if (model_params.is_anomaly) {
    for (size_t i = 0; i < L - 2; i++) {
      if (i == 0) {
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_z0, i, sigma_x, i + 1);
        mpo_gen.AddTerm(-0.5 * model_params.omega_1, {sigma_z0, sigma_x, sigma_z}, {i, i + 1, i + 2});
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_x, i + 1, sigma_z, i + 2);
      } else if (i == L - 3) {
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_z, i, sigma_x, i + 1);
        mpo_gen.AddTerm(-0.5 * model_params.omega_1, {sigma_z, sigma_x, sigma_z_end}, {i, i + 1, i + 2});
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_x, i + 1, sigma_z_end, i + 2);
      } else {
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_z, i, sigma_x, i + 1);
        mpo_gen.AddTerm(-0.5 * model_params.omega_1, {sigma_z, sigma_x, sigma_z}, {i, i + 1, i + 2});
        mpo_gen.AddTerm(0.5 * model_params.omega_0, sigma_x, i + 1, sigma_z, i + 2);
      }
      mpo_gen.AddTerm(0.5 * model_params.omega_1, sigma_x, i + 1);
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
        Tensor op1, op2;
        if (i == 0) {
          op1 = sigma_z0;
        } else {
          op1 = sigma_z;
        }
        if (j == L - 1) {
          op2 = sigma_z_end;
        } else {
          op2 = sigma_z;
        }
        mpo_gen.AddTerm(model_params.J / 2, op1, i, op2, j);
      }
    }
  } else if (model_params.decay_level == -1) {
    for (size_t i = 0; i < L; i++) {
      for (size_t j = i + 2; j < L; j++) {
        Tensor op1, op2;
        if (i == 0) {
          op1 = sigma_z0;
        } else {
          op1 = sigma_z;
        }
        if (j == L - 1) {
          op2 = sigma_z_end;
        } else {
          op2 = sigma_z;
        }
        mpo_gen.AddTerm(model_params.J / 2 / (j - i), op1, i, op2, j);
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

  auto mpo = mpo_gen.GenMatReprMPO(true);
//  }

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
  FiniteMPST mps(sites);

  if (params.Threads == 0) {
    const size_t max_threads = std::thread::hardware_concurrency();
    params.Threads = max_threads;
  }

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  std::vector<long unsigned int> stat_labs(L);
  for (size_t i = 0; i < L; i++) {
    stat_labs[i] = (i % 2);
  }

  std::srand(std::time(nullptr));
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(stat_labs.begin(), stat_labs.end(), g);
  stat_labs[0] = boundary_condition % 2;
  stat_labs[L - 1] = (boundary_condition / 2) % 2;
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
  for (size_t i = 0; i < DMRG_time; i++) {
    std::cout << "D_max = " << Dmax_set[i] << std::endl;
    qlmps::FiniteVMPSSweepParams sweep_params(
        params.Sweeps,
        Dmin_set[i], Dmax_set[i], params.CutOff,
        qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
    );
    e0 = qlmps::FiniteDMRG(mps, mpo, sweep_params);
    mps_path = sweep_params.mps_path;
  }

  std::cout << "E0/site: " << e0 / L << std::endl;
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

  mps.Load(mps_path);
  auto ee_list = mps.GetEntanglementEntropy(1);
  std::copy(ee_list.begin(), ee_list.end(), std::ostream_iterator<double>(std::cout, " "));

  std::string file_name = "ee";
  std::ofstream ofs(file_name, std::ofstream::binary);
  ofs.write((const char *) ee_list.data(), ee_list.size() * sizeof(double));
  ofs << std::endl;
  ofs.close();

  return 0;
}


