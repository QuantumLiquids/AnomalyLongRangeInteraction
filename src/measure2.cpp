//
// Created by haoxinwang on 27/11/2023.
//


/**
 * 2 processes parallel
 */
#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <ctime>
#include "gqdouble.h"
//#include "myutil.h"
//#include "operators.h"
#include "params_case.h"
#include <random>
#include "dmrg_my_measure.h"


using namespace gqmps2;
using namespace gqten;
using namespace std;

int main(int argc, char *argv[]) {
//  namespace mpi = boost::mpi;
//  mpi::environment env;
//  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t L = params.L;
  LongRangInteractionModelParamters model_params(params);
  model_params.Print();

  clock_t startTime, endTime;
  startTime = clock();

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);


  Tensor sz = Tensor({pb_in, pb_out});
  Tensor sp = Tensor({pb_in, pb_out});
  Tensor sm = Tensor({pb_in, pb_out});
  Tensor id = Tensor({pb_in, pb_out});
  sz({0, 0}) = 0.5;
  sz({1, 1}) = -0.5;
  sp({0, 1}) = 1.0;
  sm({1, 0}) = 1.0;
  id({0, 0}) = 1.0;
  id({1, 1}) = 1.0;
  const SiteVec<TenElemT, U1QN> sites = SiteVec<TenElemT, U1QN>(L, pb_out);
  using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  std::vector<std::vector<size_t>> sites_set;
  for (size_t i = L / 4 + 1; i < std::min(3 * L / 4 + 2, L); i++) {
    sites_set.push_back({L / 4, i});
  }
  mps.Load();
  // Measure the S^2 operator

  Timer one_site_timer("measure  two site operators");
  MeasureTwoSiteOp(mps, {sz, sz}, id, sites_set, "zz");
  MeasureTwoSiteOp(mps, {sp, sm}, id, sites_set, "pm");
  MeasureTwoSiteOp(mps, {sm, sm}, id, sites_set, "mm");
  MeasureTwoSiteOp(mps, {sm, sp}, id, sites_set, "mp");
  MeasureTwoSiteOp(mps, {sp, sp}, id, sites_set, "pp");
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}