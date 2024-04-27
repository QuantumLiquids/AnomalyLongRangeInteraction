//
// Created by haoxinwang on 27/11/2023.
//


/**
 * 2 processes parallel
 */
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "gqdouble.h"
//#include "myutil.h"
//#include "operators.h"
#include "params_case.h"
#include <random>
#include "dmrg_my_measure.h"



using namespace qlmps;
using namespace qlten;
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

  qlten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);


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
  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1QN>;
  FiniteMPST mps(sites);

  Timer one_site_timer("measure  one site operators");
  MeasureOneSiteOp(mps, kMpsPath, sz, "sz");
  MeasureOneSiteOp(mps, kMpsPath, sp, "sp");
  MeasureOneSiteOp(mps, kMpsPath, sm, "sm");
  cout << "measured one point function.<====" << endl;
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}