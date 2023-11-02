#include "gqmps2/gqmps2.h"
#include "myutil.h"
#include <random>

using gqmps2::kMpsPath;
using gqmps2::kMpsTenBaseName;
using gqmps2::kGQTenFileSuffix;

//number of mps file in default mps path("./mps")
size_t GetNumofMps() {
  size_t NumberOfMpsFile = 0;
  for (NumberOfMpsFile = 0; NumberOfMpsFile < 1e5; NumberOfMpsFile++) {
    std::string file;
    file = kMpsPath + "/" + kMpsTenBaseName + std::to_string(NumberOfMpsFile) + "." + kGQTenFileSuffix;
    std::ifstream ifs(file, std::ifstream::binary);
    if (ifs.good()) {
      ifs.close();
    } else {
      break;
    }
  }
  return NumberOfMpsFile;
}

void Show(std::vector<size_t> v) {
  for (auto iter = v.begin(); iter < v.end(); iter++) {
    std::cout << *iter << ",";
  }
  std::cout << '\b' << std::endl;
}

bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t> &D_set) {
  int nOptionIndex = 1;
  std::string D_string;
  std::string arguement1 = "--D=";
  bool has_D_parameter(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      D_string = &argv[nOptionIndex][arguement1.size()];
      has_D_parameter = true;
    }
    nOptionIndex++;
  }

  //split thread num list
  const char *split = ",";
  char *p;
  const size_t MAX_CHAR_LENTH = 1000;
  char D_char[MAX_CHAR_LENTH];
  for (size_t i = 0; i < MAX_CHAR_LENTH; i++) {
    D_char[i] = 0;
  }

  strcpy(D_char, D_string.c_str());

  p = strtok(D_char, split);
  while (p != nullptr) {
    D_set.push_back(atoi(p));
    p = strtok(nullptr, split);
  }

  return has_D_parameter;
}

/// for translation symmetry
std::vector<size_t> GenerateDirectStateLabel(const size_t ly,
                                             const size_t lx,
                                             const size_t num_hole,
                                             const size_t ky_int) {
  const size_t hole_lab = 3;
  const size_t double_occupy_lab = 0;
  size_t n = lx * ly;
  size_t num_ele = n - num_hole;

  assert(ky_int < ly);
  assert(num_hole < n);
  std::vector<size_t> state_label(n, hole_lab); // all holes;
  std::vector<size_t> filled_sites(num_ele);
  size_t total_ky_int = 0;

  std::srand(std::time(nullptr));
  std::random_device rd;
  std::mt19937 g(rd());

  if (num_hole >= ly) {
    for (size_t i = 0; i < num_ele - 1; i++) {
      filled_sites[i] = i;
      total_ky_int += i;
    }
    total_ky_int = total_ky_int % ly;
    size_t remain_ky_int = (ky_int - total_ky_int + ly) % ly;
    filled_sites.back() = remain_ky_int + (lx - 1) * ly;
    total_ky_int += remain_ky_int;
    assert(total_ky_int % ly == ky_int);

    for (size_t i = 0; i < num_ele; i++) {
      size_t site = filled_sites[i];
      state_label[site] = i % 2 + 1; //spin up/down
    }
  } else if (num_hole == 0) { //half filling
    for (size_t i = 0; i < n; i++) {
      state_label[i] = i % 2 + 1; //spin up/down
      total_ky_int += i;
    }
    total_ky_int = total_ky_int % ly;
    if (total_ky_int != ky_int) {
      size_t remain_ky_int = (ky_int - total_ky_int + ly) % ly;
      state_label[0] = hole_lab;
      state_label[remain_ky_int] = double_occupy_lab;
    }
  } else if (num_hole <= lx) { // 0 < num_hole < ly
    size_t sz_lable = 0;
    for (size_t i = 0; i < n; i++) {
      size_t y = i % ly;
      size_t x = i / ly;
      if (y == 0 && x < num_hole) {
        state_label[i] = hole_lab;
      } else {
        state_label[i] = sz_lable % 2 + 1;
        total_ky_int += i;
      }
    }
    total_ky_int = total_ky_int % ly;
    size_t over_ky_int = (total_ky_int + ly - ky_int) % ly;
    std::swap(state_label[0], state_label[over_ky_int]);
  } else { //lx < num_hole < ly
    size_t sz_lable = 0;
    for (size_t i = 0; i < n; i++) {
      size_t y = i % ly;
      size_t x = i / ly;
      if (y == 0) {
        state_label[i] = hole_lab;
      } else if (x == 0 && y < (num_hole - lx + 1)) {
        state_label[i] = hole_lab;
      } else {
        state_label[i] = sz_lable % 2 + 1;
        total_ky_int += i;
      }
    }
    total_ky_int = total_ky_int % ly;
    size_t over_ky_int = (total_ky_int + ly - ky_int) % ly;
    std::swap(state_label[ly], state_label[over_ky_int + ly]);
  }
  if (num_hole != 0) {
    for (size_t y = 0; y < ly; y++) {
      std::vector<size_t> labs_on_leg(lx);
      for (size_t x = 0; x < lx; x++) {
        labs_on_leg[x] = state_label[x * ly + y];
      }
      std::shuffle(labs_on_leg.begin(), labs_on_leg.end(), g);
      for (size_t x = 0; x < lx; x++) {
        state_label[x * ly + y] = labs_on_leg[x];
      }
    }
  }
  std::copy(state_label.begin(), state_label.end(), std::ostream_iterator<int>(std::cout, " "));
  return state_label;
}

/**
 * Generate all the pairs of momenta (k1,k2,k3,k4)
 * which satisfy (k1 + k2 - k3 - k4) % ly
 * and 0 <= k1 < k2 <= ly - 1, 0 <= k3 < k4 <= ly - 1.
 *
 * This function can used in measuring the four-fermion correlation functions.
 *
 * @param ly
 * @return
 */
std::vector<std::vector<size_t>> GenAllOrderedMomentumPairs(const size_t ly) {
  std::vector<std::vector<size_t>> momentum_pairs;
  for (size_t k1 = 0; k1 < ly - 1; k1++) {
    for (size_t k2 = k1 + 1; k2 < ly; k2++) {
      for (size_t k3 = 0; k3 < ly - 1; k3++) {
        for (size_t k4 = k3 + 1; k4 < ly; k4++) {
          if ((k1 + k2 - k3 - k4 + 2 * ly) % ly == 0) {
            momentum_pairs.push_back({k1, k2, k3, k4});
          }
        }
      }
    }
  }
  return momentum_pairs;
}


std::vector<std::vector<size_t>> GenAllMomentumPairs(const size_t ly) {
  std::vector<std::vector<size_t>> momentum_pairs;
  for (size_t k1 = 0; k1 < ly ; k1++) {
    for (size_t k2 = 0; k2 < ly; k2++) {
      for (size_t k3 = 0; k3 < ly ; k3++) {
        for (size_t k4 = 0; k4 < ly; k4++) {
          if ((k1 + k2 - k3 - k4 + 2 * ly) % ly == 0) {
            momentum_pairs.push_back({k1, k2, k3, k4});
          }
        }
      }
    }
  }
  return momentum_pairs;
}


/**
 * Generate all the pairs of momenta (k1,k2,k3)
 * which satisfy (2 * k1 - k2 - k3 ) % ly
 * and 0 <= k1 <= ly - 1, 0 <= k2 < k3 <= ly - 1.
 */
std::vector<std::vector<size_t>> GenAllMomentumPairs3(const size_t ly) {
  std::vector<std::vector<size_t>> momentum_pairs;
  for (size_t k1 = 0; k1 < ly; k1++) {
    for (size_t k3 = 0; k3 < ly - 1; k3++) {
      for (size_t k4 = k3 + 1; k4 < ly; k4++) {
        if ((2 * k1 - k3 - k4 + 2 * ly) % ly == 0) {
          momentum_pairs.push_back({k1, k3, k4});
        }
      }
    }
  }
  return momentum_pairs;
}

bool ParserMeasureSite(const int argc, char *argv[],
            size_t &start,
            size_t &end) {
  int nOptionIndex = 1;

  std::string arguement1 = "--start=";
  std::string arguement2 = "--end=";
  bool start_argument_has(false), end_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      start = atoi(para_string.c_str());
      start_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement2.size()];
      end = atoi(para_string.c_str());
      end_argument_has = true;
    }
    nOptionIndex++;
  }

  if (start_argument_has != end_argument_has) {
    std::cout << "Only setting one start/end argument, exit(1)." << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    std::cout << "Note: no start/end argument, set it by default (L/4, 3*L/4+2)." << std::endl;
  }

  return start_argument_has;
}