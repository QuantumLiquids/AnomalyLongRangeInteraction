#include "gqdouble.h"

Tensor sz = Tensor({pb_in, pb_out});
Tensor sp = Tensor({pb_in, pb_out});
Tensor sm = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});

Tensor sigma_z = Tensor({pb_in, pb_out});
Tensor sigma_x = Tensor({pb_in, pb_out});

void OperatorInitial() {
  static bool initialized = false;
  if (!initialized) {
    sz({0, 0}) = 0.5;
    sz({1, 1}) = -0.5;
    sp({0, 1}) = 1.0;
    sm({1, 0}) = 1.0;
    id({0, 0}) = 1;
    id({1, 1}) = 1;

    sigma_z = 2.0 * sz;
    sigma_x = sp + sm;

    initialized = true;
  }
}
