#include "qldouble.h"

Tensor sz = Tensor({pb_in, pb_out});
Tensor sp = Tensor({pb_in, pb_out});
Tensor sm = Tensor({pb_in, pb_out});
Tensor id = Tensor({pb_in, pb_out});

Tensor sigma_z = Tensor({pb_in, pb_out});
Tensor sigma_x = Tensor({pb_in, pb_out});

Tensor sigma_z0 = Tensor({pb_in_front, pb_out_front});
Tensor sigma_z_end = Tensor({pb_in_end, pb_out_end});

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

    sigma_z0({0, 0}) = 1.0;
    sigma_z0({1, 1}) = -1.0;

    sigma_z_end({0, 0}) = 1.0;
    sigma_z_end({1, 1}) = -1.0;

    initialized = true;
  }
}
