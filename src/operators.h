
#ifndef SPIN_ONE_HALF_SRC_OPERATORS_H
#define SPIN_ONE_HALF_SRC_OPERATORS_H

#include "qlten/qlten.h"
#include "qldouble.h"

extern Tensor sz, sp, sm, id;
extern Tensor sigma_z, sigma_x;

extern Tensor sigma_z0, sigma_z_end;
void OperatorInitial();

#endif // SPIN_ONE_HALF_SRC_OPERATORS_H