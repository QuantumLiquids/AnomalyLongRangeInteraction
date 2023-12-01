#ifndef SPIN_ONE_HALF_SRC_GQDOUBLE_H
#define SPIN_ONE_HALF_SRC_GQDOUBLE_H

#include "boost/serialization/complex.hpp"
#include "gqten/gqten.h"

using gqten::QNCard;
using gqten::U1QNVal;
using gqten::GQTenIndexDirType;

//using TenElemT = gqten::GQTEN_Complex;
using TenElemT = gqten::GQTEN_Double;
using U1QN = gqten::special_qn::U1QN;
using Tensor = gqten::GQTensor<TenElemT, U1QN>;

using QNSctT = gqten::QNSector<U1QN>;
using IndexT = gqten::Index<U1QN>;

const auto qn0 = U1QN(0);

const IndexT pb_out = IndexT({QNSctT(U1QN(0), 2)},
                             GQTenIndexDirType::OUT
);

const auto pb_in = gqten::InverseIndex(pb_out);

#endif // SPIN_ONE_HALF_SRC_GQDOUBLE_H
