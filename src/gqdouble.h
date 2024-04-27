#ifndef SPIN_ONE_HALF_SRC_GQDOUBLE_H
#define SPIN_ONE_HALF_SRC_GQDOUBLE_H

#include "boost/serialization/complex.hpp"
#include "qlten/qlten.h"

using qlten::QNCard;
using qlten::U1QNVal;
using qlten::TenIndexDirType;

//using TenElemT = qlten::QLTEN_Complex;
using TenElemT = qlten::QLTEN_Double;
using U1QN = qlten::special_qn::U1QN;
using Tensor = qlten::QLTensor<TenElemT, U1QN>;

using QNSctT = qlten::QNSector<U1QN>;
using IndexT = qlten::Index<U1QN>;

const auto qn0 = U1QN(0);

const IndexT pb_out = IndexT({QNSctT(U1QN(0), 2)},
                             TenIndexDirType::OUT
);

const auto pb_in = qlten::InverseIndex(pb_out);

#endif // SPIN_ONE_HALF_SRC_GQDOUBLE_H
