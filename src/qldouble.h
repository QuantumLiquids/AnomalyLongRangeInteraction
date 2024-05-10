#ifndef SPIN_ONE_HALF_SRC_GQDOUBLE_H
#define SPIN_ONE_HALF_SRC_GQDOUBLE_H

#include "boost/serialization/complex.hpp"
#include "qlten/qlten.h"
#include "z2z2qn.h"

using qlten::QNCard;
using qlten::TenIndexDirType;

//using TenElemT = qlten::QLTEN_Complex;
using TenElemT = qlten::QLTEN_Double;
using qlten::special_qn::Z2Z2QN;
using Tensor = qlten::QLTensor<TenElemT, Z2Z2QN>;

using QNSctT = qlten::QNSector<Z2Z2QN>;
using IndexT = qlten::Index<Z2Z2QN>;

const auto qn0 = Z2Z2QN(0, 0);

const IndexT pb_out = IndexT({QNSctT(Z2Z2QN(0, 0), 2)},
                             TenIndexDirType::OUT
);

const auto pb_in = qlten::InverseIndex(pb_out);

const IndexT pb_out_front = IndexT({QNSctT(Z2Z2QN(1, 0), 1),
                                    QNSctT(Z2Z2QN(0, 0), 1)},
                                   TenIndexDirType::OUT
);
const auto pb_in_front = qlten::InverseIndex(pb_out_front);

const IndexT pb_out_end = IndexT({QNSctT(Z2Z2QN(0, 1), 1),
                                  QNSctT(Z2Z2QN(0, 0), 1)},
                                 TenIndexDirType::OUT);
const auto pb_in_end = qlten::InverseIndex(pb_out_end);

#endif // SPIN_ONE_HALF_SRC_GQDOUBLE_H
