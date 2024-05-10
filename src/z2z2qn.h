// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <haoxinwang@cuhk.edu.hk>
* Creation Date: 2024-May-10
*
* Description: QuantumLiquids/tensor project. `Z2 \cross Z2` QN class.
*/

/**
@file z2z2qn.h
*/

#ifndef QLTEN_QLTENSOR_SPECIAL_QN_Z2Z2QN_H
#define QLTEN_QLTENSOR_SPECIAL_QN_Z2Z2QN_H

#include "qlten/framework/vec_hash.h"   //_HASH_XXPRIME_1...
#include "qlten/qltensor/qn/qn.h"       //QNCardVec
#include "qlten/qltensor/qn/qnval.h"    //QNVal
#include "qlten/qltensor/qn/qnval_u1.h" //Z2Z2QNVal
#include <boost/serialization/serialization.hpp>

namespace qlten {
namespace special_qn {

class Z2Z2QN : public Showable {
 public:
  Z2Z2QN(void);
  Z2Z2QN(const int z2val1, const int z2val2);
  Z2Z2QN(const Z2Z2QN &);

  Z2Z2QN &operator=(const Z2Z2QN &);

  ~Z2Z2QN(void);

  Z2Z2QN operator-(void) const;
  Z2Z2QN &operator+=(const Z2Z2QN &);

  Z2Z2QN operator+(const Z2Z2QN &rhs) const;
  Z2Z2QN operator-(const Z2Z2QN &rhs) const;

  size_t dim(void) const { return 1; }

  bool operator==(const Z2Z2QN &rhs) const {
    return hash_ == rhs.hash_;
  }

  bool operator!=(const Z2Z2QN &rhs) const {
    return !(*this == rhs);
  }

  //Hashable
  size_t Hash() const { return hash_; }

  void StreamRead(std::istream &);
  void StreamWrite(std::ostream &) const;

  void Show(const size_t indent_level = 0) const override;

  static Z2Z2QN Zero(void) {
    return Z2Z2QN(0, 0);
  }

 private:
  size_t CalcHash_(void) const;

  int z2vals_[2];
  size_t hash_;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & z2vals_;
    ar & hash_;
  }
};

inline Z2Z2QN::Z2Z2QN(void) : z2vals_{0, 0}, hash_(CalcHash_()) {}

inline Z2Z2QN::Z2Z2QN(const int z2val1, const int z2val2) :
    z2vals_{z2val1, z2val2}, hash_(CalcHash_()) {}

inline Z2Z2QN::Z2Z2QN(const Z2Z2QN &rhs) :
    z2vals_{rhs.z2vals_[0], rhs.z2vals_[1]}, hash_(rhs.hash_) {}

inline Z2Z2QN::~Z2Z2QN() {}

inline Z2Z2QN &Z2Z2QN::operator=(const Z2Z2QN &rhs) {
  for (size_t i = 0; i < 2; i++) {
    z2vals_[i] = rhs.z2vals_[i];
  }
  hash_ = rhs.hash_;
  return *this;
}

inline Z2Z2QN Z2Z2QN::operator-() const {
  return *this;
}

inline Z2Z2QN &Z2Z2QN::operator+=(const Z2Z2QN &rhs) {
  z2vals_[0] = (z2vals_[0] + rhs.z2vals_[0]) % 2;
  z2vals_[1] = (z2vals_[1] + rhs.z2vals_[1]) % 2;
  hash_ = CalcHash_();
  return *this;
}

inline Z2Z2QN Z2Z2QN::operator+(const Z2Z2QN &rhs) const {
  Z2Z2QN res(*this);
  res += rhs;
  return res;
}

inline Z2Z2QN Z2Z2QN::operator-(const Z2Z2QN &rhs) const {
  Z2Z2QN res(*this);
  res += (-rhs);
  return res;
}

inline
void Z2Z2QN::StreamRead(std::istream &is) {
  is >> z2vals_[0];
  is >> z2vals_[1];
  is >> hash_;
}
inline
void Z2Z2QN::StreamWrite(std::ostream &os) const {
  os << z2vals_[0] << "\n"
     << z2vals_[1] << "\n"
     << hash_ << "\n";
}
inline
void Z2Z2QN::Show(const size_t indent_level) const {
  std::cout << IndentPrinter(indent_level)
            << "Z2Z2QN:  ("
            << z2vals_[0]
            << ", "
            << z2vals_[1]
            << ")"
            << "\n";
}
inline
size_t Z2Z2QN::CalcHash_() const {
  size_t hash_val = z2vals_[1] * 2 + z2vals_[0];
  return ((hash_val << 10) | (hash_val >> 54)); // To avoid collide of QNSector
}

inline std::istream &operator>>(std::istream &is, Z2Z2QN &qn) {
  qn.StreamRead(is);
  return is;
}

inline std::ostream &operator<<(std::ostream &os, const Z2Z2QN &qn) {
  qn.StreamWrite(os);
  return os;
}

inline size_t Hash(const Z2Z2QN &qn) { return qn.Hash(); }

}//special_qn
}//qlten

#endif //QLTEN_QLTENSOR_SPECIAL_QN_Z2Z2QN_H
