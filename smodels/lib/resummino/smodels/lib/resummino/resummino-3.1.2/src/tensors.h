#ifndef TENSORS_H_
#define TENSORS_H_

#include "Fastor/Fastor.h"
#include "clooptools.h"
using namespace Fastor;

extern Tensor<ComplexType, 4, 4, 3> eps_x;
const int set_eps();
// const int ret = set_eps();

enum { I, J, K, L, M, O };

Tensor<ComplexType, 4> eps_mult(const Tensor<ComplexType, 3> &preFact,
                                const Tensor<ComplexType, 4> &integral);

template <size_t... Rest, FASTOR_INDEX... All>
Tensor<ComplexType, 4, Rest...> eps_mult(const Tensor<ComplexType, 3> &preFact,
                                         const Tensor<ComplexType, 4, Rest...> &integral);

/**
 *  Passarino Veltman access through Fastor tensors
 */

template <int i> ComplexType A(const double m1s);

// TODO no idea why this can't be moved to tensors.c like GetB/C
template <int i> Tensor<ComplexType, 4> GetA(const double m1s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = A<i>(m1s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = A<i>(m1s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = A<i>(m1s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = A<i>(m1s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}

template <int i> ComplexType B(const double p1s, const double m1s, const double m2s);
template <int i, int j> ComplexType B(const double p1s, const double m1s, const double m2s);

template <int i> Tensor<ComplexType, 4> GetB(const double p1s, const double m1s, const double m2s);

template <int i, int j>
Tensor<ComplexType, 4> GetB(const double p1s, const double m1s, const double m2s);

Tensor<ComplexType, 2, 4> GetBi(const double p1s, const double m1s, const double m2s);

template <int i> ComplexType DB(const double p1s, const double m1s, const double m2s);
template <int i, int j> ComplexType DB(const double p1s, const double m1s, const double m2s);

template <int i> Tensor<ComplexType, 4> GetDB(const double p1s, const double m1s, const double m2s);

template <int i, int j>
Tensor<ComplexType, 4> GetDB(const double p1s, const double m1s, const double m2s);
template <int i>
ComplexType C(const double p1s, const double p2s, const double p3s, const double m1s,
              const double m2s, const double m3s);
template <int i, int j>
ComplexType C(const double p1s, const double p2s, const double p3s, const double m1s,
              const double m2s, const double m3s);
template <int i>
ComplexType _C(const double p1s, const double p2s, const double p3s, const double m1s,
               const double m2s, const double m3s) {
  return C<i>(p1s, p2s, p3s, m1s, m2s, m3s);
};
template <int i, int j>
ComplexType _C(const double p1s, const double p2s, const double p3s, const double m1s,
               const double m2s, const double m3s) {
  return C<i, j>(p1s, p2s, p3s, m1s, m2s, m3s);
};

template <int i>
Tensor<ComplexType, 4> GetC(const double p1s, const double p2s, const double p3s, const double m1s,
                            const double m2s, const double m3s);
template <int i, int j>
Tensor<ComplexType, 4> GetC(const double p1s, const double p2s, const double p3s, const double m1s,
                            const double m2s, const double m3s);

Tensor<ComplexType, 3, 4> GetCi(const double p1s, const double p2s, const double p3s,
                                const double m1s, const double m2s, const double m3s);
Tensor<ComplexType, 4, 4> GetCij(const double p1s, const double p2s, const double p3s,
                                 const double m1s, const double m2s, const double m3s);
Tensor<ComplexType, 7, 4> GetCall(const double p1s, const double p2s, const double p3s,
                                  const double m1s, const double m2s, const double m3s);

template <int i>
ComplexType D(const double p1s, const double p2s, const double p3s, const double p4s,
              const double p1p2s, const double p2p3s, const double m1s, const double m2s,
              const double m3s, const double m4s);
template <int i, int j>
ComplexType D(const double p1s, const double p2s, const double p3s, const double p4s,
              const double p1p2s, const double p2p3s, const double m1s, const double m2s,
              const double m3s, const double m4s);

// template<int i>
// Tensor<ComplexType,4> GetD(const double p1s, const double p2s, const double
// p3s, const double p4s,const double p1p2s,const double p2p3s, const double
// m1s, const double m2s, const double m3s, const double m4s); template<int
// i,int j> Tensor<ComplexType,4> GetD(const double p1s, const double p2s, const
// double p3s, const double p4s,const double p1p2s,const double p2p3s, const
// double m1s, const double m2s, const double m3s, const double m4s);

template <int i, int j>
Tensor<ComplexType, 4> GetB(const double p1s, const double m1s, const double m2s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = B<i, j>(p1s, m1s, m2s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = B<i, j>(p1s, m1s, m2s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = B<i, j>(p1s, m1s, m2s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = B<i, j>(p1s, m1s, m2s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}

template <int i>
Tensor<ComplexType, 4> GetDB(const double p1s, const double m1s, const double m2s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = DB<i>(p1s, m1s, m2s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = DB<i>(p1s, m1s, m2s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = DB<i>(p1s, m1s, m2s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = DB<i>(p1s, m1s, m2s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}
template <int i, int j>
Tensor<ComplexType, 4> GetDB(const double p1s, const double m1s, const double m2s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = DB<i, j>(p1s, m1s, m2s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = DB<i, j>(p1s, m1s, m2s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = DB<i, j>(p1s, m1s, m2s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = DB<i, j>(p1s, m1s, m2s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}
template <int i>
Tensor<ComplexType, 4> GetD(const double p1s, const double p2s, const double p3s, const double p4s,
                            const double p1p2s, const double p2p3s, const double m1s,
                            const double m2s, const double m3s, const double m4s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = D<i>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = D<i>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = D<i>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = D<i>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}
template <int i, int j>
Tensor<ComplexType, 4> GetD(const double p1s, const double p2s, const double p3s, const double p4s,
                            const double p1p2s, const double p2p3s, const double m1s,
                            const double m2s, const double m3s, const double m4s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = D<i, j>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = D<i, j>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = D<i, j>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = D<i, j>(p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}

template <int i>
Tensor<ComplexType, 4> UVGetB(const double p1s, const double m1s, const double m2s);

template <int i> Tensor<ComplexType, 4> UVGetA(const double m2s);

template <size_t N> Tensor<ComplexType, N> real(Tensor<ComplexType, N> in) {
  return (in + conj(in)) / 2;
}
template <size_t N> Tensor<ComplexType, N> imag(Tensor<ComplexType, N> in) {
  return (in - conj(in)) / 2;
}

#endif