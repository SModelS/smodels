#include "tensors.h"

Tensor<ComplexType, 4, 4, 3> eps_x = {
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 0}},
    {{0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 0, 0}},
    {{0, 0, 0}, {0, 0, 0}, {0, 1, 0}, {1, 0, 0}}};
const int set_eps() {
  eps_x.zeros();

  eps_x(0, 0, 0) = 1;
  eps_x(1, 1, 0) = 1;
  eps_x(2, 2, 0) = 1;

  eps_x(0, 1, 1) = 1;
  eps_x(1, 2, 1) = 1;

  eps_x(0, 2, 2) = 1;

  eps_x(3, 2, 1) = 1;
  eps_x(3, 3, 0) = 1;
  return 0;
}

Tensor<ComplexType, 4> eps_mult(const Tensor<ComplexType, 3> &preFact,
                                const Tensor<ComplexType, 4> &integral) {
  // set_eps();
  return einsum<Index<L>, Index<M, I, L>, Index<M>>(preFact, eps_x, integral);
}

template <size_t... Rest, FASTOR_INDEX... All>
Tensor<ComplexType, 4, Rest...>
eps_mult(const Tensor<ComplexType, 3> &preFact,
         const Tensor<ComplexType, 4, Rest...> &integral) {
  return einsum<Index<L>, Index<M, I, L>, Index<M, All...>>(preFact, eps_x,
                                                            integral);
}

template <> ComplexType A<0>(const double m1s) { return A0(m1s); }

template <>
ComplexType B<0>(const double p1s, const double m1s, const double m2s) {
  return B0(p1s, m1s, m2s);
}
template <>
ComplexType B<1>(const double p1s, const double m1s, const double m2s) {
  return B1(p1s, m1s, m2s);
}

template <>
ComplexType B<0, 0>(const double p1s, const double m1s, const double m2s) {
  return B0i(bb00, p1s, m1s, m2s);
}
template <>
ComplexType B<1, 1>(const double p1s, const double m1s, const double m2s) {
  return B0i(bb11, p1s, m1s, m2s);
}
template <>
ComplexType DB<0>(const double p1s, const double m1s, const double m2s) {
  return DB0(p1s, m1s, m2s);
}
template <>
ComplexType DB<1>(const double p1s, const double m1s, const double m2s) {
  return DB1(p1s, m1s, m2s);
}
template <>
ComplexType DB<0, 0>(const double p1s, const double m1s, const double m2s) {
  return DB00( p1s, m1s, m2s);
}
template <>
ComplexType DB<1, 1>(const double p1s, const double m1s, const double m2s) {
  return DB11( p1s, m1s, m2s);
}


template <int i>
Tensor<ComplexType, 4> GetB(const double p1s, const double m1s,
                            const double m2s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = B<i>(p1s, m1s, m2s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = B<i>(p1s, m1s, m2s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = B<i>(p1s, m1s, m2s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = B<i>(p1s, m1s, m2s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}

Tensor<ComplexType, 2, 4> GetBi(const double p1s, const double m1s,
                                const double m2s) {
  Tensor<ComplexType, 2, 4> Bi;
  Bi(0, all) = GetB<0>(p1s, m1s, m2s);
  Bi(1, all) = GetB<1>(p1s, m1s, m2s);
  return Bi;
}

template <>
ComplexType C<0>(const double p1s, const double p2s, const double p3s,
                 const double m1s, const double m2s, const double m3s) {
  return C0i(cc0, p1s, p2s, p3s, m1s, m2s, m3s);
}
template <>
ComplexType C<1>(const double p1s, const double p2s, const double p3s,
                 const double m1s, const double m2s, const double m3s) {
  return C0i(cc1, p1s, p2s, p3s, m1s, m2s, m3s);
}
template <>
ComplexType C<2>(const double p1s, const double p2s, const double p3s,
                 const double m1s, const double m2s, const double m3s) {
  return C0i(cc2, p1s, p2s, p3s, m1s, m2s, m3s);
}

template <int i>
Tensor<ComplexType, 4> GetC(const double p1s, const double p2s,
                            const double p3s, const double m1s,
                            const double m2s, const double m3s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = C<i>(p1s, p2s, p3s, m1s, m2s, m3s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = C<i>(p1s, p2s, p3s, m1s, m2s, m3s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = C<i>(p1s, p2s, p3s, m1s, m2s, m3s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = C<i>(p1s, p2s, p3s, m1s, m2s, m3s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}
template <int i, int j>
Tensor<ComplexType, 4> GetC(const double p1s, const double p2s,
                            const double p3s, const double m1s,
                            const double m2s, const double m3s) {
  ComplexType b0_IR_0, b0_IR_1, b0_IR_2, b0_UV_1;
  setlambda(0);
  b0_IR_0 = C<i, j>(p1s, p2s, p3s, m1s, m2s, m3s);
  setuvdiv(0.0);
  setlambda(-1);
  b0_IR_1 = C<i, j>(p1s, p2s, p3s, m1s, m2s, m3s);
  setuvdiv(1.0);
  setlambda(-1);
  b0_UV_1 = C<i, j>(p1s, p2s, p3s, m1s, m2s, m3s) - b0_IR_1;
  setlambda(-2);
  b0_IR_2 = C<i, j>(p1s, p2s, p3s, m1s, m2s, m3s);
  return Tensor<ComplexType, 4>({b0_IR_2, b0_IR_1, b0_IR_0, b0_UV_1});
}

Tensor<ComplexType, 3, 4> GetCi(const double p1s, const double p2s,
                                const double p3s, const double m1s,
                                const double m2s, const double m3s) {
  Tensor<ComplexType, 3, 4> Ci;
  Ci(0, all) = GetC<0>(p1s, p2s, p3s, m1s, m2s, m3s);
  Ci(1, all) = GetC<1>(p1s, p2s, p3s, m1s, m2s, m3s);
  Ci(2, all) = GetC<2>(p1s, p2s, p3s, m1s, m2s, m3s);
  return Ci;
}

template <>
ComplexType C<0, 0>(const double p1s, const double p2s, const double p3s,
                    const double m1s, const double m2s, const double m3s) {
  return C0i(cc00, p1s, p2s, p3s, m1s, m2s, m3s);
}
template <>
ComplexType C<1, 1>(const double p1s, const double p2s, const double p3s,
                    const double m1s, const double m2s, const double m3s) {
  return C0i(cc11, p1s, p2s, p3s, m1s, m2s, m3s);
}
template <>
ComplexType C<1, 2>(const double p1s, const double p2s, const double p3s,
                    const double m1s, const double m2s, const double m3s) {
  return C0i(cc12, p1s, p2s, p3s, m1s, m2s, m3s);
}
template <>
ComplexType C<2, 2>(const double p1s, const double p2s, const double p3s,
                    const double m1s, const double m2s, const double m3s) {
  return C0i(cc22, p1s, p2s, p3s, m1s, m2s, m3s);
}

Tensor<ComplexType, 4, 4> GetCij(const double p1s, const double p2s,
                                 const double p3s, const double m1s,
                                 const double m2s, const double m3s) {
  Tensor<ComplexType, 4, 4> Cij;
  Cij(0, all) = GetC<0, 0>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(1, all) = GetC<1, 1>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(2, all) = GetC<1, 2>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(3, all) = GetC<2, 2>(p1s, p2s, p3s, m1s, m2s, m3s);
  return Cij;
}

Tensor<ComplexType, 7, 4> GetCall(const double p1s, const double p2s,
                                  const double p3s, const double m1s,
                                  const double m2s, const double m3s) {
  Tensor<ComplexType, 7, 4> Cij;
  Cij(0, all) = GetC<0>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(1, all) = GetC<1>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(2, all) = GetC<2>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(3, all) = GetC<0, 0>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(4, all) = GetC<1, 1>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(5, all) = GetC<1, 2>(p1s, p2s, p3s, m1s, m2s, m3s);
  Cij(6, all) = GetC<2, 2>(p1s, p2s, p3s, m1s, m2s, m3s);
  return Cij;
}

template <>
ComplexType D<0>(const double p1s, const double p2s, const double p3s,
                 const double p4s, const double p1p2s, const double p2p3s,
                 const double m1s, const double m2s, const double m3s,
                 const double m4s) {
  return D0i(dd0, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<1>(const double p1s, const double p2s, const double p3s,
                 const double p4s, const double p1p2s, const double p2p3s,
                 const double m1s, const double m2s, const double m3s,
                 const double m4s) {
  return D0i(dd1, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<2>(const double p1s, const double p2s, const double p3s,
                 const double p4s, const double p1p2s, const double p2p3s,
                 const double m1s, const double m2s, const double m3s,
                 const double m4s) {
  return D0i(dd2, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<3>(const double p1s, const double p2s, const double p3s,
                 const double p4s, const double p1p2s, const double p2p3s,
                 const double m1s, const double m2s, const double m3s,
                 const double m4s) {
  return D0i(dd3, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}

template <>
ComplexType D<0, 0>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd00, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<1, 1>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd11, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<1, 2>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd12, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<1, 3>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd13, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<2, 2>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd22, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<2, 3>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd23, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}
template <>
ComplexType D<3, 3>(const double p1s, const double p2s, const double p3s,
                    const double p4s, const double p1p2s, const double p2p3s,
                    const double m1s, const double m2s, const double m3s,
                    const double m4s) {
  return D0i(dd33, p1s, p2s, p3s, p4s, p1p2s, p2p3s, m1s, m2s, m3s, m4s);
}

template <>
Tensor<ComplexType, 4> UVGetB<0>(const double p1s, const double m1s,
                                 const double m2s) {
  Tensor<ComplexType, 4> r = {0., 0., 0., 1.0};
  return r;
}
template <>
Tensor<ComplexType, 4> UVGetB<1>(const double p1s, const double m1s,
                                 const double m2s) {
  Tensor<ComplexType, 4> r = {0., 0., 0., -0.5};
  return r;
}
template <> Tensor<ComplexType, 4> UVGetA<0>(const double m1s) {
  Tensor<ComplexType, 4> r = {0., 0., 0., m1s};
  return r;
}