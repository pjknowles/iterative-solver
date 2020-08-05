#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H

#include <molpro/linalg/array/ArrayHandler.h>

namespace molpro {
namespace linalg {
namespace array {

template <class AL, class AR = AL>
class ArrayHandlerDistr : public ArrayHandler<AL, AR> {
public:
  using typename ArrayHandler<AL, AR>::value_type_L;
  using typename ArrayHandler<AL, AR>::value_type_R;
  using typename ArrayHandler<AL, AR>::value_type;
  using ArrayHandler<AL, AR>::copySame;
  using ArrayHandler<AL, AR>::copy;
  using ArrayHandler<AL, AR>::scal;
  using ArrayHandler<AL, AR>::fill;
  using ArrayHandler<AL, AR>::axpy;
  using ArrayHandler<AL, AR>::dot;
  using typename ArrayHandler<AL, AR>::ProxyHandle;

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

protected:
  using ArrayHandler<AL, AR>::lazy_handle;
  using ArrayHandler<AL, AR>::error;

  AR copyRR(const AR &source) override { return {source}; };
  AL copyLL(const AL &source) override { return {source}; };
  AR copyLR(const AL &source) override { return {source}; }
  AL copyRL(const AR &source) override { return {source}; }

  void scalL(value_type alpha, AL &x) override { x.scal(alpha); }
  void scalR(value_type alpha, AR &x) override { x.scal(alpha); }
  void fillL(value_type alpha, AL &x) override { x.fill(alpha); }
  void fillR(value_type alpha, AR &x) override { x.fill(alpha); }

  void axpyLL(value_type alpha, const AL &x, AL &y) override { y.axpy(alpha, x); }
  void axpyLR(value_type alpha, const AL &x, AR &y) override { y.axpy(alpha, x); }
  void axpyRL(value_type alpha, const AR &x, AL &y) override { y.axpy(alpha, x); }
  void axpyRR(value_type alpha, const AR &x, AR &y) override { y.axpy(alpha, x); }

  value_type dotLL(const AL &x, const AL &y) { return x.dot(y); }
  value_type dotLR(const AL &x, const AR &y) { return x.dot(y); }
  value_type dotRL(const AR &x, const AL &y) { return x.dot(y); }
  value_type dotRR(const AR &x, const AR &y) { return x.dot(y); }
};

} // namespace array
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERDISTR_H
