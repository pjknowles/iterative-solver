#include "QSpace.h"

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace qspace {

void merge_subspace_qq(size_t i, size_t j, double a, double b, SubspaceData& qq) {
  assert(i < j);
  for (auto d : {EqnData::S, EqnData::H}) {
    auto rows = qq[d].rows();
    auto cols = qq[d].cols();
    auto diag = a * a * qq[d](i, i) + a * b * qq[d](i, j) + a * b * qq[d](j, i) + b * b * qq[d](j, j);
    auto row_i = qq[d].slice({i, 0}, {i + 1, cols});
    auto row_j = qq[d].slice({j, 0}, {j + 1, cols});
    auto col_i = qq[d].slice({0, i}, {rows, i + 1});
    auto col_j = qq[d].slice({0, j}, {rows, j + 1});
    row_i.scal(a);
    row_i.axpy(b, row_j);
    col_i.scal(a);
    col_i.axpy(b, col_j);
    qq[d](i, i) = diag;
    qq[d].remove_row_col(j, j);
  }
}

void merge_subspace_qr(size_t i, size_t j, double a, double b, SubspaceData& qr) {
  assert(i < j);
  for (auto d : {EqnData::S, EqnData::H}) {
    auto cols = qr[d].cols();
    auto row_i = qr[d].slice({i, 0}, {i + 1, cols});
    auto row_j = qr[d].slice({j, 0}, {j + 1, cols});
    row_i.scal(a);
    row_i.axpy(b, row_j);
    qr[d].remove_row(j);
  }
}

void merge_subspace_rq(size_t i, size_t j, double a, double b, SubspaceData& rq) {
  assert(i < j);
  for (auto d : {EqnData::S, EqnData::H}) {
    auto rows = rq[d].rows();
    auto col_i = rq[d].slice({0, i}, {rows, i + 1});
    auto col_j = rq[d].slice({0, j}, {rows, j + 1});
    col_i.scal(a);
    col_i.axpy(b, col_j);
    rq[d].remove_col(j);
  }
}

void erase_subspace(size_t i, SubspaceData& qq, SubspaceData& qr, SubspaceData& rq) {
  for (auto d : {EqnData::S, EqnData::H}) {
    qq[d].remove_row_col(i, i);
    qr[d].remove_row(i);
    rq[d].remove_col(i);
  }
}

} // namespace qspace
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro