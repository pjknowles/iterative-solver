#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
namespace molpro::linalg::itsolv::subspace::xspace {
//! Stores partitioning of XSpace into P, Q and R blocks with sizes and offsets for each one
struct Dimensions {
  Dimensions() = default;
  Dimensions(size_t np, size_t nq, size_t nc) : nP(np), nQ(nq), nC(nc) {}
  size_t nP = 0;
  size_t nQ = 0;
  size_t nC = 0;
  size_t nX = nP + nQ + nC;
  size_t oP = 0;
  size_t oQ = nP;
  size_t oC = oQ + nQ;
};
} // namespace molpro::linalg::itsolv::subspace::xspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
