#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
namespace molpro::linalg::itsolv::subspace {
//! Stores partitioning of XSpace into P, Q and R blocks with sizes and offsets for each one
struct Dimensions {
  Dimensions() = default;
  Dimensions(size_t np, size_t nq, size_t nc) : nP(np), nQ(nq), nD(nc) {}
  size_t nP = 0;
  size_t nQ = 0;
  size_t nD = 0;
  size_t nX = nP + nQ + nD;
  size_t oP = 0;
  size_t oQ = nP;
  size_t oD = oQ + nQ;
  size_t nRHS = 0; //!< number of rigt-hand-side vectors in the system of linear equations
};
} // namespace molpro::linalg::itsolv::subspace

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
