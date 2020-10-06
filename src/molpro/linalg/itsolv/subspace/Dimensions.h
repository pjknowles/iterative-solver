#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
namespace xspace {
//! Stores partitioning of XSpace into P, Q and R blocks with sizes and offsets for each one
struct Dimensions {
  Dimensions() = default;
  Dimensions(size_t np, size_t nq, size_t nr, size_t nc) : nP(np), nQ(nq), nR(nr), nC(nc) {}
  size_t nP = 0;
  size_t nQ = 0;
  size_t nR = 0;
  size_t nC = 0;
  size_t nX = nP + nQ + nR + nC;
  size_t oP = 0;
  size_t oQ = nP;
  size_t oR = oQ + nQ;
  size_t oC = oR + nR;
};
}
} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_DIMENSIONS_H
