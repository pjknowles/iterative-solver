#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
namespace molpro {
namespace linalg {
namespace iterativesolver {
/*!
 * @brief Information about performance of IterativeSolver instance
 */
struct Statistics {
  int iterations = 0;
  int q_creations = 0;
  int p_creations = 0;
  int q_deletions = 0;
};
} // namespace iterativesolver
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
