#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
#include <ostream>
namespace molpro::linalg::itsolv {
/*!
 * @brief Information about performance of IterativeSolver instance
 */
struct Statistics {
  int iterations = 0;
  int r_creations = 0;
  int q_creations = 0;
  int p_creations = 0;
  int q_deletions = 0;
  int d_creations = 0;
  int best_r_creations = 0;
  int current_r_creations = 0;
  int line_searches = 0;
  int line_search_steps = 0;
};
inline std::ostream& operator<<(std::ostream& o, const Statistics& statistics) {
  if (statistics.iterations > 0)
    o << statistics.iterations << " iterations, ";
  if (statistics.r_creations > 0)
    o << statistics.r_creations << " R vectors, ";
  if (statistics.q_creations != statistics.r_creations)
    o << statistics.q_creations << " Q creations, ";
  if (statistics.q_deletions > 0)
    o << statistics.q_deletions << " Q deletions, ";
  if (statistics.p_creations > 0)
    o << statistics.p_creations << " P vectors, ";
  if (statistics.d_creations > 0)
    o << statistics.d_creations << " D vectors, ";
  if (statistics.best_r_creations > 0)
    o << statistics.best_r_creations << " best R vectors, ";
  if (statistics.current_r_creations > 0)
    o << statistics.current_r_creations << " current R vectors, ";
  if (statistics.line_searches > 0)
    o << statistics.line_searches << " line searches, ";
  if (statistics.line_search_steps > 0)
    o << statistics.line_search_steps << " line search steps, ";
  o << "\b\b";

  return o;
}
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
