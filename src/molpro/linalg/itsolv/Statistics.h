#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
#include <molpro/linalg/itsolv/ArrayHandlers.h>
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
  std::string rq_ops = "";
  std::string qr_ops = "";
  std::string rr_ops = "";
  std::string qq_ops = "";
  std::string rp_ops = "";
  std::string qp_ops = "";
};

template <typename R, typename Q, typename P>
void read_handler_counts(std::shared_ptr<Statistics> stats, std::shared_ptr<ArrayHandlers<R, Q, P>> handlers) {
  stats->rr_ops = handlers->rr().counter_to_string("R", "R");
  stats->qr_ops = handlers->qr().counter_to_string("Q", "R");
  stats->rq_ops = handlers->rq().counter_to_string("R", "Q");
  stats->qq_ops = handlers->qq().counter_to_string("Q", "Q");
  stats->rp_ops = handlers->rp().counter_to_string("R", "P");
  stats->qp_ops = handlers->qp().counter_to_string("Q", "P");
}

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
  o << statistics.rr_ops << " ";
  o << statistics.qr_ops << " ";
  o << statistics.rq_ops << " ";
  o << statistics.qq_ops << " ";
  o << statistics.rp_ops << " ";
  o << statistics.qp_ops << " ";
  o << "\b\b";

  return o;
}
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
