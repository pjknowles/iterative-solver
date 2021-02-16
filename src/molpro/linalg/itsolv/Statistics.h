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
  int r_scal = 0;
  int rr_dot = 0;
  int rr_axpy = 0;
  int rr_copy = 0;
  int rr_gemm_inner = 0;
  int rr_gemm_outer = 0;
  int q_scal = 0;
  int qr_dot = 0;
  int qr_axpy = 0;
  int qr_copy = 0;
  int qr_gemm_inner = 0;
  int qr_gemm_outer = 0;
  int rq_dot = 0;
  int rq_axpy = 0;
  int rq_copy = 0;
  int rq_gemm_inner = 0;
  int rq_gemm_outer = 0;
  int qq_dot = 0;
  int qq_axpy = 0;
  int qq_copy = 0;
  int qq_gemm_inner = 0;
  int qq_gemm_outer = 0;
  int p_scal = 0;
  int rp_dot = 0;
  int rp_axpy = 0;
  int rp_gemm_inner = 0;
  int rp_gemm_outer = 0;
  int qp_dot = 0;
  int qp_axpy = 0;
  int qp_gemm_inner = 0;
  int qp_gemm_outer = 0;
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
  if (statistics.r_scal > 0)
    o << statistics.r_scal << " scaling of R vectors, ";
  if (statistics.rr_dot > 0)
    o << statistics.rr_dot << " dot products between R vectors, ";
  if (statistics.rr_axpy > 0)
    o << statistics.rr_axpy << " axpy operations between R vectors, ";
  if (statistics.rr_copy > 0)
    o << statistics.rr_copy << " copy operations between R vectors, ";
  if (statistics.rr_gemm_inner > 0)
    o << statistics.rr_gemm_inner << " gemm_inner operations between R vectors, ";
  if (statistics.rr_gemm_outer > 0)
    o << statistics.rr_gemm_outer << " gemm_outer operations between R vectors, ";
  if (statistics.q_scal > 0)
    o << statistics.q_scal << " scaling of Q vectors, ";
  if (statistics.qq_dot > 0)
    o << statistics.qq_dot << " dot products between Q vectors, ";
  if (statistics.qq_axpy > 0)
    o << statistics.qq_axpy << " axpy operations between Q vectors, ";
  if (statistics.qq_copy > 0)
    o << statistics.qq_copy << " copy operations between Q vectors, ";
  if (statistics.qq_gemm_inner > 0)
    o << statistics.qq_gemm_inner << " gemm_inner operations between Q vectors, ";
  if (statistics.qq_gemm_outer > 0)
    o << statistics.rr_gemm_outer << " gemm_outer operations between Q vectors, ";
  if (statistics.qr_dot > 0)
    o << statistics.qr_dot << " dot products between Q and R vectors, ";
  if (statistics.qr_axpy > 0)
    o << statistics.qr_axpy << " axpy operations between Q and R vectors, ";
  if (statistics.qr_copy > 0)
    o << statistics.qr_copy << " copy operations between Q and R vectors, ";
  if (statistics.qr_gemm_inner > 0)
    o << statistics.qr_gemm_inner << " gemm_inner operations between Q and R vectors, ";
  if (statistics.qr_gemm_outer > 0)
    o << statistics.qr_gemm_outer << " gemm_outer operations between Q and R vectors, ";
  if (statistics.rq_dot > 0)
    o << statistics.rq_dot << " dot products between R and Q vectors, ";
  if (statistics.rq_axpy > 0)
    o << statistics.rq_axpy << " axpy operations between R and Q vectors, ";
  if (statistics.rq_copy > 0)
    o << statistics.rq_copy << " copy operations between R and Q vectors, ";
  if (statistics.rq_gemm_inner > 0)
    o << statistics.rq_gemm_inner << " gemm_inner operations between R and Q vectors, ";
  if (statistics.rq_gemm_outer > 0)
    o << statistics.rq_gemm_outer << " gemm_outer operations between R and Q vectors, ";
  if (statistics.p_scal > 0)
    o << statistics.p_scal << " scaling of Q vectors, ";
  if (statistics.rp_dot > 0)
    o << statistics.rp_dot << " dot products between R and P vectors, ";
  if (statistics.rp_axpy > 0)
    o << statistics.rp_axpy << " axpy operations between R and P vectors, ";
  if (statistics.rp_gemm_inner > 0)
    o << statistics.rp_gemm_inner << " gemm_inner operations between R and P vectors, ";
  if (statistics.rp_gemm_outer > 0)
    o << statistics.rp_gemm_outer << " gemm_outer operations between R and P vectors, ";
  if (statistics.qp_dot > 0)
    o << statistics.qp_dot << " dot products between Q and P vectors, ";
  if (statistics.qp_axpy > 0)
    o << statistics.qp_axpy << " axpy operations between Q and P vectors, ";
  if (statistics.qp_gemm_inner > 0)
    o << statistics.qp_gemm_inner << " gemm_inner operations between Q and P vectors, ";
  if (statistics.qp_gemm_outer > 0)
    o << statistics.qp_gemm_outer << " gemm_outer operations between Q and P vectors, ";
  o << "\b\b";

  return o;
}
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_STATISTICS_H_
