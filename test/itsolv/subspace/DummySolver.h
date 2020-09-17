#ifndef LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYSOLVER_H
#define LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYSOLVER_H
#include <molpro/linalg/itsolv/IterativeSolver.h>

using molpro::linalg::itsolv::IterativeSolver;
using molpro::linalg::itsolv::Statistics;

namespace {
template <class R, class Q, class P>
struct DummySolver : IterativeSolver<R, Q, P> {
  using typename IterativeSolver<R, Q, P>::value_type;
  using typename IterativeSolver<R, Q, P>::scalar_type;
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action) override { return 0; };
  size_t add_vector(std::vector<R>& parameters, std::vector<R>& action, std::vector<P>& parametersP) override {
    return 0;
  };
  size_t add_p(std::vector<P>& Pvectors, const value_type* PP, std::vector<R>& parameters, std::vector<R>& action,
               std::vector<P>& parametersP) override {
    return 0;
  };
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual) override{};
  void solution(const std::vector<int>& roots, std::vector<R>& parameters, std::vector<R>& residual,
                std::vector<P>& parametersP) override{};
  std::vector<size_t> suggest_p(const std::vector<R>& solution, const std::vector<R>& residual, size_t maximumNumber,
                                double threshold) override {
    return {};
  };

  const std::vector<int>& working_set() const override { return ws; };
  const std::vector<scalar_type>& errors() const override { return er; };
  const Statistics& statistics() const override { return st; };
  std::vector<int> ws;
  std::vector<scalar_type> er;
  Statistics st;
};
} // namespace

#endif // LINEARALGEBRA_TEST_ITSOLV_SUBSPACE_DUMMYSOLVER_H
