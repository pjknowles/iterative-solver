#ifndef LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#define LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>

using _Rvector_ = std::vector<double>;
class ExampleProblem : public molpro::linalg::itsolv::Problem<_Rvector_> {
protected:
  double matrix(int i, int j) const { return i == j ? i + 1 : 0.001 * ((i + j)%n); }

public:
  const size_t n;
  ExampleProblem(int n = 10) : n(n) {}

  void precondition(const VecRef<_Rvector_> &action,
                    const std::vector<typename _Rvector_::value_type> &shift) const override {
    for (int k = 0; k < action.size(); k++) {
      auto &a = action[k].get();
      for (int i = 0; i < a.size(); i++)
        a[i] /= (matrix(i, i) - shift[k] + 1e-15);
    }
  }

  double residual(const _Rvector_ &v, _Rvector_ &a) const override {
    double value = 0;
    for (int i = 0; i < a.size(); i++) {
      a[i] = 0;
      for (int j = 0; j < a.size(); j++)
        a[i] += matrix(i, j) * (v[j] - 1);
      value += 0.5 * a[i] * (v[i] - 1);
    }
    return value;
  }

  void action(const CVecRef<_Rvector_> &parameters, const VecRef<_Rvector_> &actions) const override {
    for (size_t k = 0; k < parameters.size(); k++) {
      const auto &v = parameters[k].get();
      auto &a = actions[k].get();
      for (int i = 0; i < a.size(); i++) {
        a[i] = 0;
        for (int j = 0; j < a.size(); j++)
          a[i] += matrix(i, j) * v[j];
      }
    }
  }
  friend
  std::ostream& operator<<(std::ostream& s, const ExampleProblem& problem) ;
};
std::ostream& operator<<(std::ostream& s, const ExampleProblem& problem) {
  s<<"ExampleProblem: matrix";
  for (int i=0;i<problem.n;i++) {
    s<<"\n"<<i;
    for(int j=0;j<problem.n;j++)s<<" "<<problem.matrix(i,j);
  }
  return s;
}
#endif // LINEARALGEBRA_EXAMPLES_EXAMPLEPROBLEM_H_
