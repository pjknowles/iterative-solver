#include "molpro/iostream.h"
#include "molpro/linalg/IterativeSolver.h"
#include "molpro/linalg/OpaqueArray.h"
#include "molpro/linalg/PagedArray.h"
#include "molpro/linalg/SimpleArray.h"
#include "molpro/linalg/array/ArrayHandlerIterable.h"
#include "test.h"
#include <Eigen/Dense>
#include <ctime>
#include <memory>
#include <type_traits>

namespace molpro {
namespace linalg {

template <class T>
void syncr(T& x, std::true_type) {
  if (!x.synchronised())
    x.sync();
}

template <class T>
void syncr(T& x, std::false_type) {}

/*!
 * \brief Test iterative solution of linear eigenvalue problem
 * \param dimension The dimension of the test matrix
 * \param roots How many eigensolutions to find
 * \param verbosity How much to report
 * \param problem Selects which test matrix to use
 * \tparam ptype Concrete class template that implements LinearAlgebra::vectorSet
 * \tparam scalar Type of matrix elements
 */
template <class ptype>
static void DavidsonTest(size_t dimension, size_t roots = 1, int verbosity = 0, int problem = 0) {

  using scalar = typename LinearEigensystem<ptype>::scalar_type;
  using element = typename LinearEigensystem<ptype>::value_type;
  using vectorSet = std::vector<ptype>;
  static Eigen::Matrix<element, Eigen::Dynamic, Eigen::Dynamic> testmatrix;

  static struct {
    void operator()(const vectorSet& psx, vectorSet& outputs) const {
      for (size_t k = 0; k < psx.size(); k++) {
        Eigen::Matrix<element, Eigen::Dynamic, 1> x(testmatrix.rows());
        psx[k].get(&x[0], testmatrix.rows(), 0);
        Eigen::VectorXd res = testmatrix * x;
        outputs[k].put(&res[0], testmatrix.rows(), 0);
      }
    }
  } action;

  static struct {
    void operator()(vectorSet& psc, const vectorSet& psg, std::vector<scalar> shift, bool append = true) const {
      size_t n = testmatrix.rows();
      std::vector<element> psck(n);
      std::vector<element> psgk(n);
      for (size_t k = 0; k < psc.size(); k++) {
        psg[k].get(&psgk[0], n, 0);
        if (not append)
          psc[k].scal(0);
        psc[k].get(&psck[0], n, 0);
        for (size_t l = 0; l < n; l++)
          psck[l] -= psgk[l] / (testmatrix(l, l) + shift[k]);
        psc[k].put(&psck[0], n, 0);
      }
    }
  } update;

  molpro::cout << "Test IterativeSolver::LinearEigensystem dimension=" << dimension << ", roots=" << roots
               << ", problem=" << problem << std::endl;
  testmatrix.resize(dimension, dimension);
  for (size_t k = 0; k < dimension; k++)
    for (size_t l = 0; l < dimension; l++)
      if (problem == 0)
        testmatrix(l, k) = -1;
      else if (problem == 1)
        testmatrix(l, k) = l == k ? 1 + k * 5 : l + k + 2;
      else if (problem == 2)
        testmatrix(l, k) = (k == l ? k + 1 : 1);
      else if (problem == 3)
        testmatrix(l, k) = 1;
      else
        throw std::logic_error("invalid problem in DavidsonTest");
  if (problem == 3)
    testmatrix(0, 1) = testmatrix(1, 0) = 1;

  LinearEigensystem<ptype> d;
  molpro::linalg::array::ArrayHandlerIterable<ptype> handler{};
  d.m_roots = roots;
  d.m_verbosity = verbosity;
  d.m_maxIterations = dimension;
  vectorSet x;
  vectorSet g;
  for (size_t root = 0; root < (size_t)d.m_roots; root++) {
    x.emplace_back(dimension);
    g.emplace_back(dimension);
    x.back().scal(0);
    element one = 1;
    x.back().put(&one, 1, root);
  }

  for (size_t iteration = 0; iteration < dimension + 1; iteration++) {
    action(x, g);
    auto nwork = d.addVector(x, g);
    d.report();
    //    for (auto root=0; root < nwork; root++) {
    //      std::vector<double> xxx(dimension);
    //      x[root].get(xxx.data(),dimension,0);
    //      std::cout << "New working vector";
    //      for (const auto& e : xxx) std::cout << " "<<e;std::cout<<std::endl;
    //    }

    for (size_t root = 0; root < (size_t)d.m_roots; root++) {
      syncr(x[root], std::is_same<ptype, linalg::PagedArray<double>>{});
      syncr(g[root], std::is_same<ptype, linalg::PagedArray<double>>{});
    }
    std::vector<scalar> shift;
    for (size_t root = 0; root < (size_t)d.m_roots; root++)
      shift.push_back(-d.eigenvalues()[root] + 1e-14);
    update(x, g, shift);
    //    auto newp = d.suggestP(x, g, 3);
    // if (d.endIteration(x, g)) break;
    //    bool upd = d.endIteration(x, g);
    for (size_t root = 0; root < (size_t)d.m_roots; root++) {
      syncr(x[root], std::is_same<ptype, linalg::PagedArray<double>>{});
      syncr(g[root], std::is_same<ptype, linalg::PagedArray<double>>{});
    }
    if (nwork == 0)
      break;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>> es(testmatrix);
  //    molpro::cout << "true eigenvalues: "<<es.eigenvalues().head(d.m_roots).transpose()<<std::endl;
  //    molpro::cout << "true eigenvectors:\n"<<es.eigenvectors().leftCols(d.m_roots).transpose()<<std::endl;

  auto ev = d.eigenvalues();
  molpro::cout << "Eigenvalues: ";
  size_t root = 0;
  for (const auto& e : ev)
    molpro::cout << " " << e << "(error=" << e - es.eigenvalues()(root++) << ")";
  molpro::cout << std::endl;
  molpro::cout << "Reported errors: ";
  for (const auto& e : d.errors())
    molpro::cout << " " << e;
  molpro::cout << std::endl;

  std::vector<int> rootlist;
  for (size_t root = 0; root < (size_t)d.m_roots; root++)
    rootlist.push_back(root);
  d.solution(rootlist, x, g);
  for (size_t root = 0; root < (size_t)d.m_roots; root++) {
    syncr(x[root], std::is_same<ptype, linalg::PagedArray<double>>{});
  }
  action(x, g);
  for (size_t root = 0; root < (size_t)d.m_roots; root++) {
    syncr(x[root], std::is_same<ptype, linalg::PagedArray<double>>{});
    syncr(g[root], std::is_same<ptype, linalg::PagedArray<double>>{});
  }
  std::vector<scalar> errors;
  for (size_t root = 0; root < (size_t)d.m_roots; root++) {
    //g[root].axpy(-ev[root], x[root]);
    handler.axpy(-ev[root], x[root], g[root]);
    //errors.push_back(g[root].dot(g[root]));
    errors.push_back(handler.dot(g[root], g[root]));
  }
  //   molpro::cout << "Square residual norms: "; for (typename std::vector<T>::const_iterator e=errors.begin();
  //   e!=errors.end(); e++) molpro::cout<<" "<<*e;molpro::cout<<std::endl;
  molpro::cout << "Square residual norms: ";
  for (const auto& e : errors)
    molpro::cout << " " << e;
  molpro::cout << std::endl;
  // be noisy about obvious problems
  if (*std::max_element(errors.begin(), errors.end()) > 1e-7)
    throw std::runtime_error("IterativeSolver::LinearEigensystem has failed tests");
}

/*!
 * \brief Test the correct operation of the non-linear equation solver. If an error is found, an exception is thrown.
 * \param verbosity How much to print.
 * - -1 Nothing at all is printed.
 * - 0 (default) Just a message that the test is taking place.
 * - 1, 2, 3,... more detail.
 * \param maxDim Maximum DIIS dimension allowed
 * \param svdThreshold Residual threshold for inclusion of a vector in the DIIS state.
 * \param mode Whether to perform DIIS, KAIN, or nothing.
 * \param difficulty Level of numerical challenge, ranging from 0 to 1.
 */
template <class ptype>
void DIISTest(int verbosity = 0, size_t maxDim = 6, double svdThreshold = 1e-10, double difficulty = 0.1) {
  using scalar = typename DIIS<ptype>::scalar_type;
  static struct {
    void operator()(const ptype& psx, ptype& outputs) const {
      size_t n = 2;
      std::vector<scalar> psxk(n);
      std::vector<scalar> output(n);

      psx.get(&(psxk[0]), n, 0);
      output[0] = (2 * psxk[0] - 2 + 400 * psxk[0] * (psxk[0] * psxk[0] - psxk[1]));
      output[1] = (200 * (psxk[1] - psxk[0] * psxk[0])); // Rosenbrock
      outputs.put(&(output[0]), n, 0);
    }
  } _Rosenbrock_residual;

  static struct {
    void operator()(ptype& psc, const ptype& psg, std::vector<scalar> shift, bool append = true) const {
      size_t n = 2;
      std::vector<scalar> psck(n);
      std::vector<scalar> psgk(n);
      psg.get(&psgk[0], n, 0);
      if (append) {
        psc.get(&psck[0], n, 0);
        psck[0] -= psgk[0] / 700;
        psck[1] -= psgk[1] / 200;
      } else {
        psck[0] = -psgk[0] / 700;
        psck[1] = -psgk[1] / 200;
      }
      psc.put(&psck[0], n, 0);
      //    molpro::cout << "Rosenbrock updater, new psc="<<psc<<std::endl;
    }
  } _Rosenbrock_updater;
  ptype x(2);
  ptype g(2);
  DIIS<ptype> d;
  d.m_maxQ = maxDim;
  d.m_svdThreshold = svdThreshold;

  if (verbosity >= 0)
    molpro::cout << "Test DIIS::iterate, difficulty=" << difficulty << std::endl;
  d.m_verbosity = verbosity - 1;
  //  d.m_options["weight"]=2;
  //  return;
  std::vector<scalar> xxx(2);
  xxx[0] = xxx[1] = 1 - difficulty; // initial guess
  x.put(&xxx[0], 2, 0);
  //  molpro::cout << "initial guess " << x << std::endl;
  bool converged = false;
  for (int iteration = 1; iteration < 1000 && not converged; iteration++) {
    //   molpro::cout <<"start of iteration "<<iteration<<std::endl;
    _Rosenbrock_residual(x, g);
    //    molpro::cout << "residual: " << g;
    d.addVector(x, g);
    std::vector<scalar> shift;
    shift.push_back(1e-10);
    _Rosenbrock_updater(x, g, shift);
    converged = d.endIteration(x, g);
    x.get(&xxx[0], 2, 0);
    //    if (verbosity > 2)
    //      molpro::cout << "new x after iterate " << x.front() << std::endl;
    if (verbosity >= 0)
      molpro::cout << "iteration "
                   << iteration
                   //                   << ", Residual norm = " << std::sqrt(d.fLastResidual())
                   << ", Distance from solution = "
                   << std::sqrt((xxx[0] - 1) * (xxx[0] - 1) + (xxx[1] - 1) * (xxx[1] - 1))
                   << ", error = " << d.errors().front() << ", converged? " << converged << std::endl;
    //   molpro::cout <<"end of iteration "<<iteration<<std::endl;
  }

  x.get(&xxx[0], 2, 0);
  molpro::cout << "Distance from solution = " << std::sqrt((xxx[0] - 1) * (xxx[0] - 1) + (xxx[1] - 1) * (xxx[1] - 1))
               << std::endl;
}

#include <cstdlib>
// struct anharmonic {
//  Eigen::MatrixXd m_F;
//  double m_gamma;
//  size_t m_n;
//  anharmonic(){}
//  void set(size_t n, double alpha, double gamma)
//  {
//    m_gamma=gamma;
//    m_n=n;

//    m_F.resize(n,n);
//    for (size_t j=0; j<n; j++) {
//        for (size_t i=0; i<n; i++)
//          m_F(i,j)=-0.5 + (((double)rand())/RAND_MAX);
//        m_F(j,j) += (j*alpha+0.5);
//      }
//  }
//  ptype guess()
//  {
//    std::vector<T> r(m_n);
//    ptype result(m_n);
//    double value=0.3;
//    for (size_t k=0; k<m_n; k++) {
//        r[k]=value;
//        value=-value;
//      }
//    result.put(&r[0],m_n,0);
//    return result;
//  }
//};

// static anharmonic instance;

// static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<T> & psx, vectorSet<T> & outputs, std::vector<T> shift=std::vector<T>(), bool
//  append=false) const override {
//    std::vector<T> psxk(instance.m_n);
//    std::vector<T> output(instance.m_n);
//    psx.front()->get(&(psxk[0]),instance.m_n,0);
//    if (append)
//      outputs.front()->get(&(output[0]),instance.m_n,0);
//    else
//      outputs.front()->scal(0);

//    for (size_t i=0; i<instance.m_n; i++) {
//        output[i] = instance.m_gamma*psxk[i];
//        for (size_t j=0; j<instance.m_n; j++)
//          output[i] += instance.m_F(j,i)*psxk[j];
//      }
//    outputs.front()->put(&output[0],instance.m_n,0);
//  }
//} _anharmonic_residual;
// static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<T> & psg, vectorSet<T> & psc, std::vector<T> shift=std::vector<T>(), bool
//  append=false) const override {
//    std::vector<T> psck(instance.m_n);
//    std::vector<T> psgk(instance.m_n);
//    psg.front()->get(&psgk[0],instance.m_n,0);
//    if (append) {
//        psc.front()->get(&psck[0],instance.m_n,0);
//        for (size_t i=0; i<instance.m_n; i++)
//          psck[i] -= psgk[i]/instance.m_F(i,i);
//      } else {
//        for (size_t i=0; i<instance.m_n; i++)
//          psck[i] =- psgk[i]/instance.m_F(i,i);
//      }
//    psc.front()->put(&psck[0],instance.m_n,0);
//  }
//} _anharmonic_preconditioner;
// void DIIS::randomTest(size_t sample, size_t n, double alpha, double gamma, DIISmode_type mode)
//{

//  int nfail=0;
//  unsigned int iterations=0, maxIterations=0;
//  for (size_t repeat=0; repeat < sample; repeat++) {
//      instance.set(n,alpha,gamma);
//      DIIS d(_anharmonic_residual,_anharmonic_preconditioner);
//      d.setMode(mode);
//      d.m_verbosity=-1;
//      d.m_maxIterations=100000;
//      ptype gg(n); vectorSet<T> g; g.push_back(std::shared_ptr<ptype>(&gg));
//      ptype xx=instance.guess(); vectorSet<T> x; x.push_back(std::shared_ptr<ptype>(&xx));
//      if (not d.solve(g,x)) nfail++;
//      iterations+=d.iterations();
//      if (maxIterations<d.iterations())
//        maxIterations=d.iterations();
//    }
//  molpro::cout << "sample="<<sample<<", n="<<n<<", alpha="<<alpha<<", gamma="<<gamma<<", average
//  iterations="<<iterations/sample<<", maximum iterations="<<maxIterations<<", nfail="<<nfail<<std::endl;
//}

#include <cstdlib>
template <class ptype, class scalar = double>
void RSPTTest(size_t n, double alpha) { // TODO conversion not finished
  using vectorSet = std::vector<ptype>;
  static struct rsptpot {
    Eigen::MatrixXd m_F;
    size_t m_n;
    rsptpot() {}
    size_t m_reference;
    void set(size_t n, double alpha) {
      //        molpro::cout<<"rsptpot set"<<n<<std::endl;
      m_n = n;
      m_reference = 0; // asserting that m_F(0,0) is the lowest

      m_F.resize(n, n);
      for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < j; i++)
          m_F(i, j) = m_F(j, i) = -0.5 + (((double)rand()) / RAND_MAX);
        m_F(j, j) = (j * alpha - 1);
      }
      //      molpro::cout << "m_F:"<<std::endl<<m_F<<std::endl;
    }
    ptype guess() {
      std::vector<scalar> r(m_n);
      ptype result(m_n);
      for (size_t k = 0; k < m_n; k++)
        r[k] = 0;
      r[m_reference] = 1;
      result.put(&r[0], m_n, 0);
      return result;
    }

  } instance;

  static struct {
    void operator()(const vectorSet& psx, vectorSet& outputs, std::vector<scalar> shift = std::vector<scalar>(),
                    bool append = false) const {
      //        molpro::cout << "rsptpot_residual"<<std::endl;
      //        molpro::cout << "input "<<psx<<std::endl;
      std::vector<scalar> psxk(instance.m_n);
      std::vector<scalar> output(instance.m_n);
      psx.front()->get(&(psxk[0]), instance.m_n, 0);
      if (append)
        outputs.front()->get(&(output[0]), instance.m_n, 0);
      else
        outputs.front()->scal(0);
      for (size_t i = 0; i < instance.m_n; i++) {
        output[i] = 0;
        for (size_t j = 0; j < instance.m_n; j++) {
          output[i] += instance.m_F(j, i) * psxk[j];
        }
      }
      outputs.front()->put(&(output[0]), instance.m_n, 0);
      //        molpro::cout << "output "<<outputs<<std::endl;
    }
  } _rsptpot_residual;
  static struct {
    void operator()(vectorSet& psc, const vectorSet& psg, std::vector<scalar> shift = std::vector<scalar>(),
                    bool append = false) const {
      //        molpro::cout << "preconditioner input="<<psg<<std::endl;
      //      if (shift.front()==0)
      //          molpro::cout << "H0 not resolvent"<<std::endl;
      std::vector<scalar> psck(instance.m_n);
      std::vector<scalar> psgk(instance.m_n);
      psg.front()->get(&psgk[0], instance.m_n, 0);
      if (shift.front() == 0)
        for (size_t i = 0; i < instance.m_n; i++)
          psck[i] = psgk[i] * instance.m_F(i, i);
      else if (append) {
        psc.front()->get(&psck[0], instance.m_n, 0);
        //          molpro::cout << "resolvent action append "<<shift.front()<<shift.front()-1<<std::endl;
        //        molpro::cout << "initial psc="<<psc<<std::endl;
        for (size_t i = 0; i < instance.m_n; i++)
          if (i != instance.m_reference)
            psck[i] -= psgk[i] / (instance.m_F(i, i) + shift.front());
      } else {
        //          molpro::cout << "resolvent action replace "<<shift.front()<<std::endl;
        for (size_t i = 0; i < instance.m_n; i++)
          psck[i] = -psgk[i] / (instance.m_F(i, i) + shift.front());
        psck[instance.m_reference] = 0;
      }
      psc.front()->put(&psck[0], instance.m_n, 0);
      //        molpro::cout << "preconditioner output="<<psc<<std::endl;
    }
  } _rsptpot_updater;

  int nfail = 0;
  unsigned int iterations = 0, maxIterations = 0;
  size_t sample = 1;
  for (size_t repeat = 0; repeat < sample; repeat++) {
    instance.set(n, alpha);
    LinearEigensystem<ptype> d;
    d.m_verbosity = -1;
    d.m_rspt = true;
    d.m_minIterations = 50;
    d.m_thresh = 1e-5;
    d.m_maxIterations = 1000;
    //      ptype gg(n);
    vectorSet g;
    g.push_back(std::make_shared<ptype>(n));
    //      ptype xx=instance.guess();
    vectorSet x;
    x.push_back(std::make_shared<ptype>(instance.guess()));
    bool converged = false;
    for (int iteration = 1; (iteration < d.m_maxIterations && not converged) || iteration < d.m_minIterations;
         iteration++) {
      molpro::cout << "start of iteration " << iteration << std::endl;
      _rsptpot_residual(x, g);
      d.addVector(x, g);
      std::vector<scalar> shift;
      shift.push_back(1e-10);
      _rsptpot_updater(x, g, shift);
      converged = d.endIteration(x, g);
      molpro::cout << "end of iteration " << iteration << std::endl;
    }
    if (std::fabs(d.energy(d.m_minIterations) - d.eigenvalues().front()) > 1e-10)
      nfail++;
    molpro::cout << "Variational eigenvalue " << d.eigenvalues().front() << std::endl;
    for (size_t k = 0; k <= d.iterations(); k++) {
      molpro::cout << "E(" << k << ") = " << d.incremental_energies()[k] << ", cumulative=" << d.energy(k)
                   << ", error=" << d.energy(k) - d.eigenvalues()[0] << std::endl;
    }
    iterations += d.iterations();
    if (maxIterations < d.iterations())
      maxIterations = d.iterations();
  }
  molpro::cout << "sample=" << sample << ", n=" << n << ", alpha=" << alpha
               << ", average iterations=" << iterations / sample << ", maximum iterations=" << maxIterations
               << ", nfail=" << nfail << std::endl;
}
} // namespace linalg
} // namespace molpro

#ifdef ITERATIVESOLVER_FORTRAN
extern "C" {
void IterativeSolverFTest();
}
#endif
static std::unique_ptr<std::ofstream> out;
TEST(IterativeSolver_test, old) {
  if (true) {
    using namespace molpro::linalg;
    //  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
    //  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
    //  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
    //  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
    //  IterativeSolver::DIIS<double>::randomTest(100,100,0.1,3.0);
//    DIISTest<PagedArray<double>>(2, 6, 1e-10, 0.0002);
//    DIISTest<SimpleArray<double>>(2, 6, 1e-10, 0.0002);
    //  MPI_Abort(MPI_COMM_WORLD,1);
    //  DIISTest<LinearAlgebra::PagedVector<double> >(1,6,1e-10,IterativeSolver::DIIS<LinearAlgebra::PagedVector<double>
    //  >::DIISmode,0.2); DIISTest<LinearAlgebra::PagedVector<double>
    //  >(1,6,1e-3,IterativeSolver::DIIS<LinearAlgebra::PagedVector<double> >::disabled,0.0002);
    //   DavidsonTest<LinearAlgebra::PagedVector<double> >(2,2,2,2,false);
    if (true) {

      DavidsonTest<SimpleArray<double>>(3, 3, 1, 2);
//      DavidsonTest<PagedArray<double>>(3, 3, 1, 2);
      DavidsonTest<SimpleArray<double>>(3, 2, 1, 2);
//      DavidsonTest<PagedArray<double>>(3, 2, 1, 2);
      DavidsonTest<SimpleArray<double>>(9, 1, 1, 2);
//      DavidsonTest<PagedArray<double>>(9, 1, 1, 2);
      //      DavidsonTest<LinearAlgebra::PagedVector<double> >(9, 1, 1, 2, false);
      DavidsonTest<SimpleArray<double>>(9, 9, 1, 1);
//      DavidsonTest<PagedArray<double>>(9, 9, 1, 1);
      //      DavidsonTest<LinearAlgebra::PagedVector<double> >(9, 1, 1, 1, false);
      DavidsonTest<SimpleArray<double>>(9, 1, 1, 1);
//      DavidsonTest<PagedArray<double>>(9, 1, 1, 1);
      DavidsonTest<SimpleArray<double>>(9, 1, 1, 2);
//      DavidsonTest<PagedArray<double>>(9, 1, 1, 2);
      DavidsonTest<SimpleArray<double>>(9, 2, 1, 2);
//      DavidsonTest<PagedArray<double>>(9, 2, 1, 2);
      DavidsonTest<SimpleArray<double>>(100, 1, 1, 2);
//      DavidsonTest<PagedArray<double>>(100, 1, 1, 2);
      //      DavidsonTest<LinearAlgebra::PagedVector<double> >(100, 3, 1, 2, false);
//      DavidsonTest<PagedArray<double>>(100, 3, 1, 2);
      DavidsonTest<SimpleArray<double>>(100, 3, 1, 2);
      DavidsonTest<OpaqueArray<double>>(100, 3, 1, 2);
    }
//  DavidsonTest<LinearAlgebra::PagedVector<double> >(600,3,1,2,true);
//  RSPTTest<LinearAlgebra::PagedVector<double> ,double>(100,2e0);
#ifdef ITERATIVESOLVER_FORTRAN
    IterativeSolverFTest();
#endif
  }
}
