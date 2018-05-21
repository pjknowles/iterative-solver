#include <ctime>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "IterativeSolver/ISDiis.h"
#include "IterativeSolver.h"
#include "PagedVector.h"
namespace LinearAlgebra{
  /*!
   * \brief Test iterative solution of linear eigenvalue problem
   * \param dimension The dimension of the test matrix
   * \param roots How many eigensolutions to find
   * \param verbosity How much to report
   * \param problem Selects which test matrix to use
   * \param orthogonalize Whether to orthogonalize expansion vectors
   * \tparam ptype Concrete class template that implements LinearAlgebra::vectorSet
   * \tparam scalar Type of matrix elements
   */
 template <class ptype, class scalar=double>
  static void DavidsonTest(size_t dimension, size_t roots=1, int verbosity=0, int problem=0, bool orthogonalize=true)
  {

   static Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> testmatrix;

   static struct {
    void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs) const {
     for (size_t k=0; k<psx.size(); k++) {
      Eigen::VectorXd x(testmatrix.rows());
      if (psx[k]->size() != (size_t)testmatrix.rows()) throw std::logic_error("psx wrong size");
      psx[k]->get(&x[0],testmatrix.rows(),0);
      Eigen::VectorXd res = testmatrix * x;
      outputs[k]->put(&res[0],testmatrix.rows(),0);
     }
    }
   } action;

   static struct {
    void operator()(vectorSet<scalar> & psc, const vectorSet<scalar> & psg, std::vector<scalar> shift, bool append=true) const {
     size_t n=testmatrix.rows();
     std::vector<scalar> psck(n);
     std::vector<scalar> psgk(n);
     for (size_t k=0; k<psc.size(); k++) {
      psg[k]->get(&psgk[0],n,0);
      if (not append) psc.zero();
      psc[k]->get(&psck[0],n,0);
      for (size_t l=0; l<n; l++)  psck[l] -= psgk[l]/(testmatrix(l,l)+shift[k]);
      psc[k]->put(&psck[0],n,0);
     }
    }
   } update;

   xout << "Test IterativeSolver::Davidson dimension="<<dimension<<", roots="<<roots<<", problem="<<problem<<", orthogonalize="<<orthogonalize<<std::endl;
   testmatrix.resize(dimension,dimension);
   for (size_t k=0; k<dimension; k++)
    for (size_t l=0; l<dimension; l++)
     if (problem==0)
      testmatrix(l,k)=-1;
     else if (problem==1)
      testmatrix(l,k)=l+k+2;
     else if (problem==2)
      testmatrix(l,k)=( k==l ? k+1 : 1);
     else if (problem==3)
      testmatrix(l,k)=( k==l ? 1 : 1);
     else
      throw std::logic_error("invalid problem in Davidson::test");
   if (problem==3) testmatrix(0,1)=testmatrix(1,0)=1;

   LinearEigensystem<scalar> d;
   d.m_roots=roots;
   d.m_verbosity=verbosity;
   d.m_maxIterations=dimension;
   d.m_orthogonalize=orthogonalize;
   vectorSet<scalar> x;
   vectorSet<scalar> g;
   for (size_t root=0; root<(size_t)d.m_roots; root++) {
    auto xx=std::make_shared<ptype>(dimension);
    xx->zero();
//    xout << "zeroed xx="<<xx<<std::endl;
    scalar one=1; xx->put(&one,1,root);
//    xout << "constructed xx="<<xx<<std::endl;
    x.push_back_clone(xx);
    auto gg=std::make_shared<ptype>(dimension);
    g.push_back_clone(gg);
   }

   for (size_t iteration=0; iteration<dimension; iteration++) {
//      for (size_t kkk=0; kkk<x.size(); kkk++)
//          xout << "before action x: "<<x[kkk]<<std::endl;
//      for (size_t kkk=0; kkk<g.size(); kkk++)
//          xout << "before action g: "<<g[kkk]<<std::endl;
    action(x,g);
//      for (size_t kkk=0; kkk<x.size(); kkk++)
//          xout << "after action x: "<<x[kkk]<<std::endl;
//      for (size_t kkk=0; kkk<g.size(); kkk++)
//          xout << "after action g: "<<g[kkk]<<std::endl;
    d.addVector(x,g);
    std::vector<scalar> shift;
    for (size_t root=0; root<(size_t)d.m_roots; root++) shift.push_back(-d.eigenvalues()[root]+1e-14);
    update(x,g,shift);
    if (d.finalize(x,g)) break;
   }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>> es(testmatrix);
//    xout << "true eigenvalues: "<<es.eigenvalues().head(d.m_roots).transpose()<<std::endl;
//    xout << "true eigenvectors:\n"<<es.eigenvectors().leftCols(d.m_roots).transpose()<<std::endl;


   std::vector<scalar> ev=d.eigenvalues();
   xout << "Eigenvalues: "; size_t root=0; for (const auto& e : ev) xout<<" "<<e<<"(error="<<e-es.eigenvalues()(root++)<<")";xout<<std::endl;
   xout << "Reported errors: "; for (const auto& e: d.errors()) xout<<" "<<e;xout<<std::endl;

   action(x,g);
   std::vector<scalar> errors;
   for (size_t root=0; root<(size_t)d.m_roots; root++) {
    g[root]->axpy(-ev[root],*x[root]);
    errors.push_back(g[root]->dot(*g[root]));
   }
//   xout << "Square residual norms: "; for (typename std::vector<scalar>::const_iterator e=errors.begin(); e!=errors.end(); e++) xout<<" "<<*e;xout<<std::endl;
   xout << "Square residual norms: "; for (const auto& e: errors) xout<<" "<<e;xout<<std::endl;
   // be noisy about obvious problems
   if (*std::max_element(errors.begin(),errors.end())>1e-7) throw std::runtime_error("IterativeSolver::Davidson has failed tests");

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
  template<class ptype, class scalar=double>
    void DIISTest(int verbosity=0,
                     size_t maxDim=6,
                     double svdThreshold=1e-10,
                     enum DIIS<scalar>::DIISmode_type mode=DIIS<scalar>::DIISmode,
                     double difficulty=0.1)
{
static struct {
  void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs) const {
    size_t n=2;
    std::vector<scalar> psxk(n);
    std::vector<scalar> output(n);

    psx.front()->get(&(psxk[0]),n,0);
    output[0]=(2*psxk[0]-2+400*psxk[0]*(psxk[0]*psxk[0]-psxk[1])); output[1]=(200*(psxk[1]-psxk[0]*psxk[0])); // Rosenbrock
    outputs.front()->put(&(output[0]),n,0);
  }
} _Rosenbrock_residual;

static struct {
  void operator()(vectorSet<scalar> & psc, const vectorSet<scalar> & psg, std::vector<scalar> shift, bool append=true) const {
    size_t n=2;
    std::vector<scalar> psck(n);
    std::vector<scalar> psgk(n);
    psg.front()->get(&psgk[0],n,0);
    if (append) {
        psc.front()->get(&psck[0],n,0);
        psck[0] -=psgk[0]/700;
        psck[1] -=psgk[1]/200;
      } else
      {
        psck[0] =-psgk[0]/700;
        psck[1] =-psgk[1]/200;
      }
    psc.front()->put(&psck[0],n,0);
//    xout << "Rosenbrock updater, new psc="<<psc<<std::endl;
  }
} _Rosenbrock_updater;
  vectorSet<scalar> x; x.push_back(std::shared_ptr<ptype>(new ptype(2)));
  vectorSet<scalar> g; g.push_back(std::shared_ptr<ptype>(new ptype(2)));
  DIIS<scalar> d;
  d.m_maxDim=maxDim;
  d.m_svdThreshold=svdThreshold;
  d.setMode(mode);

  if (verbosity>=0) xout << "Test DIIS::iterate, difficulty="<<difficulty<<std::endl;
  d.Reset();
  d.m_verbosity=verbosity-1;
  std::vector<scalar> xxx(2);
  xxx[0]=xxx[1]=1-difficulty; // initial guess
  x.front()->put(&xxx[0],2,0);
  xout << "initial guess"<<x[0]<<std::endl;
  xout << "initial guess"<<x<<std::endl;
  bool converged=false;
  for (int iteration=1; iteration < 1000 && not converged; iteration++) {
   xout <<"start of iteration "<<iteration<<std::endl;
      _Rosenbrock_residual(x,g);
      optionMap o; //o["weight"]=2;
      d.addVector(x,g,o);
      std::vector<scalar> shift; shift.push_back(1e-10);
      _Rosenbrock_updater(x,g,shift);
      converged = d.finalize(x,g);
      x.front()->get(&xxx[0],2,0);
      if (verbosity>2)
        xout << "new x after iterate "<<x.front()<<std::endl;
      if (verbosity>=0)
        xout << "iteration "<<iteration<<", Residual norm = "<<std::sqrt(d.fLastResidual())
             << ", Distance from solution = "<<std::sqrt((xxx[0]-1)*(xxx[0]-1)+(xxx[1]-1)*(xxx[1]-1))
          <<", error = "<<d.errors().front()
            <<", converged? "<<converged
           <<std::endl;
   xout <<"end of iteration "<<iteration<<std::endl;
    }

  x.front()->get(&xxx[0],2,0);
  xout   << "Distance from solution = "<<std::sqrt((xxx[0]-1)*(xxx[0]-1)+(xxx[1]-1)*(xxx[1]-1))
           <<std::endl;

}

#include <cstdlib>
//struct anharmonic {
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
//    std::vector<scalar> r(m_n);
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

//static anharmonic instance;

//static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
//    std::vector<scalar> psxk(instance.m_n);
//    std::vector<scalar> output(instance.m_n);
//    psx.front()->get(&(psxk[0]),instance.m_n,0);
//    if (append)
//      outputs.front()->get(&(output[0]),instance.m_n,0);
//    else
//      outputs.front()->zero();

//    for (size_t i=0; i<instance.m_n; i++) {
//        output[i] = instance.m_gamma*psxk[i];
//        for (size_t j=0; j<instance.m_n; j++)
//          output[i] += instance.m_F(j,i)*psxk[j];
//      }
//    outputs.front()->put(&output[0],instance.m_n,0);
//  }
//} _anharmonic_residual;
//static struct : IterativeSolverBase::ParameterSetTransformation {
//  void operator()(const vectorSet<scalar> & psg, vectorSet<scalar> & psc, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const override {
//    std::vector<scalar> psck(instance.m_n);
//    std::vector<scalar> psgk(instance.m_n);
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
//void DIIS::randomTest(size_t sample, size_t n, double alpha, double gamma, DIISmode_type mode)
//{

//  int nfail=0;
//  unsigned int iterations=0, maxIterations=0;
//  for (size_t repeat=0; repeat < sample; repeat++) {
//      instance.set(n,alpha,gamma);
//      DIIS d(_anharmonic_residual,_anharmonic_preconditioner);
//      d.setMode(mode);
//      d.m_verbosity=-1;
//      d.m_maxIterations=100000;
//      ptype gg(n); vectorSet<scalar> g; g.push_back(std::shared_ptr<ptype>(&gg));
//      ptype xx=instance.guess(); vectorSet<scalar> x; x.push_back(std::shared_ptr<ptype>(&xx));
//      if (not d.solve(g,x)) nfail++;
//      iterations+=d.iterations();
//      if (maxIterations<d.iterations())
//        maxIterations=d.iterations();
//    }
//  xout << "sample="<<sample<<", n="<<n<<", alpha="<<alpha<<", gamma="<<gamma<<", average iterations="<<iterations/sample<<", maximum iterations="<<maxIterations<<", nfail="<<nfail<<std::endl;
//}


#include <cstdlib>
  template <class ptype, class scalar=double>
void RSPTTest(size_t n, double alpha)
{
  static struct rsptpot {
    Eigen::MatrixXd m_F;
    size_t m_n;
    rsptpot(){}
    size_t m_reference;
    void set(size_t n, double alpha)
    {
//        xout<<"rsptpot set"<<n<<std::endl;
      m_n=n;
      m_reference=0; // asserting that m_F(0,0) is the lowest

      m_F.resize(n,n);
      for (size_t j=0; j<n; j++) {
          for (size_t i=0; i<j; i++)
            m_F(i,j)=m_F(j,i)= -0.5 + (((double)rand())/RAND_MAX);
          m_F(j,j) = (j*alpha-1);
        }
//      xout << "m_F:"<<std::endl<<m_F<<std::endl;
    }
    ptype guess()
    {
      std::vector<scalar> r(m_n);
      ptype result(m_n);
      for (size_t k=0; k<m_n; k++)
        r[k]=0;
      r[m_reference]=1;
      result.put(&r[0],m_n,0);
      return result;
    }

  } instance;


    static struct {
    void operator()(const vectorSet<scalar> & psx, vectorSet<scalar> & outputs, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const {
//        xout << "rsptpot_residual"<<std::endl;
//        xout << "input "<<psx<<std::endl;
      std::vector<scalar> psxk(instance.m_n);
      std::vector<scalar> output(instance.m_n);
      psx.front()->get(&(psxk[0]),instance.m_n,0);
      if (append)
        outputs.front()->get(&(output[0]),instance.m_n,0);
      else
        outputs.front()->zero();
      for (size_t i=0; i<instance.m_n; i++) {
          output[i] = 0;
          for (size_t j=0; j<instance.m_n; j++) {
            output[i] += instance.m_F(j,i)*psxk[j];
          }
        }
        outputs.front()->put(&(output[0]),instance.m_n,0);
//        xout << "output "<<outputs<<std::endl;
    }
    } _rsptpot_residual;
    static struct {
    void operator()(vectorSet<scalar> & psc, const vectorSet<scalar> & psg, std::vector<scalar> shift=std::vector<scalar>(), bool append=false) const {
//        xout << "preconditioner input="<<psg<<std::endl;
//      if (shift.front()==0)
//          xout << "H0 not resolvent"<<std::endl;
      std::vector<scalar> psck(instance.m_n);
      std::vector<scalar> psgk(instance.m_n);
      psg.front()->get(&psgk[0],instance.m_n,0);
      if (shift.front()==0)
          for (size_t i=0; i<instance.m_n; i++)
            psck[i] = psgk[i]*instance.m_F(i,i);
      else if (append) {
          psc.front()->get(&psck[0],instance.m_n,0);
//          xout << "resolvent action append "<<shift.front()<<shift.front()-1<<std::endl;
//        xout << "initial psc="<<psc<<std::endl;
          for (size_t i=0; i<instance.m_n; i++)
              if (i != instance.m_reference)
                  psck[i] -= psgk[i]/(instance.m_F(i,i)+shift.front());
        } else {
//          xout << "resolvent action replace "<<shift.front()<<std::endl;
          for (size_t i=0; i<instance.m_n; i++)
            psck[i] =- psgk[i]/(instance.m_F(i,i)+shift.front());
          psck[instance.m_reference]=0;
        }
      psc.front()->put(&psck[0],instance.m_n,0);
//        xout << "preconditioner output="<<psc<<std::endl;
    }
    } _rsptpot_updater;

  int nfail=0;
  unsigned int iterations=0, maxIterations=0;
  size_t sample=1;
  for (size_t repeat=0; repeat < sample; repeat++) {
      instance.set(n,alpha);
      RSPT<scalar> d;
      d.m_verbosity=-1;
      d.m_minIterations=50;
      d.m_thresh=1e-5;
      d.m_maxIterations=1000;
//      ptype gg(n);
      vectorSet<scalar> g; g.push_back(std::make_shared<ptype>(n));
//      ptype xx=instance.guess();
      vectorSet<scalar> x; x.push_back(std::make_shared<ptype>(instance.guess()));
 bool converged=false;
for (int iteration=1; (iteration < d.m_maxIterations && not converged) || iteration < d.m_minIterations; iteration++) {
   xout <<"start of iteration "<<iteration<<std::endl;
      _rsptpot_residual(x,g);
      d.addVector(x,g);
      std::vector<scalar> shift; shift.push_back(1e-10);
      _rsptpot_updater(x,g,shift);
      converged = d.finalize(x,g);
   xout <<"end of iteration "<<iteration<<std::endl;
    }
      if (std::fabs(d.energy(d.m_minIterations)-d.eigenvalues().front()) > 1e-10) nfail++;
      xout << "Variational eigenvalue "<<d.eigenvalues().front()<<std::endl;
      for (size_t k=0; k<=d.iterations(); k++) {
          xout << "E("<<k<<") = "<<d.incremental_energies()[k]<<", cumulative="<<d.energy(k)<<", error="<<d.energy(k)-d.eigenvalues()[0]<<std::endl;
      }
      iterations+=d.iterations();
      if (maxIterations<d.iterations())
        maxIterations=d.iterations();
    }
  xout << "sample="<<sample<<", n="<<n<<", alpha="<<alpha<<", average iterations="<<iterations/sample<<", maximum iterations="<<maxIterations<<", nfail="<<nfail<<std::endl;
}
}


extern "C" { void IterativeSolverFTest();}
int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (rank==0) std::cout << size<<" MPI process"<<(size>1?"es ":"")<<std::endl;
//    std::cout << "MPI process number "<<rank<<std::endl;
  }
#endif

  if (true) {
   using namespace LinearAlgebra;
//  IterativeSolver::DIIS::randomTest(100,100,0.1,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.2,0.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,1.0);
//  IterativeSolver::DIIS::randomTest(100,100,0.1,2.0);
//  LinearAlgebra::DIIS<double>::randomTest(100,100,0.1,3.0);
//  DIISTest<PagedVector<double> >(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.0002);
//  DIISTest<PagedVector<double> >(1,6,1e-10,LinearAlgebra::DIIS<double>::DIISmode,0.2);
//  DIISTest<PagedVector<double> >(1,6,1e-3,LinearAlgebra::DIIS<double>::disabled,0.0002);
  DavidsonTest<PagedVector<double> >(2,1,1,2,true);
  DavidsonTest<PagedVector<double> >(9,1,1,2,true);
  DavidsonTest<PagedVector<double> >(9,1,1,2,false);
  DavidsonTest<PagedVector<double> >(9,1,1,1);
  DavidsonTest<PagedVector<double> >(9,1,1,2);
  DavidsonTest<PagedVector<double> >(9,2,1,2);
  DavidsonTest<PagedVector<double> >(100,1,1,2);
  DavidsonTest<PagedVector<double> >(100,3,1,2);
  DavidsonTest<PagedVector<double> >(100,3,1,2,false);
  DavidsonTest<PagedVector<double> >(600,3,1,2,true);
//  RSPTTest<PagedVector<double> ,double>(100,2e0);
  IterativeSolverFTest();
    }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
