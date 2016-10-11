#include "RSPT.h"
using namespace IterativeSolver;

RSPT::RSPT(const ParameterSetTransformation residualFunction, const ParameterSetTransformation preconditionerFunction)
  : IterativeSolverBase(residualFunction, preconditionerFunction)
  , m_roots(1)
{
  m_linear = true;
  m_orthogonalize = false;
  m_singularity_shift=1e-50;
  m_minIterations=10; // Ensure at least this order
}

void RSPT::extrapolate(ParameterVectorSet & residual, ParameterVectorSet & solution, ParameterVectorSet & other, const optionMap options)
{
  size_t n=m_solutions.size();
  // on entry, solution contains |n-1> and residual contains H|n-1>, already stored in m_solutions & m_residuals.
  // on exit, incremental_energies[n] contains E_n = <0 |H-H_0-E_1|n-1> = <0|H|n-1>-(E_0+E_1)<0|n-1>.
  // eventually, Wigner 2n+1 rule should be implemented.
  // on exit, residual contains
  // |d> = -(H_0-E_0)|n>
  //     = (H-H_0-E_1)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  //     = (H-E_0-E_1)|n-1> - (H_0-E_0)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  // note that <0|d>=<0|(H-H_0-E_1)|n-1> - sum_{k=0}^{n-2} E_{n-k} <0|k>. NOT RELEVANT IN GENERAL CASE
  // on exit, solution contains 0.

  {
  Eigen::MatrixXcd oldeigenvalues=m_subspaceEigenvalues;
  diagonalizeSubspaceMatrix();
  if (isnan(m_subspaceEigenvalues(0).real())) m_subspaceEigenvalues=oldeigenvalues;
  }
  if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<m_subspaceEigenvalues<<std::endl;
  if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<m_subspaceEigenvectors<<std::endl;

//      xout << "|"<<n-1<<">: "<<solution;
//      xout << "H|"<<n-1<<">: "<<residual;
  if (n == 1) {
      std::vector<double> shift; shift.push_back(0);
      m_preconditionerFunction(solution,residual,shift,false);
//      xout << "solution="<<solution;
//      xout << "residual="<<residual;
      m_E0 = solution.front()->dot(residual.front());
      solution.zero();
      m_incremental_energies.resize(2);
      m_incremental_energies[0]=m_E0;
      m_incremental_energies[1]=m_subspaceMatrix(0,0)-m_incremental_energies[0];
      residual.zero();
      residual.axpy(1,m_residuals.back());
      residual.axpy(-m_subspaceMatrix(0,0),m_solutions.back());
      m_updateShift.clear();m_updateShift.push_back(-m_E0);
      xout << "E_0="<<m_E0<<std::endl;
      xout << "E_1="<<m_incremental_energies[1]<<std::endl;
//      xout << "d(1)after incrementing solution"<<residual<<std::endl;
  } else {
      m_incremental_energies.push_back(m_subspaceMatrix(n-1,0)-m_subspaceMatrix(0,0)*m_subspaceOverlap(n-1,0));
      xout << "E_"<<n<<"="<<m_incremental_energies.back()<<", Eigenvalue="<<m_subspaceEigenvalues(0)<<std::endl;
      residual.zero();
      residual.axpy(1,m_residuals.back());
//      xout << "d(n)=g(n-1)"<<residual<<std::endl;
      residual.axpy(1,m_lastH0mE0psi);
//      xout << "d(n)=g(n-1)+d(n-1)"<<residual<<std::endl;
//      xout << "H(0,0) "<<m_subspaceMatrix(0,0)<<std::endl;
//      xout << "c(n-1) "<<m_solutions.back()<<std::endl;
      residual.axpy(-m_subspaceMatrix(0,0),m_solutions.back());
//      xout << "d(n)=g(n-1)+d(n-1)-(E0+E1)c(n-1)"<<residual<<std::endl;
      solution.zero();
  // this is structured for multistate, but not thought about yet
  for (size_t kkk=0; kkk<residual.size(); kkk++) {
      size_t l=0;
      for (size_t ll=0; ll<n-1; ll++) {
          for (size_t lll=0; lll<m_solutions[ll].size(); lll++) {
                  residual[kkk]->axpy(-m_incremental_energies[n-ll],m_solutions[ll][lll]);
//      xout << "d(n)after incrementing solution"<<residual<<std::endl;
                  l++;
            }
        }
    }

  }
  m_updateShift.resize(m_roots);
  for (size_t root=0; root<(size_t)m_roots; root++) m_updateShift[root]=m_singularity_shift-m_subspaceMatrix(root,root);
  m_lastH0mE0psi = residual; // we will need this in the next iteration // FIXME does this leak memory?
//  xout << "-(H0-E0)|"<<n<<">: "<<residual<<std::endl;
}
double RSPT::energy(size_t order, size_t state)
{
    double result;
    for (size_t k=0; k<=order; k++)
        result += m_incremental_energies[k];
    return result;
}

#include "SimpleParameterVector.h"
#include <random>
#include <chrono>
  struct rsptpot {
    Eigen::MatrixXd m_F;
    size_t m_n;
    rsptpot(){}
    int m_reference;
    void set(size_t n, double alpha)
    {
//        xout<<"rsptpot set"<<n<<std::endl;
      m_n=n;
      m_reference=0; // asserting that m_F(0,0) is the lowest
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator (seed);
      std::uniform_real_distribution<double> distribution(-0.5,0.5);

      m_F.resize(n,n);
      for (size_t j=0; j<n; j++) {
          for (size_t i=0; i<j; i++)
            m_F(i,j)=m_F(j,i)=1;//distribution(generator);
          m_F(j,j) = (j*alpha-1);
        }
//      xout << "m_F:"<<std::endl<<m_F<<std::endl;
    }
    SimpleParameterVector guess()
    {
      SimpleParameterVector result(m_n);
      for (size_t k=0; k<m_n; k++) {
        result[k]=0;
        }
      result[m_reference]=1;
      return result;
    }
  };

  static rsptpot *instance;

    static void _rsptpot_residual(const ParameterVectorSet & psx, ParameterVectorSet & outputs, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
//        xout << "rsptpot_residual"<<std::endl;
//        xout << "input "<<psx<<std::endl;
      if (not append) outputs.front()->zero();
      for (size_t i=0; i<instance->m_n; i++) {
          (*outputs.front())[i] = 0;
          for (size_t j=0; j<instance->m_n; j++) {
            (*outputs.front())[i] += instance->m_F(j,i)*(*psx.front())[j];
//              xout << "add " <<instance->m_F(j,i)<<"*"<<(*psx.front())[j]<<std::endl;
          }
//          xout << (*outputs.front())[i]<<std::endl;
        }
//        xout << "output "<<outputs<<std::endl;
    }
    static void _rsptpot_preconditioner(const ParameterVectorSet & psg, ParameterVectorSet & psc, std::vector<ParameterScalar> shift=std::vector<ParameterScalar>(), bool append=false) {
//        xout << "preconditioner input="<<psg<<std::endl;
//      if (shift.front()==0)
//          xout << "H0 not resolvent"<<std::endl;
      if (shift.front()==0)
          for (size_t i=0; i<instance->m_n; i++)
            (*psc.front())[i] = (*psg.front())[i]*instance->m_F(i,i);
      else if (append) {
//          xout << "resolvent action append "<<shift.front()<<shift.front()-1<<std::endl;
//        xout << "initial psc="<<psc<<std::endl;
          for (size_t i=0; i<instance->m_n; i++)
              if (i != instance->m_reference)
                  (*psc.front())[i] -= (*psg.front())[i]/(instance->m_F(i,i)+shift.front());
        } else {
//          xout << "resolvent action replace "<<shift.front()<<std::endl;
          for (size_t i=0; i<instance->m_n; i++)
            (*psc.front())[i] =- (*psg.front())[i]/(instance->m_F(i,i)+shift.front());
          (*psc.front())[instance->m_reference]=0;
        }
//        xout << "preconditioner output="<<psc<<std::endl;
    }
void RSPT::test(size_t n, double alpha)
{
    instance = new rsptpot();
  int nfail=0;
  unsigned int iterations=0, maxIterations=0;
  size_t sample=1;
  for (size_t repeat=0; repeat < sample; repeat++) {
      instance->set(n,alpha);
      RSPT d(&_rsptpot_residual,&_rsptpot_preconditioner);
      d.m_verbosity=-1;
      d.m_roots=1;
      d.m_order=100;
      d.m_thresh=1e-5;
      d.m_maxIterations=88;
      SimpleParameterVector gg(n); ParameterVectorSet g; g.push_back(&gg);
      SimpleParameterVector xx=instance->guess();
      ParameterVectorSet x; x.push_back(&xx);
      d.solve(g,x);
      if (std::fabs(d.energy(d.m_order)-d.eigenvalues().front()) > 1e-10) nfail++;
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
