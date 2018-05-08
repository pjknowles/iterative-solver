#ifndef RSPT_H
#define RSPT_H
#include "IterativeSolver.h"
#include <stdexcept>
#include <cmath>

namespace LinearAlgebra{

  /** @example RSPTexample.cpp */
  /*!
 * \brief A class that finds the lowest eigensolution of a matrix as a perturbation series
 *
 * Example of simplest use: @include RSPTexample.cpp
 *
 */
 template <class scalar>
  class RSPT : public IterativeSolverBase<scalar>
  {
  using IterativeSolverBase<scalar>::m_residuals;
  using IterativeSolverBase<scalar>::m_solutions;
  using IterativeSolverBase<scalar>::m_others;
  using IterativeSolverBase<scalar>::m_roots;
  using IterativeSolverBase<scalar>::m_linear;
  public:
  using IterativeSolverBase<scalar>::m_verbosity;
    /*!
     * \brief RSPT
   * \param PP The PP block of the matrix
     */
  RSPT( const Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>& PP=Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>(0,0) )
      : IterativeSolverBase<scalar>(PP)
    {
      this->m_linear = true;
      this->m_orthogonalize = false;
      this->m_minIterations=10; // Ensure at least this order
      m_roots=1;
    }
    static void test (size_t n, double alpha);
  protected:
    virtual void extrapolate(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
{
  size_t n=m_solutions.size();
  // on entry, solution contains |n-1> and residual contains H|n-1>, already stored in m_solutions & m_residuals.
  // on exit, incremental_energies[n] contains E_n = <0 |H-H_0-E_1|n-1> = <0|H|n-1>-(E_0+E_1)<0|n-1>.
  // eventually, Wigner 2n+1 rule should be implemented.
  // on exit, residual contains
  // |d> = -(H_0-E_0)|n>
  //     = (H-H_0-E_1)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  //     = (H-E_0-E_1)|n-1> - (H_0-E_0)|n-1> - sum_{k=0}^{n-2} E_{n-k} |k>.
  // on exit, solution contains 0.

  {
  Eigen::MatrixXcd oldeigenvalues=this->m_subspaceEigenvalues;
  this->diagonalizeSubspaceMatrix();
  if (std::isnan(this->m_subspaceEigenvalues(0).real())) this->m_subspaceEigenvalues=oldeigenvalues;
  }
  if (m_verbosity>1) xout << "Subspace eigenvalues"<<std::endl<<this->m_subspaceEigenvalues<<std::endl;
  if (m_verbosity>2) xout << "Subspace eigenvectors"<<std::endl<<this->m_subspaceEigenvectors<<std::endl;

//      xout << "|"<<n-1<<">: "<<solution;
//      xout << "H|"<<n-1<<">: "<<residual;
  if (n == 1) {
      // the preconditioner function must take special action when shift=0
      // and return H0.residual not -(H0-shift)^{-1}.residual
      std::vector<double> shift; shift.push_back(0);
      throw std::logic_error("RSPT coding not complete");
//      m_preconditionerFunction(solution,residual,shift,false);
//      xout << "solution="<<solution;
//      xout << "residual="<<residual;
      m_E0 = solution.front()->dot(residual.front());
  }
      solution.zero();
      residual.zero();
      residual.axpy(1,m_residuals.back());
  if (n == 1) {
      m_incremental_energies.resize(2);
      m_incremental_energies[0]=m_E0;
      m_incremental_energies[1]=this->m_subspaceMatrix(0,0)-m_incremental_energies[0];
      residual.axpy(-this->m_subspaceMatrix(0,0),m_solutions.back());
      this->m_updateShift.clear();this->m_updateShift.push_back(-(1+std::numeric_limits<double>::epsilon())*m_E0);
      if (m_verbosity>=0)
          xout << "E_0="<<m_E0<<std::endl << "E_1="<<m_incremental_energies[1]<<std::endl;
//      xout << "d(1)after incrementing solution"<<residual<<std::endl;
  } else {
      m_incremental_energies.push_back(this->m_subspaceMatrix(n-1,0)-this->m_subspaceMatrix(0,0)*this->m_subspaceOverlap(n-1,0));
      if (m_verbosity>=0)
          xout << "E_"<<n<<"="<<m_incremental_energies.back()<<", Eigenvalue="<<this->m_subspaceEigenvalues(0)<<std::endl;
//      xout << "d(n)=g(n-1)"<<residual<<std::endl;
      residual.axpy(1,m_lastH0mE0psi);
//      xout << "d(n)=g(n-1)+d(n-1)"<<residual<<std::endl;
//      xout << "H(0,0) "<<this->m_subspaceMatrix(0,0)<<std::endl;
//      xout << "c(n-1) "<<m_solutions.back()<<std::endl;
      residual.axpy(-this->m_subspaceMatrix(0,0),m_solutions.back());
//      xout << "d(n)=g(n-1)+d(n-1)-(E0+E1)c(n-1)"<<residual<<std::endl;
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
  m_lastH0mE0psi = residual; // we will need this in the next iteration // FIXME does this leak memory?
//  xout << "-(H0-E0)|"<<n<<">: "<<residual<<std::endl;
}
    virtual void extrapolate(vectorSet<scalar> & residual, vectorSet<scalar> & solution, const optionMap options=optionMap()) { vectorSet<scalar> other; extrapolate(residual,solution,other,options); }

  public:
    int m_order; ///< Up to what order of perturbation theory should the energy be obtained.
    std::vector<double> incremental_energies(size_t state=0) {return m_incremental_energies;} ///< The incremental energies order by order.
///> The total energy to a given order.
    double energy(size_t order, ///< the desired maximum order of perturbation theory
                  size_t state=0 ///< the desired state
                               )
{
    double result=0.0;
    for (size_t k=0; k<=order; k++)
        result += m_incremental_energies[k];
    return result;
}

  private:
    vectorSet<scalar> m_lastH0mE0psi;
    std::vector<double> m_incremental_energies;
    double m_E0;
  };

// testing code below here
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
      d.interpolate(x,g);
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


#endif // RSPT_H
