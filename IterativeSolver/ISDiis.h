#ifndef DIIS_H
#define DIIS_H

#include "IterativeSolver.h"

namespace LinearAlgebra {

  /** @example DIISexample.cpp */
  /*!
 * \brief A class that encapsulates accelerated convergence of non-linear equations
 * through the DIIS or related methods.
 *
 * Example of simplest use: @include DIISexample.cpp
 *
 */
 template <class scalar=double>
  class DIIS : public IterativeSolverBase<scalar>
  {
  using IterativeSolverBase<scalar>::m_residuals;
  using IterativeSolverBase<scalar>::m_solutions;
  using IterativeSolverBase<scalar>::m_others;
  public:
  using IterativeSolverBase<scalar>::m_verbosity;
    enum DIISmode_type {disabled ///< No extrapolation is performed
                        , DIISmode ///< Direct Inversion in the Iterative Subspace
                        , KAINmode ///< Krylov Accelerated Inexact Newton
                       };
    /*!
   * \brief DIIS
   * \param PP The PP block of the matrix
   */
  DIIS( const Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>& PP=Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic>(0,0) )
   : IterativeSolverBase<scalar>(PP)
  , m_svdThreshold(1e-10)
  , m_maxDim(6)
{
  this->m_orthogonalize = false;
  setMode(DIISmode);
  Reset();
}

    ~DIIS() { }


    /*!
   * \brief Set options for DIIS.
   * \param mode Whether to perform DIIS, KAIN, or nothing.
   */
    virtual void setMode(enum DIISmode_type mode=DIISmode)
    {
     m_DIISmode = mode;
     this->m_subspaceMatrixResRes = mode!=KAINmode;
//     this->m_preconditionResiduals = mode==KAINmode; // FIXME

     Reset();
     if (m_verbosity>1)
      xout << "m_DIISmode set to "<<m_DIISmode<<std::endl;
    }
    /*!
   * \brief discards previous iteration vectors, but does not clear records
   */
    void Reset()
{
  m_LastResidualNormSq=1e99; // so it can be tested even before extrapolation is done
  this->m_lastVectorIndex=0;
  while (this->m_subspaceMatrix.rows()>0) this->deleteVector(0); ;
  this->m_residuals.clear();
  this->m_solutions.clear();
  this->m_others.clear();
  this->m_Weights.clear();
}
    /*!
   * \brief Introduce a new iteration vector, and perform extrapolation
   * \param residual
   * The vector that
   * will be the one that is analysed to construct the extrapolation.
   * \param solution
   * The current solution that gave rise to residual, and which will be extrapolated to a new predicted solution.
   * \param other (optional)
   * Corresponding other vectors whose sequence will be extrapolated.
   * \param options can contain "weight=xxx" where xxx is the weight to be given to this vector. These options would normally be passed as the corresponding parameter in iterate().
   */
  protected:
    void extrapolate(vectorSet<scalar> & solution, vectorSet<scalar> & residual, vectorSet<scalar> & other, const optionMap options=optionMap())
{
  //	  xout << "Enter DIIS::extrapolate"<<std::endl;
  //	  xout << "residual : "<<residual<<std::endl;
  //	  xout << "solution : "<<solution<<std::endl;
  this->m_updateShift.clear();this->m_updateShift.push_back(-(1+std::numeric_limits<double>::epsilon())*this->m_subspaceMatrix(0,0));
  double weight=options.count("weight") ? (options.find("weight")->second) : 1.0;
  if (this->m_maxDim <= 1 || this->m_DIISmode == disabled) return;

  if (residual.size() > 1) throw std::logic_error("DIIS does not handle multiple solutions");
  this->m_lastVectorIndex=this->m_residuals.size()-1;

  //  if (m_subspaceMatrix.rows() < 9) {
  //      xout << "m_subspaceMatrix on entry to DIIS::extrapolate"<<std::endl<<m_subspaceMatrix<<std::endl;
  //  }
  size_t nDim = this->m_subspaceMatrix.rows();
  this->m_LastResidualNormSq = std::fabs(this->m_subspaceMatrix(nDim-1,nDim-1));
  //  xout << "this->m_LastResidualNormSq "<<this->m_LastResidualNormSq<<std::endl;

  this->m_Weights.push_back(weight);

  Eigen::Array<scalar,Eigen::Dynamic,Eigen::Dynamic> d = this->m_subspaceMatrix.diagonal().array().abs();
  int worst=0, best=0;
  for (size_t i=0; i<nDim; i++) {
      if (d(i) > d(worst)) worst=i;
      if (d(i) < d(best)) best=i;
    }
  double fBaseScale = std::sqrt(d(worst)*d(best));

  while (nDim > this->m_maxDim) { // prune away the worst/oldest vector. Algorithm to be done properly yet
      size_t prune = worst;
      if (true || prune == nDim-1) { // the current vector is the worst, so delete the oldest
          //          xout << "this->m_dateOfBirth: "; for (auto b=this->m_dateOfBirth.begin(); b!=this->m_dateOfBirth.end(); b++) xout <<(*b); xout<<std::endl;
          prune = std::min_element(this->m_dateOfBirth.begin(),this->m_dateOfBirth.end())-this->m_dateOfBirth.begin();
        }
      //      xout << "prune="<<prune<<std::endl;
      //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
      this->deleteVector(prune);
      m_Weights.erase(m_Weights.begin()+prune);
      nDim--;
      //  xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
    }
  if (nDim != (size_t)this->m_subspaceMatrix.rows()) throw std::logic_error("problem in pruning");
  if (m_Weights.size() != (size_t)this->m_subspaceMatrix.rows()) {
      xout << "nDim="<<nDim<<", m_Weights.size()="<<m_Weights.size()<<std::endl;
      throw std::logic_error("problem after pruning weights");
    }



  //  if (m_subspaceMatrix.rows() < 9) {
  //      xout << "m_subspaceMatrix"<<std::endl<<m_subspaceMatrix<<std::endl;
  //  }

  // build actual DIIS system for the subspace used.
  Eigen::VectorXd
      Rhs(nDim+1),
      Coeffs(nDim+1);
  Eigen::MatrixXd
      B(nDim+1, nDim+1);

  // Factor out common size scales from the residual dots.
  // This is done to increase numerical stability for the case when _all_
  // residuals are very small.
  B.block(0,0,nDim,nDim) = this->m_subspaceMatrix/fBaseScale;
  Rhs.head(nDim) = Eigen::VectorXd::Zero(nDim);

  // make Lagrange/constraint lines.
  for ( size_t i = 0; i < nDim; ++ i )
    B(i, nDim) =  B(nDim, i) = -m_Weights[i];
  B(nDim, nDim) = 0.0;
  Rhs[nDim] = -1;
  //  xout << "B:"<<std::endl<<B<<std::endl;
  //  xout << "Rhs:"<<std::endl<<Rhs<<std::endl;

  // invert the system, determine extrapolation coefficients.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(m_svdThreshold);
  Coeffs=svd.solve(Rhs).head(nDim);
  m_LastAmplitudeCoeff = Coeffs[nDim-1];
  if (m_verbosity>1) xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
  for (size_t k=0; k<(size_t)Coeffs.rows(); k++)
    if (std::isnan(Coeffs(k))) {
        xout << "B:"<<std::endl<<B<<std::endl;
        xout << "Rhs:"<<std::endl<<Rhs<<std::endl;
        xout << "Combination of iteration vectors: "<<Coeffs.transpose()<<std::endl;
        throw std::overflow_error("NaN detected in DIIS submatrix solution");
      }
  residual.zero();
  solution.zero();
  size_t k=0;
  for (size_t l=0; l<other.size(); l++) other[l]->zero();
  for (size_t kk=0; kk<m_residuals.size(); kk++) {
      if (m_residuals[kk].m_active.front()){
          residual.front()->axpy(Coeffs[k],*m_residuals[kk].front());
          solution.front()->axpy(Coeffs[k],*m_solutions[kk].front());
          for (size_t l=0; l<other.size(); l++)
            other[l]->axpy(Coeffs[k],*m_others[kk][l]);
          k++;
        }
    }
  if (m_verbosity>2) {
      xout << "DIIS.extrapolate() final extrapolated solution: "<<solution.front()<<std::endl;
      xout << "DIIS.extrapolate() final extrapolated residual: "<<residual.front()<<std::endl;
    }
}

    void extrapolate(vectorSet<scalar> & solution, vectorSet<scalar> & residual, const optionMap options=optionMap()) { vectorSet<scalar> other; extrapolate(solution,residual,other,options); }

    /*!
   * \brief Return the square L2 norm of the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
  public:
    double fLastResidual() const { return m_LastResidualNormSq; }
    /*!
   * \brief Return the coefficient of the last residual vector in the extrapolated residual from the last call to extrapolate() or iterate().
   * \return
   */
    double fLastCoeff() const { return m_LastAmplitudeCoeff; }
    unsigned int nLastDim() const {return m_residuals.size();}
    unsigned int nMaxDim() const { return m_maxDim; }
    double m_svdThreshold; ///< Threshold for singular-value truncation in linear equation solver.
    size_t m_maxDim; ///< Maximum DIIS dimension allowed.
    static void randomTest(size_t sample, size_t n=100, double alpha=0.1, double gamma=0.0, DIISmode_type mode=DIISmode);
  private:
    typedef unsigned int uint;
    enum DIISmode_type m_DIISmode;

    std::vector<double> m_Weights;

    // the following variables are kept for informative/displaying purposes
    double
    // dot(R,R) of last residual vector fed into this state.
    m_LastResidualNormSq,
    // coefficient the actual new vector got in the last DIIS step
    m_LastAmplitudeCoeff;

  };

  template <class scalar>
  class KAIN : public DIIS<scalar>
  {
  public:
    KAIN()
      : DIIS<scalar>() { setMode(this->KAINmode);}
  };
}
#endif // DIIS_H
