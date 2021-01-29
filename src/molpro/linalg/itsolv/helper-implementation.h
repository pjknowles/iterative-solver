#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <molpro/linalg/itsolv/helper.h>

namespace molpro::linalg::itsolv {

template <typename value_type>
int propose_singularity_deletion(size_t n, size_t ndim, const value_type* m, const std::vector<size_t>& candidates,
                                 double threshold) {
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> singularTester_(m, ndim, ndim);
  auto singularTester = singularTester_.block(0, 0, n, n);
  Eigen::JacobiSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> svd(singularTester, Eigen::ComputeThinV);
  //    molpro::cout << "propose_singularity_deletion threshold=" << threshold << std::endl;
  //        molpro::cout << "matrix:\n" << singularTester << std::endl;
  //        molpro::cout << "singular values:\n" << svd.singularValues().transpose() << std::endl;
  //        molpro::cout << "V:\n" << svd.matrixV() << std::endl;
  //        molpro::cout << "candidates:";
  //        for (const auto& c : candidates)
  //          molpro::cout << " " << c;
  //        molpro::cout << std::endl;
  auto sv = svd.singularValues();
  std::vector<double> svv;
  for (auto k = 0; k < n; k++)
    svv.push_back(sv(k));
  auto most_singular = std::min_element(svv.begin(), svv.end()) - svv.begin();
  //        molpro::cout << "most_singular " << most_singular << std::endl;
  if (svv[most_singular] > threshold)
    return -1;
  for (const auto& k : candidates) {
    //      if (std::fabs(svd.matrixV()(k, most_singular)) > 1e-3)
    //        molpro::cout << "taking candidate " << k << ": " << svd.matrixV()(k, most_singular) << std::endl;
    if (std::abs(svd.matrixV()(k, most_singular)) > 1e-3)
      return (int)k;
  }
  return -1;
}

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m, double threshold) {
  assert(m.size() == nrows * ncols);
  if (m.empty())
    return {};
  auto mat = Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(m.data(), nrows, ncols);
  auto svd = Eigen::JacobiSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>, Eigen::NoQRPreconditioner>(
      mat, Eigen::ComputeThinV);
  auto svd_system = std::list<SVD<value_type>>{};
  auto sv = svd.singularValues();
  for (int i = int(ncols) - 1; i >= 0; --i) {
    if (std::abs(sv(i)) < threshold) {
      auto t = SVD<value_type>{};
      t.value = sv(i);
      t.v.reserve(ncols);
      for (size_t j = 0; j < ncols; ++j) {
        t.v.emplace_back(svd.matrixV()(j, i));
      }
      svd_system.emplace_back(std::move(t));
    }
  }
  return svd_system;
}

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int>>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m, double threshold) {
  assert(false); // Complex not implemented here
  return {};
}

template <typename value_type>
void printMatrix(const std::vector<value_type>& m, size_t rows, size_t cols, std::string title, std::ostream& s) {
  s << title << "\n"
    << Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(m.data(), rows, cols) << std::endl;
}

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int>>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, const size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity) {
  assert(false); // Complex not implemented here
}

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, const size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity) {
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> HrowMajor(
      matrix.data(), dimension, dimension);
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> H(dimension, dimension);
  H = HrowMajor;
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> S(metric.data(), dimension, dimension);
  Eigen::MatrixXcd subspaceEigenvectors; // FIXME templating
  Eigen::VectorXcd subspaceEigenvalues;  // FIXME templating
  //   Eigen::GeneralizedEigenSolver<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> s(H, S);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(svdThreshold);
  //    molpro::cout << "singular values of overlap " << svd.singularValues().transpose() << std::endl;
  //   auto Hbar = svd.solve(H);
  if (verbosity > 1 && svd.rank() < S.cols())
    molpro::cout << "SVD rank " << svd.rank() << " in subspace of dimension " << S.cols() << std::endl;
  if (verbosity > 2 && svd.rank() < S.cols())
    molpro::cout << "singular values " << svd.singularValues().transpose() << std::endl;
  auto svmh = svd.singularValues().head(svd.rank()).eval();
  for (auto k = 0; k < svd.rank(); k++)
    svmh(k) = 1 / std::sqrt(svmh(k));
  auto Hbar = (svmh.asDiagonal()) * (svd.matrixU().leftCols(svd.rank()).adjoint()) * H *
              svd.matrixV().leftCols(svd.rank()) * (svmh.asDiagonal());
  //   molpro::cout << "S\n"<<S<<std::endl;
  //   molpro::cout << "S singular values"<<(Eigen::DiagonalMatrix<value_type, Eigen::Dynamic,
  //   Eigen::Dynamic>(svd.singularValues().head(svd.rank())))<<std::endl; molpro::cout << "S inverse singular
  //   values"<<Eigen::DiagonalMatrix<value_type,
  //   Eigen::Dynamic>(svd.singularValues().head(svd.rank())).inverse()<<std::endl; molpro::cout << "S singular
  //   values"<<sv<<std::endl; molpro::cout << "H\n"<<H<<std::endl; molpro::cout << "Hbar\n"<<Hbar<<std::endl;
  Eigen::EigenSolver<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> s(Hbar);
  //    molpro::cout << "s.eigenvectors()\n"<<s.eigenvectors()<<std::endl;
  subspaceEigenvalues = s.eigenvalues();
  if (s.eigenvalues().imag().norm() < 1e-10 and s.eigenvectors().imag().norm() < 1e-10) { // real eigenvectors
    subspaceEigenvectors = svd.matrixV().leftCols(svd.rank()) * svmh.asDiagonal() * s.eigenvectors().real();

  } else { // complex eigenvectors
#ifdef __INTEL_COMPILER
    molpro::cout << "Hbar\n" << Hbar << std::endl;
    molpro::cout << "Eigenvalues\n" << s.eigenvalues() << std::endl;
    molpro::cout << "Eigenvectors\n" << s.eigenvectors() << std::endl;
    throw std::runtime_error("Intel compiler does not support working with complex eigen3 entities properly");
#endif
    subspaceEigenvectors = svd.matrixV().leftCols(svd.rank()) * svmh.asDiagonal() * s.eigenvectors();
  }

  {
    // sort
    auto eigval = subspaceEigenvalues;
    auto eigvec = subspaceEigenvectors;
    std::vector<Eigen::Index> map;
    for (Eigen::Index k = 0; k < Hbar.cols(); k++) {
      Eigen::Index ll;
      for (ll = 0; std::count(map.begin(), map.end(), ll) != 0; ll++)
        ;
      for (Eigen::Index l = 0; l < Hbar.cols(); l++) {
        if (std::count(map.begin(), map.end(), l) == 0) {
          if (eigval(l).real() < eigval(ll).real())
            ll = l;
        }
      }
      map.push_back(ll);
      subspaceEigenvalues(k) = eigval(ll);
      //    molpro::cout << "new sorted eigenvalue "<<k<<", "<<ll<<", "<<eigval(ll)<<std::endl;
      //    molpro::cout << eigvec.col(ll)<<std::endl;
      subspaceEigenvectors.col(k) = eigvec.col(ll);
    }
  }
  //   molpro::cout << "sorted eigenvalues\n"<<subspaceEigenvalues<<std::endl;
  //   molpro::cout << "sorted eigenvectors\n"<<subspaceEigenvectors<<std::endl;
  if (!hermitian) {
    Eigen::MatrixXcd ovlTimesVec(subspaceEigenvectors.cols(), subspaceEigenvectors.rows()); // FIXME templating
    for (auto repeat = 0; repeat < 3; ++repeat)
      for (Eigen::Index k = 0; k < subspaceEigenvectors.cols(); k++) {
        if (std::abs(subspaceEigenvalues(k)) <
            1e-12) { // special case of zero eigenvalue -- make some real non-zero vector definitely in the null space
          subspaceEigenvectors.col(k).real() += double(0.3256897) * subspaceEigenvectors.col(k).imag();
          subspaceEigenvectors.col(k).imag().setZero();
        }
        if (hermitian)
          for (Eigen::Index l = 0; l < k; l++) {
            //        auto ovl =
            //            (subspaceEigenvectors.col(l).adjoint() * m_subspaceOverlap * subspaceEigenvectors.col(k))(
            //            0, 0); (ovlTimesVec.row(l) * subspaceEigenvectors.col(k))(0,0);
            //            ovlTimesVec.row(l).dot(subspaceEigenvectors.col(k));
            //        auto norm =
            //            (subspaceEigenvectors.col(l).adjoint() * subspaceOverlap * subspaceEigenvectors.col(l))(
            //                0,
            //                0);
            //      molpro::cout << "k="<<k<<", l="<<l<<", ovl="<<ovl<<" norm="<<norm<<std::endl;
            //      molpro::cout << subspaceEigenvectors.col(k).transpose()<<std::endl;
            //      molpro::cout << subspaceEigenvectors.col(l).transpose()<<std::endl;
            subspaceEigenvectors.col(k) -= subspaceEigenvectors.col(l) * // ovl;// / norm;
                                           ovlTimesVec.row(l).dot(subspaceEigenvectors.col(k));
            //        molpro::cout<<"immediately after projection " << k<<l<<" "<<
            //        (subspaceEigenvectors.col(l).adjoint() * subspaceOverlap * subspaceEigenvectors.col(k))( 0,
            //        0)<<std::endl;
          }
        //      for (Eigen::Index l = 0; l < k; l++) molpro::cout<<"after projection loop " << k<<l<<" "<<
        //      (subspaceEigenvectors.col(l).adjoint() * subspaceOverlap * subspaceEigenvectors.col(k))( 0,
        //      0)<<std::endl; molpro::cout <<
        //      "eigenvector"<<std::endl<<subspaceEigenvectors.col(k).adjoint()<<std::endl;
        auto ovl =
            //          (subspaceEigenvectors.col(k).adjoint() * subspaceOverlap *
            //          subspaceEigenvectors.col(k))(0,0);
            subspaceEigenvectors.col(k).adjoint().dot(S * subspaceEigenvectors.col(k));
        subspaceEigenvectors.col(k) /= std::sqrt(ovl.real());
        ovlTimesVec.row(k) = subspaceEigenvectors.col(k).adjoint() * S;
        //      for (Eigen::Index l = 0; l < k; l++)
        //      molpro::cout<<"after normalisation " << k<<l<<" "<< (subspaceEigenvectors.col(l).adjoint() *
        //      subspaceOverlap * subspaceEigenvectors.col(k))( 0, 0)<<std::endl; molpro::cout <<
        //      "eigenvector"<<std::endl<<subspaceEigenvectors.col(k).adjoint()<<std::endl;
        // phase
        Eigen::Index lmax = 0;
        for (Eigen::Index l = 0; l < subspaceEigenvectors.rows(); l++) {
          if (std::abs(subspaceEigenvectors(l, k)) > std::abs(subspaceEigenvectors(lmax, k)))
            lmax = l;
        }
        if (subspaceEigenvectors(lmax, k).real() < 0)
          subspaceEigenvectors.col(k) = -subspaceEigenvectors.col(k);
        //      for (Eigen::Index l = 0; l < k; l++)
        //      molpro::cout << k<<l<<" "<<
        //                       (subspaceEigenvectors.col(l).adjoint() * subspaceOverlap *
        //                       subspaceEigenvectors.col(k))( 0, 0)<<std::endl;
      }
  } // if (!hermitian)
  //     molpro::cout << "eigenvalues"<<std::endl<<subspaceEigenvalues<<std::endl;
  //     molpro::cout << "eigenvectors"<<std::endl<<subspaceEigenvectors<<std::endl;
  // TODO complex should be implemented with a specialised function
  // TODO real should be implemented with always-executed runtime assertion that eigensolution turns out to be real
  auto complex_error_vecs = (subspaceEigenvectors - subspaceEigenvectors.real()).norm();
  auto complex_error_vals = (subspaceEigenvalues - subspaceEigenvalues.real()).norm();
  assert(complex_error_vecs < 1e-12);
  assert(complex_error_vals < 1e-12);
  eigenvectors.resize(dimension * Hbar.cols());
  eigenvalues.resize(Hbar.cols());
  //    if constexpr (std::is_class<value_type>::value) {
  Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(eigenvectors.data(), dimension, Hbar.cols()) =
      subspaceEigenvectors.real();
  Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, 1>> ev(eigenvalues.data(), Hbar.cols());
  ev = subspaceEigenvalues.real();
  //    } else {
  //      Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(m_evec_xx.data(), dimension, dimension)
  //      =
  //          subspaceEigenvectors;
  //      Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, 1>>(eigenvalues.data(), dimension) = subspaceEigenvalues;
  //    }
}

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int>>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, const size_t dimension, size_t nroot,
                           double augmented_hessian, double svdThreshold, int verbosity) {
  assert(false); // Complex not implemented here
}

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
void solve_LinearEquations(std::vector<value_type>& solution, std::vector<value_type>& eigenvalues,
                           const std::vector<value_type>& matrix, const std::vector<value_type>& metric,
                           const std::vector<value_type>& rhs, const size_t dimension, size_t nroot,
                           double augmented_hessian, double svdThreshold, int verbosity) {
  const Eigen::Index nX = dimension;
  solution.resize(nX * nroot);
  //  std::cout << "augmented_hessian "<<augmented_hessian<<std::endl;
  if (augmented_hessian > 0) { // Augmented hessian
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> subspaceMatrix;
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> subspaceOverlap;
    subspaceMatrix.conservativeResize(nX + 1, nX + 1);
    subspaceOverlap.conservativeResize(nX + 1, nX + 1);
    subspaceMatrix.block(0, 0, nX, nX) =
        Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(matrix.data(), nX, nX);
    subspaceOverlap.block(0, 0, nX, nX) =
        Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(metric.data(), nX, nX);
    eigenvalues.resize(nroot);
    for (size_t root = 0; root < nroot; root++) {
      for (Eigen::Index i = 0; i < nX; i++) {
        subspaceMatrix(i, nX) = subspaceMatrix(nX, i) = -augmented_hessian * rhs[i + nX * root];
        subspaceOverlap(i, nX) = subspaceOverlap(nX, i) = 0;
      }
      subspaceMatrix(nX, nX) = 0;
      subspaceOverlap(nX, nX) = 1;
      //      std::cout << "subspace augmented hessian subspaceMatrix\n"<<subspaceMatrix<<std::endl;
      //      std::cout << "subspace augmented hessian subspaceOverlap\n"<<subspaceOverlap<<std::endl;

      Eigen::GeneralizedEigenSolver<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> s(subspaceMatrix,
                                                                                                 subspaceOverlap);
      auto eval = s.eigenvalues();
      auto evec = s.eigenvectors();
      Eigen::Index imax = 0;
      for (Eigen::Index i = 0; i < nX + 1; i++)
        if (eval(i).real() < eval(imax).real())
          imax = i;
      eigenvalues[root] = eval(imax).real();
      auto Solution = evec.col(imax).real().head(nX) / (augmented_hessian * evec.real()(nX, imax));
      for (auto k = 0; k < nX; k++)
        solution[k + nX * root] = Solution(k);
      //      std::cout << "subspace augmented hessian solution\n"<<Solution<<std::endl;
    }
  } else { // straight solution of linear equations
    Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> subspaceMatrixR(
        matrix.data(), nX, nX);
    Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> RHS_R(rhs.data(), nX,
                                                                                                       nroot);
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> subspaceMatrix = subspaceMatrixR;
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> RHS = RHS_R;
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> Solution;
    Solution = subspaceMatrix.householderQr().solve(RHS);
    //    std::cout << "subspace linear equations solution\n"<<Solution<<std::endl;
    for (size_t root = 0; root < nroot; root++)
      for (auto k = 0; k < nX; k++)
        solution[k + nX * root] = Solution(k, root);
  }
}

template <typename value_type>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, const size_t dimension,
                double svdThreshold, int verbosity) {
  auto nQ = dimension - 1;
  solution.resize(nQ + 1);
  if (nQ > 0) {
    //    Eigen::VectorXd Rhs(nQ), Coeffs(nQ);
    //    Eigen::MatrixXd B(nQ, nQ);
    Eigen::Matrix<value_type, Eigen::Dynamic, 1> Rhs(nQ), Coeffs(nQ);
    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> B(nQ, nQ);

    Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> subspaceMatrix(matrix.data(), nQ + 1,
                                                                                               nQ + 1);
    B.block(0, 0, nQ, nQ) = subspaceMatrix.block(0, 0, nQ, nQ);
    Rhs = -subspaceMatrix.block(0, nQ, nQ, 1);

    //    molpro::cout << "B:" << std::endl << B << std::endl;
    //    molpro::cout << "Rhs:" << std::endl << Rhs << std::endl;

    // invert the system, determine extrapolation coefficients.
    Eigen::JacobiSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> svd(B, Eigen::ComputeThinU |
                                                                                           Eigen::ComputeThinV);
    svd.setThreshold(svdThreshold);
    //    molpro::cout << "svdThreshold "<<svdThreshold<<std::endl;
    //    molpro::cout << "U\n"<<svd.matrixU()<<std::endl;
    //    molpro::cout << "V\n"<<svd.matrixV()<<std::endl;
    //    molpro::cout << "singularValues\n"<<svd.singularValues()<<std::endl;
    Coeffs = svd.solve(Rhs).head(nQ);
    if (verbosity > 1)
      molpro::cout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
    for (size_t k = 0; k < (size_t)Coeffs.rows(); k++) {
      if (std::isnan(std::abs(Coeffs(k)))) {
        molpro::cout << "B:" << std::endl << B << std::endl;
        molpro::cout << "Rhs:" << std::endl << Rhs << std::endl;
        molpro::cout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
        throw std::overflow_error("NaN detected in DIIS submatrix solution");
      }
      solution[k] = Coeffs(k);
    }
  }
  solution[nQ] = 1;
}
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
