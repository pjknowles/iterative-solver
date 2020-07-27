#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#include <Eigen/Dense>
#include <cmath>
#include <molpro/linalg/iterativesolver/helper.h>

template <typename scalar_type>
int molpro::linalg::iterativesolver::helper<scalar_type>::propose_singularity_deletion(
    size_t n, size_t ndim, const scalar_type* m, const std::vector<int>& candidates, double threshold) {
  Eigen::Map<const Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>> singularTester_(m, ndim, ndim);
  auto singularTester = singularTester_.block(0, 0, n, n);
  Eigen::JacobiSVD<Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>> svd(singularTester, Eigen::ComputeThinV);
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
      return k;
  }
  return -1;
}

template <typename scalar_type>
void molpro::linalg::iterativesolver::helper<scalar_type>::printMatrix(const std::vector<scalar_type> m, size_t rows,
                                                                       size_t cols, std::string title,
                                                                       std::ostream& s) {
  s << title << "\n"
    << Eigen::Map<const Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>>(m.data(), rows, cols) << std::endl;
}

template <>
void molpro::linalg::iterativesolver::helper<std::complex<double>>::eigenproblem(std::vector<std::complex<double>>& eigenvectors, std::vector<std::complex<double>>& eigenvalues, const std::vector<std::complex<double>> matrix, const std::vector<std::complex<double>> metric, const size_t dimension, bool hermitian, double svdThreshold, int verbosity) {}

template <typename T>
void molpro::linalg::iterativesolver::helper<T>::eigenproblem(std::vector<T>& eigenvectors, std::vector<T>& eigenvalues, const std::vector<T> matrix, const std::vector<T> metric, const size_t dimension, bool hermitian, double svdThreshold, int verbosity) {
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> H(matrix.data(), dimension, dimension);
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> S(metric.data(), dimension, dimension);
  Eigen::MatrixXcd subspaceEigenvectors; // FIXME templating
  Eigen::VectorXcd subspaceEigenvalues;  // FIXME templating
  //    molpro::cout << "diagonalizeSubspaceMatrix H:\n" << H.format(Eigen::FullPrecision) << std::endl;
  //    molpro::cout << "diagonalizeSubspaceMatrix S:\n" << S.format(Eigen::FullPrecision) << std::endl;
  //   Eigen::GeneralizedEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> s(H, S);
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
  //   molpro::cout << "S singular values"<<(Eigen::DiagonalMatrix<T, Eigen::Dynamic,
  //   Eigen::Dynamic>(svd.singularValues().head(svd.rank())))<<std::endl; molpro::cout << "S inverse singular
  //   values"<<Eigen::DiagonalMatrix<T, Eigen::Dynamic>(svd.singularValues().head(svd.rank())).inverse()<<std::endl;
  //   molpro::cout << "S singular values"<<sv<<std::endl;
  //   molpro::cout << "H\n"<<H<<std::endl;
  //   molpro::cout << "Hbar\n"<<Hbar<<std::endl;
  Eigen::EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> s(Hbar);
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
    std::vector<size_t> map;
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
      subspaceEigenvalues[k] = eigval(ll);
      //    molpro::cout << "new sorted eigenvalue "<<k<<", "<<ll<<", "<<eigval(ll)<<std::endl;
      //    molpro::cout << eigvec.col(ll)<<std::endl;
      subspaceEigenvectors.col(k) = eigvec.col(ll);
    }
  }
  //   molpro::cout << "sorted eigenvalues\n"<<subspaceEigenvalues<<std::endl;
  //   molpro::cout << "sorted eigenvectors\n"<<subspaceEigenvectors<<std::endl;
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
  //     molpro::cout << "eigenvalues"<<std::endl<<subspaceEigenvalues<<std::endl;
  //     molpro::cout << "eigenvectors"<<std::endl<<subspaceEigenvectors<<std::endl;
  eigenvectors.resize(dimension * dimension);
  eigenvalues.resize(dimension);
  // TODO complex should be implemented with a specialised function
  static_assert(not std::is_class_v<T>, "Complex not implemented here");
  // TODO real should be implemented with always-executed runtime assertion that eigensolution turns out to be real
  assert((subspaceEigenvectors - subspaceEigenvectors.real()).norm() < 1e-12);
  assert((subspaceEigenvalues - subspaceEigenvalues.real()).norm() < 1e-12);
  //    if constexpr (std::is_class<T>::value) {
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(eigenvectors.data(), dimension, dimension) =
      subspaceEigenvectors.real();
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> ev(eigenvalues.data(), dimension);
  ev = subspaceEigenvalues.real();
  //    } else {
  //      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(m_evec_xx.data(), dimension, dimension) =
  //          subspaceEigenvectors;
  //      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(eigenvalues.data(), dimension) = subspaceEigenvalues;
  //    }
}

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
