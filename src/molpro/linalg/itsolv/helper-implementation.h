#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <molpro/Profiler.h>
#include <molpro/lapacke.h>
#include <molpro/linalg/itsolv/helper.h>

namespace molpro::linalg::itsolv {

template <typename value_type>
std::list<SVD<value_type>> svd_eigen_jacobi(size_t nrows, size_t ncols, const array::Span<value_type>& m,
                                            double threshold) {
  auto mat = Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(m.data(), nrows, ncols);
  auto svd = Eigen::JacobiSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>, Eigen::NoQRPreconditioner>(
      mat, Eigen::ComputeThinV);
  auto svd_system = std::list<SVD<value_type>>{};
  auto sv = svd.singularValues();
  for (int i = int(ncols) - 1; i >= 0; --i) {
    if (std::abs(sv(i)) < threshold) { // TODO: This seems to discard values ABOVE the threshold, not below it. it's
      auto t = SVD<value_type>{}; // also not scaling this threshold relative to the max singular value - find out why
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

template <typename value_type>
std::list<SVD<value_type>> svd_eigen_bdcsvd(size_t nrows, size_t ncols, const array::Span<value_type>& m,
                                            double threshold) {
  auto mat = Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(m.data(), nrows, ncols);
  auto svd = Eigen::BDCSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>>(mat, Eigen::ComputeThinV);
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

#if defined HAVE_CBLAS
template <typename value_type>
std::list<SVD<value_type>> svd_lapacke_dgesdd(size_t nrows, size_t ncols, const array::Span<value_type>& mat,
                                              double threshold) {
  int info;
  int m = nrows;
  int n = ncols;
  int sdim = std::min(m, n);
  std::vector<double> sv(sdim), u(nrows * nrows), v(ncols * ncols);
  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A', int(nrows), int(ncols), const_cast<double*>(mat.data()), int(ncols),
                        sv.data(), u.data(), int(nrows), v.data(), int(ncols));
  auto svd_system = std::list<SVD<value_type>>{};
  for (int i = int(ncols) - 1; i >= 0; --i) {
    if (std::abs(sv[i]) < threshold) {
      auto t = SVD<value_type>{};
      t.value = sv[i];
      t.v.reserve(ncols);
      for (size_t j = 0; j < ncols; ++j) {
        t.v.emplace_back(v[i * ncols + j]);
      }
      svd_system.emplace_back(std::move(t));
    }
  }
  return svd_system;
}

template <typename value_type>
std::list<SVD<value_type>> svd_lapacke_dgesvd(size_t nrows, size_t ncols, const array::Span<value_type>& mat,
                                              double threshold) {
  int info;
  int m = nrows;
  int n = ncols;
  int sdim = std::min(m, n);
  std::vector<double> sv(sdim), u(nrows * nrows), v(ncols * ncols);
  double superb[sdim - 1];
  info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'A', int(nrows), int(ncols), const_cast<double*>(mat.data()), int(ncols),
                        sv.data(), u.data(), int(nrows), v.data(), int(ncols), superb);
  auto svd_system = std::list<SVD<value_type>>{};
  for (int i = int(ncols) - 1; i >= 0; --i) {
    if (std::abs(sv[i]) < threshold) {
      auto t = SVD<value_type>{};
      t.value = sv[i];
      t.v.reserve(ncols);
      for (size_t j = 0; j < ncols; ++j) {
        t.v.emplace_back(v[i * ncols + j]);
      }
      svd_system.emplace_back(std::move(t));
    }
  }
  return svd_system;
}

#endif

#ifdef MOLPRO
extern "C" int dsyev_c(char, char, int, double*, int, double*);
#endif

/**
 * A wrapper function for lapacke_dsyev (linear eigensystem solver) from the lapack C interface (lapacke.h).
 * @param[in] matrix the input matrix (will not be altered!). Must be square. Must be dimension*dimension elements long.
 * @param[out] eigenvectors the matrix of eigenvectors, row major ordering. Must be square and the same size as matrix.
 * @param[out] eigenvalues a list of eigenvalues. Must be the length specified by the 'dimension' parameter.
 * @param[in] dimension length of one axis of the matrix.
 * \returns status. If 0, successful exit. If -i, the ith argument had an illegal value. If i, the algorithm failed to
 * converge.
 */
int eigensolver_lapacke_dsyev(const std::vector<double>& matrix, std::vector<double>& eigenvectors,
                              std::vector<double>& eigenvalues, const size_t dimension) {

  // validate input
  if (eigenvectors.size() != matrix.size()) {
    throw std::runtime_error("Matrix of eigenvectors and input matrix are not the same size!");
  }

  if (eigenvectors.size() != dimension * dimension || eigenvalues.size() != dimension) {
    throw std::runtime_error("Size of eigenvectors/eigenvlaues do not match dimension!");
  }

  // magic letters
  //  static const char compute_eigenvalues = 'N';
  static const char compute_eigenvalues_eigenvectors = 'V';
  //  static const char store_upper_triangle = 'U';
  static const char store_lower_triangle = 'L';

  // copy input matrix (lapack overwrites)
  std::copy(matrix.begin(), matrix.end(), eigenvectors.begin());

  // set lapack vars
  lapack_int status;
  lapack_int leading_dimension = dimension;
  lapack_int order = dimension;

  // call to lapack
#ifdef MOLPRO
  status = dsyev_c(compute_eigenvalues_eigenvectors, store_lower_triangle, order, eigenvectors.data(),
                   leading_dimension, eigenvalues.data());
#else
  status = LAPACKE_dsyev(LAPACK_COL_MAJOR, compute_eigenvalues_eigenvectors, store_lower_triangle, order,
                         eigenvectors.data(), leading_dimension, eigenvalues.data());
#endif

  return status;
}

/**
 * A wrapper function for lapacke_dsyev (linear eigensystem solver) from the lapack C interface (lapacke.h).
 * @param[in] matrix the input matrix (will not be altered!). Must be square. Must be dimension*dimension elements long.
 * @param[in] dimension length of one axis of the matrix.
 * \returns a std::list of instances of SVD, a struct containing one eigenvalue and one eigenvector. For a real,
 * symmetric matrix, these are equivalent to singular values, and S/D (which are both the same).
 */
std::list<SVD<double>> eigensolver_lapacke_dsyev(size_t dimension, std::vector<double>& matrix) {
  std::vector<double> eigvecs(dimension * dimension);
  std::vector<double> eigvals(dimension);

  // call to lapack
  int success = eigensolver_lapacke_dsyev(matrix, eigvecs, eigvals, dimension);
  if (success < 0) {
    throw std::invalid_argument("Invalid argument of eigensolver_lapacke_dsyev: ");
  }
  if (success > 0) {
    throw std::runtime_error("Lapacke_dsyev (eigensolver) failed to converge. "
                             " elements of an intermediate tridiagonal form did not converge to zero.");
  }

  auto eigensystem = std::list<SVD<double>>{};

  // populate eigensystem
  for (int i = dimension - 1; i >= 0;
       i--) { // note: flipping this axis gives parity with results of eigen::jacobiSVD
    auto temp_eigenproblem = SVD<double>{};
    temp_eigenproblem.value = eigvals[i];
    for (size_t j = 0; j < dimension; j++) {
      temp_eigenproblem.v.emplace_back(eigvecs[j + (dimension * i)]);
    }
    eigensystem.emplace_back(temp_eigenproblem);
  }

  return eigensystem;
}

/**
 * A wrapper function for lapacke_dsyev (linear eigensystem solver) from the lapack C interface (lapacke.h).
 * @param[in] dimension length of one axis of the matrix.
 * @param[in] matrix a span wrapping some data structure containing the elements of the matrix. Should be
 * dimension*dimension in length.
 * \returns a std::list of instances of SVD, a struct containing one eigenvalue and one eigenvector. For a real,
 * symmetric matrix, these are equivalent to singular values, and S/D (which are both the same).
 */
std::list<SVD<double>> eigensolver_lapacke_dsyev(size_t dimension,
                                                 const molpro::linalg::array::span::Span<double>& matrix) {
  // TODO: this should be the other way around, eigensolver_lapacke_dsyev should take a span by default and this should
  // wrap it with a vector
  std::vector<double> v;
  v.insert(v.begin(), matrix.begin(), matrix.end());
  return eigensolver_lapacke_dsyev(dimension, v);
}

/**
 * Get the rank of some matrix, given a threshold.
 * @param[in] eigenvalues the matrix, as a vector.
 * @param[in] threshold the threshold. Note that this is the normalised threshold, a value between 0 and 1, relative to
 * the largest element in the matrix.
 * \returns the rank. For an empty matrix, returns 0.
 */
template <typename value_type>
size_t get_rank(std::vector<value_type> eigenvalues, value_type threshold) {
  if (eigenvalues.size() == 0) {
    return 0;
  }
  value_type max = *max_element(eigenvalues.begin(), eigenvalues.end());
  value_type threshold_scaled = threshold * max;
  size_t count =
      std::count_if(eigenvalues.begin(), eigenvalues.end(), [&](auto const& val) { return val >= threshold_scaled; });
  return count;
}

/**
 * Get the rank of some matrix, given a threshold.
 * @param[in] svd_system a std::list containing SVD objects, such as those created by itsolv::svd_system
 * @param[in] threshold the threshold. Note that this is the normalised threshold, a value between 0 and 1, relative to
 * the largest element in the matrix.
 * \returns the rank. For an empty matrix, returns 0.
 */
template <typename value_type>
size_t get_rank(std::list<SVD<value_type>> svd_system, value_type threshold) {
  // compute max
  value_type max_value = 0;
  std::list<SVD<double>>::iterator it;
  for (it = svd_system.begin(); it != svd_system.end(); it++) {
    if (it->value > max_value) {
      max_value = it->value;
    }
  }
  // scale threshold
  value_type threshold_scaled = threshold * max_value;

  size_t rank = 0;
  // get rank
  for (it = svd_system.begin(); it != svd_system.end(); it++) {
    if (it->value > threshold_scaled) {
      rank += 1;
    }
  }
  return rank;
}

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m, double threshold,
                                      bool hermitian, bool reduce_to_rank) {
  std::list<SVD<value_type>> svds;
  assert(m.size() == nrows * ncols);
  if (m.empty())
    return {};
  if (hermitian) {
    assert(nrows == ncols);
    svds = eigensolver_lapacke_dsyev(nrows, m);
    for (auto s = svds.begin(); s != svds.end();)
      if (s->value > threshold)
        s = svds.erase(s);
      else
        ++s;
  } else {
    //#if defined HAVE_LAPACKE
    //    if (nrows > 16)
    //      return svd_lapacke_dgesdd<value_type>(nrows, ncols, m, threshold);
    //    return svd_lapacke_dgesvd<value_type>(nrows, ncols, m, threshold);
    //#endif
    svds = svd_eigen_jacobi<value_type>(nrows, ncols, m, threshold);
    // return svd_eigen_bdcsvd<value_type>(nrows, ncols, m, threshold);
  }

  // reduce to rank
  if (reduce_to_rank) {
    int rank = get_rank(svds, threshold);
    for (int i = ncols; i > rank; i--) {
      svds.pop_back();
    }
  }
  return svds;
}

template <typename value_type, typename std::enable_if_t<is_complex<value_type>{}, int>>
std::list<SVD<value_type>> svd_system(size_t nrows, size_t ncols, const array::Span<value_type>& m, double threshold,
                                      bool hermitian, bool reduce_to_rank) {
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
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity, bool condone_complex) {
  assert(false); // Complex not implemented here
}

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
void eigenproblem(std::vector<value_type>& eigenvectors, std::vector<value_type>& eigenvalues,
                  const std::vector<value_type>& matrix, const std::vector<value_type>& metric, size_t dimension,
                  bool hermitian, double svdThreshold, int verbosity, bool condone_complex) {
  auto prof = molpro::Profiler::single();
  prof->start("itsolv::eigenproblem");
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> HrowMajor(
      matrix.data(), dimension, dimension);
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> H(dimension, dimension);
  H = HrowMajor;
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> S(metric.data(), dimension, dimension);
  Eigen::MatrixXcd subspaceEigenvectors; // FIXME templating
  Eigen::VectorXcd subspaceEigenvalues;  // FIXME templating
  // Eigen::GeneralizedEigenSolver<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> s(H, S);

  // initialisation of variables
  Eigen::VectorXd singularValues;
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixV;
  Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixU;
  std::vector<double> eigvecs;
  std::vector<double> eigvals;
  int rank;

  // if the matrix is hermitian, we can use lapacke_dsyev
  if (hermitian) {
    eigvecs.resize(dimension * dimension);
    eigvals.resize(dimension);
    int success = eigensolver_lapacke_dsyev(metric, eigvecs, eigvals, dimension);
    if (success != 0) {
      throw std::runtime_error("Eigensolver did not converge");
    }
    singularValues = Eigen::Map<Eigen::VectorXd>(eigvals.data(), dimension);
    matrixV = Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(
        eigvecs.data(), dimension, dimension);
    matrixU = Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(
        eigvecs.data(), dimension, dimension);
    rank = get_rank(eigvals, svdThreshold);
  } else {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    singularValues = svd.singularValues();
    matrixV = svd.matrixV();
    matrixU = svd.matrixU();
    rank = svd.rank();
  }

  // svd.setThreshold(svdThreshold);
  //     molpro::cout << "singular values of overlap " << svd.singularValues().transpose() << std::endl;
  //    auto Hbar = svd.solve(H);
  if (verbosity > 1 && rank < S.cols())
    molpro::cout << "SVD rank " << rank << " in subspace of dimension " << S.cols() << std::endl;
  if (verbosity > 2 && rank < S.cols())
    molpro::cout << "singular values " << singularValues.transpose() << std::endl;
  auto svmh = singularValues.head(rank);
  for (auto k = 0; k < rank; k++)
    svmh(k) = svmh(k) > 1e-14 ? 1 / std::sqrt(svmh(k)) : 0;
  auto Hbar =
      (svmh.asDiagonal()) * (matrixU.leftCols(rank).adjoint()) * H * matrixV.leftCols(rank) * (svmh.asDiagonal());
  // std::cout << "\n\nHbar: \n" << Hbar << "\n\n";
  //    molpro::cout << "S\n"<<S<<std::endl;
  //    molpro::cout << "S singular values"<<(Eigen::DiagonalMatrix<value_type, Eigen::Dynamic,
  //    Eigen::Dynamic>(svd.singularValues().head(svd.rank())))<<std::endl; molpro::cout << "S inverse singular
  //    values"<<Eigen::DiagonalMatrix<value_type,
  //    Eigen::Dynamic>(svd.singularValues().head(svd.rank())).inverse()<<std::endl; molpro::cout << "S singular
  //    values"<<sv<<std::endl; molpro::cout << "H\n"<<H<<std::endl; molpro::cout << "Hbar\n"<<Hbar<<std::endl;
  Eigen::EigenSolver<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> s(Hbar);
  //      molpro::cout << "s.eigenvectors()\n"<<s.eigenvectors()<<std::endl;
  subspaceEigenvalues = s.eigenvalues();
  if (s.eigenvalues().imag().norm() < 1e-10) { // real eigenvalues
                                               //    molpro::cout << "eigenvalues near-enough real" << std::endl;
    subspaceEigenvalues = subspaceEigenvalues.real();
    subspaceEigenvectors = s.eigenvectors();
    // complex eigenvectors need to be rotated
    // assume that they come in consecutive pairs
    for (int i = 0; i < subspaceEigenvectors.cols(); i++) {
      if (subspaceEigenvectors.col(i).imag().norm() > 1e-10) {
        int j = i + 1;
        if (std::abs(subspaceEigenvalues(i) - subspaceEigenvalues(j)) < 1e-10 and
            subspaceEigenvectors.col(j).imag().norm() > 1e-10) {
          subspaceEigenvectors.col(j) = subspaceEigenvectors.col(i).imag() / subspaceEigenvectors.col(i).imag().norm();
          subspaceEigenvectors.col(i) = subspaceEigenvectors.col(i).real() / subspaceEigenvectors.col(i).real().norm();
        }
      }
    }
    subspaceEigenvectors = matrixV.leftCols(rank) * svmh.asDiagonal() * subspaceEigenvectors;
  } else { // complex eigenvectors
//    molpro::cout << "eigenvalues not near-enough real"<<std::endl;
//    molpro::cout << "s.eigenvalues() "<< s.eigenvalues().transpose()<<std::endl;
#ifdef __INTEL_COMPILER
    molpro::cout << "Hbar\n" << Hbar << std::endl;
    molpro::cout << "Eigenvalues\n" << s.eigenvalues() << std::endl;
    molpro::cout << "Eigenvectors\n" << s.eigenvectors() << std::endl;
    throw std::runtime_error("Intel compiler does not support working with complex eigen3 entities properly");
#endif
    subspaceEigenvectors = matrixV.leftCols(rank) * svmh.asDiagonal() * s.eigenvectors();
    //    std::cout << "subspaceEigenvectors\n" << subspaceEigenvectors << std::endl;
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
      double maxcomp =0;
      for (Eigen::Index l = 0; l < Hbar.cols(); l++) {
          if (std::abs(subspaceEigenvectors.col(k)[l].real()) > std::abs(subspaceEigenvectors.col(k)[maxcomp].real()))
              maxcomp = l;
      }
      if (subspaceEigenvectors.col(k)[maxcomp].real() < 0)
          subspaceEigenvectors.col(k) = - subspaceEigenvectors.col(k);
    }
  }

  // TODO: Need to address the case of near-zero eigenvalues (as below for non-hermitian case) and clean-up
  //  non-hermitian case

  //   molpro::cout << "sorted eigenvalues\n"<<subspaceEigenvalues<<std::endl;
  //   molpro::cout << "sorted eigenvectors\n"<<subspaceEigenvectors<<std::endl;
  //  molpro::cout << "hermitian="<<hermitian<<std::endl;
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
  }
  //   if (!hermitian) {
  //     molpro::cout << "eigenvalues"<<std::endl<<subspaceEigenvalues<<std::endl;
  //     molpro::cout << "eigenvectors" << std::endl << subspaceEigenvectors << std::endl;
  //   }
  if (condone_complex) {
    for (Eigen::Index root = 0; root < Hbar.cols(); ++root) {
      if (subspaceEigenvalues(root).imag() != 0) {
        //         molpro::cout << "complex eigenvalues: " << subspaceEigenvalues(root) << ", " <<
        //         subspaceEigenvalues(root + 1)
        //                      << std::endl;
        subspaceEigenvalues(root) = subspaceEigenvalues(root + 1) = subspaceEigenvalues(root).real();
        subspaceEigenvectors.col(root) = subspaceEigenvectors.col(root).real();
        subspaceEigenvectors.col(root + 1) = subspaceEigenvectors.col(root + 1).imag();
        ++root;
      }
    }
  }
  if ((subspaceEigenvectors - subspaceEigenvectors.real()).norm() > 1e-10 or
      (subspaceEigenvalues - subspaceEigenvalues.real()).norm() > 1e-10) {
    throw std::runtime_error("unexpected complex solution found");
  }
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
  prof->stop();
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

template <typename value_type, typename std::enable_if_t<!is_complex<value_type>{}, std::nullptr_t>>
void solve_DIIS(std::vector<value_type>& solution, const std::vector<value_type>& matrix, const size_t dimension,
                double svdThreshold, int verbosity) {
  auto nAug = dimension + 1;
  //  auto nQ = dimension - 1;
  solution.resize(dimension);
  //  if (nQ > 0) {
  Eigen::VectorXd Rhs(nAug), Coeffs(nAug);
  Eigen::MatrixXd BAug(nAug, nAug);
  //    Eigen::Matrix<value_type, Eigen::Dynamic, 1> Rhs(nQ), Coeffs(nQ);
  //    Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> B(nQ, nQ);
  //
  Eigen::Map<const Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> subspaceMatrix(matrix.data(), dimension,
                                                                                             dimension);
  BAug.block(0, 0, dimension, dimension) = subspaceMatrix;
  for (size_t i = 0; i < dimension; ++i) {
    BAug(dimension, i) = BAug(i, dimension) = -1;
    Rhs(i) = 0;
  }
  BAug(dimension, dimension) = 0;
  Rhs(dimension) = -1;
  //
  //        molpro::cout << "BAug:" << std::endl << BAug << std::endl;
  //        molpro::cout << "Rhs:" << std::endl << Rhs << std::endl;

  // invert the system, determine extrapolation coefficients.
  Eigen::JacobiSVD<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> svd(BAug, Eigen::ComputeThinU |
                                                                                            Eigen::ComputeThinV);

  //  std::cout << "svd thresholds " << svdThreshold << "," << svd.singularValues().maxCoeff() << std::endl;
  //  std::cout << "singular values " << svd.singularValues().transpose() << std::endl;
  svd.setThreshold(svdThreshold * svd.singularValues().maxCoeff() * 0);
  //    molpro::cout << "svdThreshold "<<svdThreshold<<std::endl;
  //    molpro::cout << "U\n"<<svd.matrixU()<<std::endl;
  //    molpro::cout << "V\n"<<svd.matrixV()<<std::endl;
  //    molpro::cout << "singularValues\n"<<svd.singularValues()<<std::endl;
  Coeffs = svd.solve(Rhs).head(dimension);
  //  Coeffs = BAug.fullPivHouseholderQr().solve(Rhs);
  //  molpro::cout << "Coeffs "<<Coeffs.transpose()<<std::endl;
  if (verbosity > 1)
    molpro::cout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
  for (size_t k = 0; k < (size_t)Coeffs.rows(); k++) {
    if (std::isnan(std::abs(Coeffs(k)))) {
      molpro::cout << "B:" << std::endl << BAug << std::endl;
      molpro::cout << "Rhs:" << std::endl << Rhs << std::endl;
      molpro::cout << "Combination of iteration vectors: " << Coeffs.transpose() << std::endl;
      throw std::overflow_error("NaN detected in DIIS submatrix solution");
    }
    solution[k] = Coeffs(k);
  }
}
} // namespace molpro::linalg::itsolv

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
