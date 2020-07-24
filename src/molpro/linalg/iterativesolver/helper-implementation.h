#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
#include <Eigen/Dense>
#include <molpro/linalg/iterativesolver/helper.h>
#include <cmath>

template <typename scalar_type>
int molpro::linalg::iterativesolver::helper<scalar_type>::propose_singularity_deletion(
    size_t n, size_t ndim, const scalar_type* m,
    const std::vector<int>& candidates, double threshold) {
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

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITERATIVESOLVER_HELPER_IMPLEMENTATION_H_
