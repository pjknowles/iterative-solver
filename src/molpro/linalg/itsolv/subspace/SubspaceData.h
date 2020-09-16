#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACEDATA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACEDATA_H
#include <map>
#include <molpro/linalg/itsolv/subspace/Matrix.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {
enum class EqnData { H, S, rhs };

using SubspaceData = std::map<EqnData, Matrix<double>>;

template <EqnData... DataTypes>
auto null_data() {
  return SubspaceData{std::make_pair<EqnData, Matrix<double>>(DataTypes, {})...};
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_SUBSPACEDATA_H