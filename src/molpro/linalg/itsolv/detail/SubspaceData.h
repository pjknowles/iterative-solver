#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DETAIL_EQNDATA_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DETAIL_EQNDATA_H
#include <map>
#include <molpro/linalg/itsolv/detail/Matrix.h>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace detail {
enum class EqnData { H, S, rhs };

template <EqnData... DataTypes>
using SubspaceData = std::map<EqnData, Matrix<double>>;

template <EqnData... DataTypes>
auto null_data() {
  return SubspaceData<DataTypes...>{std::make_pair<EqnData, Matrix<double>>(DataTypes, {})...};
}
} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_DETAIL_EQNDATA_H
