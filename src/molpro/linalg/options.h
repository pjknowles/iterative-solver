#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_OPTIONS_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_OPTIONS_H_
#include <memory>
#include <molpro/Options.h>

namespace molpro::linalg {

/*!
 * @brief Get the Options object associated with LinearAlgebra
 * @return The current option set. If none has been assigned, in the Molpro context the "LINEARALGEBRA" option set is
 * used, otherwise empty.
 */
const std::shared_ptr<const molpro::Options> options();
/*!
 * @brief Set the options for LinearAlgebra
 * @param options
 */
void set_options(const molpro::Options& options);
} // namespace molpro::linalg

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_OPTIONS_H_
