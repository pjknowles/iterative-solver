#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERWITHCOPYFUNC_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERWITHCOPYFUNC_H
#include <functional>
#include <molpro/linalg/array/ArrayHandler.h>
#include <utility>

namespace molpro::linalg::array {

/*!
 * @brief Base for handlers whose left arrays are created by a function given by the user.
 *
 * Disk arrays may need information beyond the array being copied in order to be created, e.g. where to put the
 * file, so that construction of a new left array as a copy of a right array is delegated to copy_func.
 */
template <class AL, class AR = AL>
class ArrayHandlerWithCopyFunc : public ArrayHandler<AL, AR> {
public:
  using copy_func_type = std::function<AL(const AR &)>;

  /*!
   * @brief Constructor taking a copy function.
   * @param copy_func returns a new array, which is a copy of its argument
   */
  explicit ArrayHandlerWithCopyFunc(copy_func_type copy_func)
      : ArrayHandler<AL, AR>(), m_copy_func(std::move(copy_func)) {}

  AL copy(const AR &source) override {
    this->m_counter->copy++;
    return m_copy_func(source);
  };
  using ArrayHandler<AL, AR>::copy;

protected:
  copy_func_type m_copy_func; //!< function for making a copy
};

} // namespace molpro::linalg::array
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ARRAY_ARRAYHANDLERWITHCOPYFUNC_H
