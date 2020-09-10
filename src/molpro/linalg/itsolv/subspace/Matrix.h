#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
#include <stddef.h>
#include <utility>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

/*!
 * @brief Matrix container that allows simple data access, slicing, copying and resizing without loosing data.
 * @tparam T
 */
template <typename T>
class Matrix {
protected:
  class Slice;
  class CSlice;

public:
  using value_type = T;
  using index_type = size_t;
  using coord_type = std::pair<size_t, size_t>;
  index_type rows() const {}
  index_type cols() const {}
  Slice slice(coord_type upper_left, coord_type bottom_right) {
    return Slice(*this, std::move(upper_left), std::move(bottom_right));
  }
  CSlice slice(coord_type upper_left, coord_type bottom_right) const {
    return CSlice(*this, std::move(upper_left), std::move(bottom_right));
  }
  void resize(const coord_type& row_col) {}
  coord_type dimension() const {}
  size_t size() const {}

protected:
  /*!
   * @brief Proxy mapping to a rectangular slice of the matrix data. Implements simple assignment.
   */
  class Slice {
  public:
    Slice(Matrix<T>& matrix, coord_type upper_left, coord_type bottom_right)
        : m(matrix), upl(std::move(upper_left)), btr(std::move(bottom_right)) {}
    Slice() = delete;
    virtual ~Slice() = default;
    Slice(const Slice&) = delete;
    Slice(Slice&&) = default;
    Slice& operator=(Slice&&) noexcept = default;

    Slice& operator=(const Slice& right) {}
    Slice& operator=(const CSlice& right) {}

    Slice& operator=(const Matrix<T>& right) {
      *this = right.slice({0, 0}, right.dimension());
      return *this;
    }

  protected:
    Matrix<T>& m;   //!< matrix being sliced
    coord_type upl; //!< upper left corner
    coord_type btr; //!< bottom right corner
  };

  //! Constant slice that cannot be assigned to
  class CSlice : public Slice {
  public:
    CSlice(const Matrix<T>& matrix, coord_type upper_left, coord_type bottom_right)
        : Slice(const_cast<Matrix<T>&>(matrix), std::move(upper_left), std::move(bottom_right)) {}
    CSlice() = delete;
    ~CSlice() override = default;
    CSlice(const CSlice&) = delete;
    CSlice(CSlice&&) = default;
    CSlice& operator=(CSlice&&) noexcept = default;
  };
};

} // namespace detail
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
