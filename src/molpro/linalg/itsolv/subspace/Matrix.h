#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace molpro {
namespace linalg {
namespace itsolv {
namespace subspace {

// FIXME Since we are only supporting a few underlying data types Matrix should not be header only
// FIXME Usage of Eigen for matrix would not be as catastrophic now that the headers are separated
/*!
 * @brief Matrix container that allows simple data access, slicing, copying and resizing without loosing data.
 * @note Row-major order.
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

public:
  explicit Matrix(coord_type dims) : m_rows(dims.first), m_cols(dims.second), m_buffer(m_rows * m_cols) {}
  Matrix() = default;
  ~Matrix() = default;
  Matrix(const Matrix<T>&) = default;
  Matrix(Matrix<T>&&) noexcept = default;
  Matrix<T>& operator=(const Matrix<T>&) = default;
  Matrix<T>& operator=(Matrix<T>&&) noexcept = default;

  //! Access element @param i row index @param j column index
  T& operator()(index_type i, index_type j) { return m_buffer[i * m_cols + j]; }

  //! Access element @param i row index @param j column index
  T operator()(index_type i, index_type j) const { return m_buffer[i * m_cols + j]; }

  //! Access the underlying data buffer
  const std::vector<T>& data() const { return m_buffer; }

  //! Sets all elements of matrix to value
  void fill(T value) { std::fill(begin(m_buffer), end(m_buffer), value); }

  //! Access a rectangular slice of the matrix.
  Slice slice(coord_type upper_left, coord_type bottom_right) {
    return Slice(*this, std::move(upper_left), std::move(bottom_right));
  }

  //! Access a constant rectangular slice of the matrix.
  CSlice slice(coord_type upper_left, coord_type bottom_right) const {
    return CSlice(*this, std::move(upper_left), std::move(bottom_right));
  }

  //! Resize the matrix. The old data is preserved and any new rows/cols are zeroed. @param dims new dimensions
  void resize(const coord_type& dims) {
    if (dims == dimensions())
      return;
    auto m = Matrix<T>(dims);
    auto upper_left = coord_type{0, 0};
    auto bottom_right = coord_type{std::min(rows(), m.rows()), std::min(cols(), m.cols())};
    slice(upper_left, bottom_right) = m.slice(upper_left, bottom_right);
    std::swap(*this, m);
  }

  index_type rows() const { return m_rows; }
  index_type cols() const { return m_cols; }
  coord_type dimensions() const { return coord_type{m_rows, m_cols}; }
  size_t size() const { return m_rows * m_cols; }

protected:
  index_type m_rows = 0;   //!< number of rows
  index_type m_cols = 0;   //!< number of columns
  std::vector<T> m_buffer; //!< data buffer

  /*!
   * @brief Proxy mapping to a rectangular slice of the matrix data. Implements simple assignment.
   */
  class Slice {
  public:
    Slice(Matrix<T>& matrix, coord_type upper_left, coord_type bottom_right)
        : mat(matrix), upl(std::move(upper_left)), btr(std::move(bottom_right)) {}
    Slice() = delete;
    ~Slice() = default;
    Slice(const Slice&) = delete;
    Slice(Slice&&) noexcept = default;
    Slice& operator=(Slice&&) noexcept = default;

    Slice& operator=(const Slice& right) {}
    Slice& operator=(const CSlice& right) {}

    Slice& operator=(const Matrix<T>& right) {
      *this = right.slice({0, 0}, right.dimensions());
      return *this;
    }
    const Matrix<T>& M() const { return mat; }

  protected:
    Matrix<T>& mat; //!< matrix being sliced
    coord_type upl; //!< upper left corner
    coord_type btr; //!< bottom right corner
  };

  //! Constant slice that cannot be assigned to
  class CSlice {
  public:
    CSlice(const Matrix<T>& matrix, coord_type upper_left, coord_type bottom_right)
        : m_slice(const_cast<Matrix<T>&>(matrix), std::move(upper_left), std::move(bottom_right)) {}
    CSlice() = delete;
    ~CSlice() = default;
    CSlice(const CSlice&) = delete;
    CSlice(CSlice&&) noexcept = default;
    CSlice& operator=(CSlice&&) noexcept = default;

    const Matrix<T>& M() const { return m_slice.mat(); }

  protected:
    Slice m_slice;
  };
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
