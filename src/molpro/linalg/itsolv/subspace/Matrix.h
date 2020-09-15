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
 *
 * This is not meant to be a general matrix, or to be very efficient. This is just a utility for internal implementation
 * of IterativeSolver and it only supports operations that are necessary. Memory management is not done very
 * efficiently, because for our use subspace matrices are relatively small.
 *
 * @note Row-major order.
 *
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

  //! Converts index of 1D data to matrix coordinate
  coord_type to_coord(size_t ind) const {
    if (ind >= size())
      throw std::out_of_range("index is larger than size");
    auto row = ind / m_cols;
    auto col = ind - row * m_cols;
    return {row, col};
  }

  //! Sets all elements of matrix to value
  void fill(T value) { std::fill(begin(m_buffer), end(m_buffer), value); }

  /*!
   * @brief Access a rectangular slice of the matrix.
   * @param upper_left indices of upper left corner
   * @param bottom_right past the end indices of bottom right corner, so that dimensions of slice were
   *                     (bottom_right-upper_left)
   */
  Slice slice(coord_type upper_left, coord_type bottom_right) {
    return Slice(*this, std::move(upper_left), std::move(bottom_right));
  }

  //! Access the whole matrix as a slice
  Slice slice() { return slice({0, 0}, dimensions()); }

  /*!
   * @brief Access a constant rectangular slice of the matrix.
   * @param upper_left indices of upper left corner
   * @param bottom_right past the end indices of bottom right corner, so that dimensions of slice were
   *                     (bottom_right-upper_left)
   */
  CSlice slice(coord_type upper_left, coord_type bottom_right) const {
    return CSlice(*this, std::move(upper_left), std::move(bottom_right));
  }

  //! Access the whole matrix as a slice
  CSlice slice() const { return slice({0, 0}, dimensions()); }

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

  //! removes a row from the matrix @param row index of the row to remove
  void remove_row(index_type row) {
    if (row >= m_rows)
      throw std::runtime_error("row is out of range");
    auto m = Matrix<T>({m_rows - 1, m_cols});
    m.slice({0, 0}, {row, m_cols}) = slice({0, 0}, {row, m_cols});
    m.slice({row, 0}, m.dimensions()) = slice({row + 1, 0}, dimensions());
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
        : mat(matrix), upl(std::move(upper_left)), btr(std::move(bottom_right)) {
      if (upl.first > btr.first || upl.second > btr.second)
        throw std::runtime_error("specified incorrect corners");
      if (upl.first < 0 || upl.second < 0)
        throw std::runtime_error("slice is out of range");
      if (btr.first > matrix.rows() || btr.second > matrix.cols())
        throw std::runtime_error("slice is out of range");
    }
    Slice() = delete;
    ~Slice() = default;
    Slice(const Slice&) = delete;
    Slice(Slice&&) noexcept = default;
    Slice& operator=(Slice&&) noexcept = default;

    Slice& operator=(const Slice& right) {
      if (dimensions() != right.dimensions())
        throw std::runtime_error("attempting to copy slices of different dimensions");
      for (size_t i = 0; i < dimensions().first; ++i) {
        for (size_t j = 0; j < dimensions().second; ++j) {
          mat(upl.first + i, upl.second + j) = right.mat(right.upl.first + i, right.upl.second + j);
        }
      }
      return *this;
    }

    Slice& axpy(T a, const Slice& x) {
      if (dimensions() != x.dimensions())
        throw std::runtime_error("attempting to copy slices of different dimensions");
      for (size_t i = 0; i < dimensions().first; ++i) {
        for (size_t j = 0; j < dimensions().second; ++j) {
          mat(upl.first + i, upl.second + j) += a * x.mat(x.upl.first + i, x.upl.second + j);
        }
      }
      return *this;
    }

    Slice& axpy(T a, const CSlice& x) {
      axpy(a, x.m_slice);
      return *this;
    }

    Slice& operator=(const CSlice& right) {
      *this = right.m_slice;
      return *this;
    }

    Slice& operator=(const Matrix<T>& right) {
      *this = right.slice({0, 0}, right.dimensions());
      return *this;
    }

    coord_type dimensions() const { return {btr.first - upl.first, btr.second - upl.second}; }

  protected:
    Matrix<T>& mat; //!< matrix being sliced
    coord_type upl; //!< upper left corner
    coord_type btr; //!< bottom right corner
  };

  //! Constant slice that cannot be assigned to
  class CSlice {
    friend Slice;

  public:
    CSlice(const Matrix<T>& matrix, coord_type upper_left, coord_type bottom_right)
        : m_slice(const_cast<Matrix<T>&>(matrix), std::move(upper_left), std::move(bottom_right)) {}
    CSlice() = delete;
    ~CSlice() = default;
    CSlice(const CSlice&) = delete;
    CSlice(CSlice&&) noexcept = default;
    CSlice& operator=(CSlice&&) noexcept = default;

  protected:
    Slice m_slice;
  };
};

} // namespace subspace
} // namespace itsolv
} // namespace linalg
} // namespace molpro

#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_SUBSPACE_MATRIX_H
