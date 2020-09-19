#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/subspace/Matrix.h>
#include <numeric>

using molpro::linalg::itsolv::subspace::Matrix;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Eq;
using ::testing::Ne;
using ::testing::Pointwise;

TEST(Matrix, default_constructor) {
  auto m = Matrix<double>();
  ASSERT_EQ(m.rows(), 0);
  ASSERT_EQ(m.cols(), 0);
  ASSERT_EQ(m.size(), 0);
  ASSERT_EQ(m.dimensions(), (Matrix<double>::coord_type{0, 0}));
}

TEST(Matrix, constructor) {
  const size_t r = 3, c = 5;
  auto m = Matrix<double>({r, c});
  ASSERT_EQ(m.rows(), r);
  ASSERT_EQ(m.cols(), c);
  ASSERT_EQ(m.size(), r * c);
  ASSERT_EQ(m.dimensions(), (Matrix<double>::coord_type{r, c}));
}

struct MatrixF : ::testing::Test {
  MatrixF() : m({r, c}){};

  void iota() {
    for (size_t i = 0, ij = 0; i < m.rows(); ++i)
      for (size_t j = 0; j < m.cols(); ++j, ++ij)
        m(i, j) = ij;
  }

  using coord_type = Matrix<double>::coord_type;
  const size_t r = 3;
  const size_t c = 4;
  Matrix<double> m;
};

TEST_F(MatrixF, to_coord) {
  ASSERT_EQ(m.to_coord(0), (coord_type{0, 0}));
  ASSERT_EQ(m.to_coord(1), (coord_type{0, 1}));
  ASSERT_EQ(m.to_coord(c), (coord_type{1, 0}));
  ASSERT_EQ(m.to_coord(c + 1), (coord_type{1, 1}));
  ASSERT_EQ(m.to_coord(2 * c + 2), (coord_type{2, 2}));
  ASSERT_EQ(m.to_coord(m.size() - 1), (coord_type{r - 1, c - 1}));
  ASSERT_THROW(m.to_coord(m.size()), std::out_of_range);
}

TEST_F(MatrixF, fill) {
  const double value = 3.14;
  m.fill(value);
  ASSERT_THAT(m.data(), Each(Eq(value)));
}

TEST_F(MatrixF, assignment) {
  iota();
  auto reference = std::vector<double>(m.size());
  std::iota(begin(reference), end(reference), 0.);
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
}

TEST_F(MatrixF, remove_row) {
  iota();
  auto ref_remove = std::array<std::vector<double>, 3>{
      {{4, 5, 6, 7, 8, 9, 10, 11}, {0, 1, 2, 3, 8, 9, 10, 11}, {0, 1, 2, 3, 4, 5, 6, 7}}};
  auto dimensions = Matrix<double>::coord_type{r - 1, c};
  for (size_t i = 0; i < 3; ++i) {
    auto m0 = m;
    m0.remove_row(i);
    ASSERT_EQ(m0.dimensions(), dimensions);
    ASSERT_THAT(m0.data(), Pointwise(Eq(), ref_remove[i])) << "removed row " << i;
  }
}

TEST_F(MatrixF, remove_col) {
  iota();
  const auto columns = std::vector<size_t>{0, 2, 3};
  //@formatter:off
  auto ref_remove = std::array<std::vector<double>, 3>{
      {{1, 2, 3,
        5, 6, 7,
        9, 10, 11},
       {0, 1, 3,
        4, 5, 7,
        8, 9, 11},
       {0, 1, 2,
        4, 5, 6,
        8, 9, 10}}};
  //@formatter:on
  auto dimensions = Matrix<double>::coord_type{r, c - 1};
  for (size_t i = 0; i < columns.size(); ++i) {
    auto m0 = m;
    m0.remove_col(columns[i]);
    ASSERT_EQ(m0.dimensions(), dimensions);
    ASSERT_THAT(m0.data(), Pointwise(Eq(), ref_remove[i])) << "removed column " << columns[i];
  }
}

TEST_F(MatrixF, remove_row_col) {
  iota();
  const auto row_col = std::vector<std::pair<size_t, size_t>>{{0, 0}, {1, 2}, {2, 3}};
  //@formatter:off
  auto ref_remove = std::array<std::vector<double>, 3>{
      {{5, 6, 7,
        9, 10,11},
       {0, 1, 3,
        8, 9, 11},
       {0, 1, 2,
        4, 5, 6,}}};
  //@formatter:on
  auto dimensions = Matrix<double>::coord_type{r-1, c - 1};
  for (size_t i = 0; i < row_col.size(); ++i) {
    auto m0 = m;
    m0.remove_row_col(row_col[i].first, row_col[i].second);
    ASSERT_EQ(m0.dimensions(), dimensions);
    ASSERT_THAT(m0.data(), Pointwise(Eq(), ref_remove[i]))
        << "removed row " << row_col[i].first << " and column " << row_col[i].second;
  }
}

TEST_F(MatrixF, slice_constructor) {
  ASSERT_NO_THROW((m.slice({0, 0}, {0, 0})));
  ASSERT_NO_THROW((m.slice({0, 0}, {r, c})));
  ASSERT_THROW((m.slice({0, 0}, {r + 1, c + 1})), std::runtime_error);
  ASSERT_THROW((m.slice({-1, -1}, {r, c})), std::runtime_error);
  ASSERT_THROW((m.slice({-1, 0}, {r, c + 1})), std::runtime_error);
}

TEST_F(MatrixF, slice_copy_empty) {
  auto m_right = m;
  m.fill(0);
  m_right.fill(1);
  const auto reference = m.data();
  m.slice({0, 0}, {0, 0}) = m_right.slice({0, 0}, {0, 0});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
  m.slice({0, 0}, {0, c}) = m_right.slice({0, 0}, {0, c});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
  m.slice({0, 0}, {r, 0}) = m_right.slice({0, 0}, {r, 0});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
  m.slice({r, 0}, {r, c}) = m_right.slice({r, 0}, {r, c});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
  m.slice({0, c}, {r, c}) = m_right.slice({0, c}, {r, c});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
  m.slice({r, c}, {r, c}) = m_right.slice({r, c}, {r, c});
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
}

TEST_F(MatrixF, slice_no_params) {
  auto m_right = m;
  const double value = 0.1;
  m_right.fill(value);
  ASSERT_THAT(m.data(), Pointwise(Ne(), m_right.data()));
  m.slice() = m_right;
  ASSERT_THAT(m.data(), Pointwise(Eq(), m_right.data()));
}

TEST_F(MatrixF, slice_copy_full_matrix) {
  auto m_right = m;
  const double value = 3.14;
  m.fill(0.);
  m_right.fill(value);
  ASSERT_THAT(m.data(), Pointwise(Ne(), m_right.data()));
  ASSERT_THROW((m.slice({0, 0}, {1, 1}) = m_right.slice({0, 0}, {2, 2})), std::runtime_error)
      << "incompatible dimensions";
  m.slice({0, 0}, m.dimensions()) = m_right.slice({0, 0}, m.dimensions());
  ASSERT_THAT(m.data(), Pointwise(Eq(), m_right.data()));
}

TEST_F(MatrixF, cslice_copy_full_matrix) {
  auto m_right = m;
  const double value = 3.14;
  m.fill(0.);
  m_right.fill(value);
  ASSERT_THAT(m.data(), Pointwise(Ne(), m_right.data()));
  const auto const_mat = m_right;
  m.slice({0, 0}, m.dimensions()) = const_mat.slice({0, 0}, m.dimensions());
  ASSERT_THAT(m.data(), Pointwise(Eq(), m_right.data()));
}

TEST_F(MatrixF, slice_copy_block) {
  auto block_row = r / 2;
  auto block_col = c / 2;
  auto m_right = Matrix<double>({block_row, block_col});
  const double value = 3.14;
  m.fill(value);
  auto reference = std::vector<double>(m.size(), value);
  for (size_t i = 0, ij = 0; i < block_row; ++i)
    for (size_t j = 0; j < block_col; ++j, ++ij) {
      m_right(i, j) = ij;
      reference[i * block_col + j] = ij;
    }
  m.slice({0, 0}, {block_row, block_col}) = m_right.slice({0, 0}, m_right.dimensions());
  ASSERT_THAT(m.data(), Pointwise(Eq(), reference));
}

TEST_F(MatrixF, slice_axpy) {
  auto m_right = m;
  const double alpha = 3.14;
  const double beta = 1.1;
  const double a = -0.5;
  m.fill(alpha);
  m_right.fill(beta);
  ASSERT_THAT(m.data(), Pointwise(Ne(), m_right.data()));
  m.slice().axpy(a, m_right.slice());
  ASSERT_THAT(m.data(), Each(Eq(alpha + a * beta)));
}