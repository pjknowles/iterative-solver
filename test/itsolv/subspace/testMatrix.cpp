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

  const size_t r = 3;
  const size_t c = 4;
  Matrix<double> m;
};

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