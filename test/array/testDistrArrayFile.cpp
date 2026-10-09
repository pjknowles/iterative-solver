#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <deque>
#include <vector>

#include "parallel_util.h"

#include <molpro/linalg/array/DistrArrayFile.h>
#include <molpro/linalg/array/DistrArraySpan.h>
#include <molpro/linalg/array/default_handler.h>
#include <molpro/linalg/array/util.h>
#include <molpro/linalg/array/util/BufferManager.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <molpro/linalg/array/util/gemm.h>
#include <molpro/linalg/options.h>

#ifdef LINEARALGEBRA_ARRAY_MPI3
#include <molpro/linalg/array/DistrArrayMPI3.h>
using molpro::linalg::array::DistrArrayMPI3;
#endif

using molpro::linalg::array::ArrayHandlerDDiskDistr;
using molpro::linalg::array::default_handler;
using molpro::linalg::array::DistrArrayFile;
using molpro::linalg::array::DistrArraySpan;
using molpro::linalg::array::Span;
using molpro::linalg::array::util::LockMPI3;
using molpro::linalg::array::util::make_distribution_spread_remainder;
using molpro::linalg::array::util::ScopeLock;
using molpro::linalg::test::mpi_comm;

using ::testing::ContainerEq;
using ::testing::DoubleEq;
using ::testing::Each;
using ::testing::Pointwise;

class DistrArrayFile_Fixture : public ::testing::Test {
public:
  DistrArrayFile_Fixture() : a(DistrArrayFile(size, mpi_comm)) {}
  void SetUp() override {
    auto dist = a.distribution();
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    MPI_Comm_size(mpi_comm, &mpi_size);
    left = dist.range(mpi_rank).first;
    right = dist.range(mpi_rank).second;
    for (int i = 0; i < mpi_size; i++) {
      displs.push_back(dist.range(i).first);
      chunks.push_back(dist.range(i).second - dist.range(i).first);
    }
  };
  void TearDown() override{};
  const size_t size = 1200;
  int mpi_size, mpi_rank;
  int left, right;
  std::vector<int> chunks, displs;
  DistrArrayFile a;
};

TEST(DistrArrayFile, constructor_size) {
  {
    auto a = DistrArrayFile(100, mpi_comm);
    LockMPI3 lock{mpi_comm};
    {
      auto l = lock.scope();
      EXPECT_EQ(a.size(), 100);
    }
  }
}

TEST_F(DistrArrayFile_Fixture, constructor_copy) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<double> w(size, 0), x(size, 0);
  auto b = a;
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
  //  b = a; // TODO operator=() not working yet
  //  b.get(dist.range(mpi_rank).first, dist.range(mpi_rank).second, v.data());
  //  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

#ifdef LINEARALGEBRA_ARRAY_MPI3
TEST(DistrArrayFile, constructor_copy_from_distr_array) {
  const double val = 0.5;
  auto a_mem = molpro::linalg::array::DistrArrayMPI3(100, mpi_comm);
  a_mem.fill(val);
  auto a_disk = DistrArrayFile{a_mem};
  LockMPI3 lock{mpi_comm};
  auto vec = a_disk.vec();
  EXPECT_THAT(vec, Each(DoubleEq(val)));
  {
    auto l = lock.scope();
    EXPECT_EQ(a_disk.communicator(), a_mem.communicator());
    EXPECT_EQ(a_disk.size(), a_mem.size());
    EXPECT_TRUE(a_disk.distribution().compatible(a_mem.distribution()));
  }
}
#endif

TEST_F(DistrArrayFile_Fixture, handler_copies) {
  using molpro::mpi::comm_global;

  auto handler = ArrayHandlerDDiskDistr<DistrArrayFile, DistrArraySpan>{};
  size_t n = 3, dim = size;
  std::vector<std::vector<double>> vx(n, std::vector<double>(dim)), vy(n, std::vector<double>(dim));
  std::vector<DistrArraySpan> mem_vecs;
  std::vector<DistrArrayFile> file_vecs;
  mem_vecs.reserve(n);
  file_vecs.reserve(n);
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  for (size_t i = 0; i < n; i++) {
    std::iota(vx[i].begin(), vx[i].end(), i + 0.5);
    auto crange = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
    auto clength = crange.second - crange.first;
    mem_vecs.emplace_back(dim, Span<DistrArraySpan::value_type>(&vx[i][crange.first], clength), comm_global());
    file_vecs.emplace_back(handler.copy(mem_vecs.back()));
  }
  for (size_t i = 0; i < n; i++) {
    file_vecs[i].get(left, right, &(*(vy[i].begin() + left)));
    // vy[i] = file_vecs[i].vec();
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, vy[i].data(), chunks.data(), displs.data(), MPI_DOUBLE,
                   mpi_comm);
    EXPECT_THAT(vx[i], Pointwise(DoubleEq(), vy[i]));
  }
}

TEST_F(DistrArrayFile_Fixture, constructor_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  DistrArrayFile b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, assignment_move) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  auto b = DistrArrayFile(1);
  b = std::move(a);
  std::vector<double> w(size, 0);
  b.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_EQ(b.size(), size);
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST(DistrArrayFile, compatible) {
  auto a = DistrArrayFile{100, mpi_comm};
  auto b = DistrArrayFile{1000, mpi_comm};
  ScopeLock l{mpi_comm};
  EXPECT_TRUE(a.compatible(a));
  EXPECT_TRUE(b.compatible(b));
  EXPECT_FALSE(a.compatible(b));
  EXPECT_EQ(a.compatible(b), b.compatible(a));
}

TEST_F(DistrArrayFile_Fixture, writeread) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0.5);
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<double> w(size, 0);
  a.get(left, right, &(*(w.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(v, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, accumulate) {
  std::vector<double> v(size), w(size), x(size), y(size);
  std::iota(v.begin(), v.end(), 0.5);
  std::iota(w.begin(), w.end(), 0.5);
  std::transform(v.begin(), v.end(), w.begin(), x.begin(), [](auto& l, auto& r) { return l + r; });
  a.put(left, right, &(*(v.cbegin() + left)));
  a.acc(left, right, &(*(w.begin() + left)));
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), x));
}

TEST_F(DistrArrayFile_Fixture, gather) {
  std::vector<double> v(size), w(size);
  int n = -2;
  std::generate(v.begin(), v.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(chunks[mpi_rank]);
  std::iota(x.begin(), x.end(), left);
  auto tmp = a.gather(x);
  for (size_t i = 0; i < x.size(); i++) {
    w[left + i] = tmp[i];
  }
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, w.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(w, Pointwise(DoubleEq(), v));
}

TEST_F(DistrArrayFile_Fixture, scatter) {
  std::vector<double> v(size, 0), w(size), y(size);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(chunks[mpi_rank]);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(chunks[mpi_rank]);
  for (size_t i = 0; i < x.size(); i++) {
    tmp[i] = w[i + left];
  }
  a.scatter(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, scatter_acc) {
  std::vector<double> v(size), w(size), y(size);
  std::iota(v.begin(), v.end(), 0);
  int n = -2;
  std::generate(w.begin(), w.end(), [&n] { return n += 2; });
  a.put(left, right, &(*(v.cbegin() + left)));
  std::vector<DistrArrayFile::index_type> x(chunks[mpi_rank]);
  std::iota(x.begin(), x.end(), left);
  std::vector<double> tmp(chunks[mpi_rank]);
  for (size_t i = 0; i < x.size(); i++) {
    tmp[i] = v[i + left];
  }
  a.scatter_acc(x, tmp);
  a.get(left, right, &(*(y.begin() + left)));
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, y.data(), chunks.data(), displs.data(), MPI_DOUBLE, mpi_comm);
  ScopeLock l{mpi_comm};
  EXPECT_THAT(y, Pointwise(DoubleEq(), w));
}

TEST_F(DistrArrayFile_Fixture, dot_DistrArray) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0);
  a.put(left, right, &(*(v.cbegin() + left)));
  const DistrArraySpan s(size, Span<double>(&(*(v.begin() + left)), right - left));
  auto ss = s.dot(s);
  auto as = a.dot(s);
  //  auto aa = a.dot(a);
  auto sa = s.dot(a);
  ScopeLock l{mpi_comm};
  EXPECT_NEAR(ss, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
  EXPECT_NEAR(as, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
  //  EXPECT_NEAR(aa,size*(size-1)*(2*size-1)/6,1e-13);
  EXPECT_NEAR(sa, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
}

TEST_F(DistrArrayFile_Fixture, dot_DistrArrayFile) {
  std::vector<double> v(size);
  std::iota(v.begin(), v.end(), 0);
  a.put(left, right, &(*(v.cbegin() + left)));
  const DistrArraySpan s(size, Span<double>(&(*(v.begin() + left)), right - left));
  const DistrArrayFile f(s);
  auto ff = f.dot(f);
  auto ff_base = f.dot(static_cast<const molpro::linalg::array::DistrArray&>(f));
  auto as = a.dot(f);
  auto sa = f.dot(a);
  ScopeLock l{mpi_comm};
  EXPECT_NEAR(ff, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
  EXPECT_NEAR(ff_base, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
  EXPECT_NEAR(as, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
  EXPECT_NEAR(sa, size * (size - 1) * (2 * size - 1) / 6, 1e-13);
}

TEST_F(DistrArrayFile_Fixture, contiguous_allocation) {

  size_t n = 10;
  size_t dim = 10;

  auto [cx, cy, cz] = molpro::linalg::test::get_contiguous(n, dim);

  auto cx_wrapped = molpro::linalg::itsolv::cwrap(cx);

  // check contiguity
  int previous_stride = 0;
  for (size_t j = 0; j < n - 1; ++j) {
    auto unique_ptr_j = cx_wrapped.at(j).get().local_buffer()->data();
    auto unique_ptr_jp1 = cx_wrapped.at(j + 1).get().local_buffer()->data();
    int stride = unique_ptr_jp1 - unique_ptr_j;
    if (j > 0) {
      if (stride != previous_stride) {
        throw std::runtime_error("yy doesn't have a consistent stride\n");
      }
    }
    previous_stride = stride;
  }
}

// A DistrArrayFile has no local buffer (issue #121): elementwise operations, also those mixing it with arrays in
// memory, page through it. With a page much smaller than the local section, every operation must give the same result
// as on DistrArraySpan.
TEST(DistrArrayFile, paged_operations_match_memory) {
  using molpro::linalg::array::DistrArray;
  using molpro::mpi::comm_global;
  const size_t dim = 23, page = 4;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm_global(), &mpi_rank);
  MPI_Comm_size(comm_global(), &mpi_size);
  const auto [lo, hi] = make_distribution_spread_remainder<size_t>(dim, mpi_size).range(mpi_rank);
  auto value = [](int series, size_t i) { return std::cos(0.7 * series + 1.3 * i) + 0.1 * series; };
  // the storage of the in-memory arrays, which must not move
  std::deque<std::vector<double>> storage;
  auto memory = [&](int series) {
    storage.emplace_back(dim);
    for (size_t i = 0; i < dim; ++i)
      storage.back()[i] = value(series, i);
    return DistrArraySpan(dim, Span<double>(&storage.back()[lo], hi - lo), comm_global());
  };
  auto disk = [&](int series) {
    DistrArrayFile a(dim);
    std::vector<double> v(hi - lo);
    for (size_t i = lo; i < hi; ++i)
      v[i - lo] = value(series, i);
    a.put(lo, hi, v.data());
    a.set_buffer_size(page);
    return a;
  };
  auto local = [&](const DistrArray& a) { return a.get(lo, hi); };
  auto expect_same = [&](const DistrArray& f, const DistrArray& m, const std::string& what) {
    EXPECT_THAT(local(f), Pointwise(DoubleEq(), local(m))) << what;
  };

  {
    auto f = disk(1);
    auto m = memory(1);
    f.fill(2.5), m.fill(2.5);
    expect_same(f, m, "fill");
    f.scal(-1.5), m.scal(-1.5);
    f.add(0.25), m.add(0.25);
    f.recip(), m.recip();
    expect_same(f, m, "scal, add, recip");
  }
  {
    auto f = disk(1), f2 = disk(2);
    auto m = memory(1), m2 = memory(2), m3 = memory(3), m3copy = memory(3);
    f.axpy(0.3, m2), m.axpy(0.3, m2);
    expect_same(f, m, "disk += memory");
    m3.axpy(-0.7, f2), m3copy.axpy(-0.7, m2);
    expect_same(m3, m3copy, "memory += disk");
    f.axpy(1.1, f2), m.axpy(1.1, m2);
    expect_same(f, m, "disk += disk");
    f.axpy(2, f), m.axpy(2, m);
    expect_same(f, m, "disk += disk (itself)");
  }
  {
    auto f = disk(1), f2 = disk(2);
    auto m = memory(1), m2 = memory(2);
    EXPECT_NEAR(f.dot(m2), m.dot(m2), 1e-13);
    EXPECT_NEAR(m2.dot(f), m2.dot(m), 1e-13);
    EXPECT_NEAR(f.dot(f2), m.dot(m2), 1e-13);
    EXPECT_NEAR(f.dot(f), m.dot(m), 1e-13);
  }
  {
    auto f = disk(1), f2 = disk(2);
    auto m = memory(1), m2 = memory(2), m4 = memory(4);
    f.copy(m4), m.copy(m4);
    expect_same(f, m, "copy disk <- memory");
    m2.copy(f2);
    expect_same(m2, f2, "copy memory <- disk");
    f.copy(f2), m.copy(m2);
    expect_same(f, m, "copy disk <- disk");
    auto g = disk(5);
    auto n = memory(5);
    g.copy_patch(m4, 3, 17), n.copy_patch(m4, 3, 17);
    expect_same(g, n, "copy_patch");
  }
  {
    auto f = disk(1), f2 = disk(2), f3 = disk(3);
    auto m = memory(1), m2 = memory(2), m3 = memory(3);
    f.times(m2), m.times(m2);
    expect_same(f, m, "times");
    f.times(f2, m3), m.times(m2, m3);
    expect_same(f, m, "times(y, z)");
    f.divide(m2, f3, 3.0, true, true), m.divide(m2, m3, 3.0, true, true);
    expect_same(f, m, "divide, appending");
    f.divide(f2, m3, 3.0), m.divide(m2, m3, 3.0);
    expect_same(f, m, "divide");
  }
  {
    auto f = disk(6);
    auto m = memory(6), m2 = memory(2);
    for (bool max : {false, true})
      for (bool ignore_sign : {false, true})
        EXPECT_EQ(f.select(5, max, ignore_sign), m.select(5, max, ignore_sign)) << max << ignore_sign;
    EXPECT_EQ(f.min_n(4), m.min_n(4));
    EXPECT_EQ(f.max_n(4), m.max_n(4));
    EXPECT_EQ(f.min_abs_n(4), m.min_abs_n(4));
    EXPECT_EQ(f.max_abs_n(4), m.max_abs_n(4));
    EXPECT_EQ(f.select_max_dot(5, m2), m.select_max_dot(5, m2));
    const DistrArray::SparseArray sparse{{0, 1.5}, {3, -2.0}, {4, 0.5}, {11, 3.0}, {12, -1.0}, {22, 0.75}};
    EXPECT_EQ(f.select_max_dot(3, sparse), m.select_max_dot(3, sparse));
    // the largest |f_i s_i| over the whole array, found directly
    auto products = std::vector<std::pair<double, size_t>>();
    for (const auto& [i, v] : sparse)
      products.emplace_back(std::abs(value(6, i) * v), i);
    std::sort(products.rbegin(), products.rend());
    auto expected = std::map<size_t, double>();
    for (size_t k = 0; k < 3; ++k)
      expected.emplace(products[k].second, products[k].first);
    const auto selected = f.select_max_dot(3, sparse);
    ASSERT_EQ(selected.size(), expected.size());
    auto ie = expected.cbegin();
    for (const auto& [index, product] : selected) {
      EXPECT_EQ(index, ie->first);
      EXPECT_NEAR(product, ie->second, 1e-14);
      ++ie;
    }
    EXPECT_NEAR(f.dot(sparse), m.dot(sparse), 1e-13);
    f.axpy(-0.5, sparse), m.axpy(-0.5, sparse);
    expect_same(f, m, "axpy sparse");
  }
  {
    // gemm with arrays on disk on both sides, and with sparse arrays, paging gemm as well
    using molpro::linalg::array::util::gemm_inner_distr_distr;
    using molpro::linalg::array::util::gemm_inner_distr_sparse;
    using molpro::linalg::array::util::gemm_outer_distr_distr;
    using molpro::linalg::array::util::gemm_outer_distr_sparse;
    using molpro::linalg::itsolv::cwrap;
    using molpro::linalg::itsolv::wrap;
    molpro::linalg::set_options(molpro::Options("ITERATIVE-SOLVER", "GEMM_PAGESIZE=4"));
    std::vector<DistrArrayFile> fx, fy;
    std::vector<DistrArraySpan> mx, my;
    for (int k = 0; k < 3; ++k) {
      fx.push_back(disk(10 + k)), fy.push_back(disk(20 + k));
      mx.push_back(memory(10 + k)), my.push_back(memory(20 + k));
    }
    const auto inner_disk = gemm_inner_distr_distr(cwrap(fy), cwrap(fx));
    const auto inner_memory = gemm_inner_distr_distr(cwrap(my), cwrap(mx));
    for (size_t k = 0; k < inner_memory.size(); ++k)
      EXPECT_NEAR(inner_disk.data()[k], inner_memory.data()[k], 1e-13) << "gemm_inner disk x disk " << k;
    auto alpha = molpro::linalg::itsolv::subspace::Matrix<double>({3, 3});
    alpha.fill(0.5);
    gemm_outer_distr_distr(alpha, cwrap(fx), wrap(fy));
    gemm_outer_distr_distr(alpha, cwrap(mx), wrap(my));
    for (int k = 0; k < 3; ++k)
      for (size_t i = 0; i < hi - lo; ++i)
        EXPECT_NEAR(local(fy[k])[i], local(my[k])[i], 1e-13) << "gemm_outer disk x disk " << k << ", " << i;
    const std::vector<DistrArray::SparseArray> sp{{{1, 0.5}, {9, -1.0}, {21, 2.0}}, {{0, 1.0}, {13, 0.25}}};
    const auto sinner_disk = gemm_inner_distr_sparse(cwrap(fx), cwrap(sp));
    const auto sinner_memory = gemm_inner_distr_sparse(cwrap(mx), cwrap(sp));
    for (size_t k = 0; k < sinner_memory.size(); ++k)
      EXPECT_NEAR(sinner_disk.data()[k], sinner_memory.data()[k], 1e-13) << "gemm_inner sparse " << k;
    auto beta = molpro::linalg::itsolv::subspace::Matrix<double>({2, 3});
    beta.fill(-0.25);
    gemm_outer_distr_sparse(beta, cwrap(sp), wrap(fx));
    gemm_outer_distr_sparse(beta, cwrap(sp), wrap(mx));
    for (int k = 0; k < 3; ++k)
      for (size_t i = 0; i < hi - lo; ++i)
        EXPECT_NEAR(local(fx[k])[i], local(mx[k])[i], 1e-13) << "gemm_outer sparse " << k << ", " << i;
    molpro::linalg::set_options(molpro::Options("ITERATIVE-SOLVER", ""));
  }
}
