#include "testDistrArray.h"

#include <molpro/linalg/array/DistrArrayMPI3.h>

using molpro::linalg::array::DistrArrayMPI3;

using ArrayTypes = ::testing::Types<DistrArrayMPI3>;
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, TestDistrArray, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistArrayInitializationF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistrArrayRangeF, ArrayTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(MPI3, DistrArrayCollectiveOpF, ArrayTypes);

struct DistributionF : public ::testing::Test, public DistrArrayMPI3 {
  DistributionF() : DistrArrayMPI3{m_dim, mpi_comm} {}

  DistributionMPI3 factory(int n_proc, size_t dim) { return DistributionMPI3(n_proc, dim); }
  static const size_t m_dim = 30;
};

TEST_F(DistributionF, constructor) {
  int np = 5;
  auto d = factory(np, 17);
  auto ref_distribution = std::vector<std::pair<unsigned long, size_t>>{{0, 4}, {4, 4}, {8, 3}, {11, 3}, {14, 3}};
  std::vector<decltype((d.range)(1))> proc_buffer;
  for (size_t i = 0; i < np; ++i)
    proc_buffer.emplace_back(d.range(i));
  ASSERT_THAT(proc_buffer, ContainerEq(ref_distribution));
}

TEST_F(DistributionF, process_map) {
  int np = 3;
  auto d = factory(np, 20);
  auto ref_lo_hi = std::vector<std::pair<unsigned long, unsigned long>>{{0, 0},   {1, 5},   {3, 7},   {5, 9},
                                                                        {13, 14}, {15, 19}, {18, 12}, {17, 6}};
  auto ref_process_pairs =
      std::vector<std::pair<int, int>>{{0, 0}, {0, 0}, {0, 1}, {0, 1}, {1, 2}, {2, 2}, {2, 1}, {2, 0}};
  auto process_pairs = std::vector<std::pair<int, int>>(ref_lo_hi.size());
  std::transform(ref_lo_hi.begin(), ref_lo_hi.end(), process_pairs.begin(),
                 [&d](auto p) { return d.locate_process(p.first, p.second); });
  ASSERT_THAT(process_pairs, ContainerEq(ref_process_pairs));
}
