#ifndef LINEARALGEBRA_TEST_ARRAY_DATA_UTIL_H
#define LINEARALGEBRA_TEST_ARRAY_DATA_UTIL_H
namespace molpro {
namespace linalg {
namespace test {
// File contents: "/dataset" where dataset is a float64 arange(0,30)
const std::string name_single_dataset{std::string(ARRAY_DATA) + "/single_dataset.hdf5"};

// File contents: "/group1/group2/dataset" where dataset is a float64 arange(0,30)
const std::string name_inner_group_dataset{std::string(ARRAY_DATA) + "/inner_group_dataset.hdf5"};

// Test files for testDistrArrayHDF5
const std::string test_file_hdf5_n1{std::string(ARRAY_DATA) + "/test_distr_array_hdf5_n1.hdf5"};
const std::string test_file_hdf5_n2{std::string(ARRAY_DATA) + "/test_distr_array_hdf5_n2.hdf5"};
} // namespace test
} // namespace linalg
} // namespace molpro
#endif // LINEARALGEBRA_TEST_ARRAY_DATA_UTIL_H
