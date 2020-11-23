#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/wrap.h>
#include <molpro/linalg/itsolv/wrap_util.h>

#include <list>

using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::cwrap_arg;
using molpro::linalg::itsolv::decay;
using molpro::linalg::itsolv::decay_t;
using molpro::linalg::itsolv::find_ref;
using molpro::linalg::itsolv::remove_elements;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;
using molpro::linalg::itsolv::wrap_arg;
using ::testing::Eq;
using ::testing::Pointwise;

struct CompareEqualRefWrap {
  template <typename T, typename R>
  bool operator()(const std::reference_wrapper<T>& l, const std::reference_wrapper<R>& r) {
    return std::addressof(l.get()) == std::addressof(r.get());
  }
};

struct CompareEqualVecRef {
  template <typename T, typename R>
  bool operator()(const VecRef<T>& l, const VecRef<R>& r) {
    return std::equal(begin(l), end(l), begin(r), CompareEqualRefWrap{});
  }
};

TEST(wrap_util, find_ref_empty) {
  auto params = std::vector<int>{};
  auto result = find_ref(wrap(params), begin(params), end(params));
  ASSERT_TRUE(result.empty());
}

TEST(wrap_util, find_ref_all) {
  auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, find_ref_all_const) {
  const auto params = std::vector<int>{1, 2, 3, 4};
  auto indices = find_ref(wrap(params), begin(params), end(params));
  auto ref_indices = std::vector<size_t>{0, 1, 2, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, find_ref_some) {
  auto params = std::vector<int>{1, 2, 3, 4, 5};
  auto wparams = wrap(params);
  wparams.erase(wparams.begin());
  wparams.erase(wparams.begin() + 1);
  wparams.erase(wparams.begin() + 2);
  auto indices = find_ref(wparams, begin(params), end(params));
  auto ref_indices = std::vector<size_t>{1, 3};
  ASSERT_FALSE(indices.empty());
  ASSERT_THAT(indices, Pointwise(Eq(), ref_indices));
}

TEST(wrap_util, remove_elements) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  auto indices = std::vector<size_t>{0, 2, 3, 5};
  auto params_ref = std::vector<int>{2, 5};
  auto result = remove_elements(params, indices);
  ASSERT_EQ(result.size(), params_ref.size());
  ASSERT_THAT(result, Pointwise(Eq(), params_ref));
}

TEST(wrap, iterable_container__list) {
  auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap(params);
  ASSERT_FALSE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(wrap, iterable_container__const_list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, iterable_container__list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = cwrap(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, iterable_container__const_list) {
  const auto params = std::list<int>{1, 2, 3, 4, 5, 6};
  auto wparams = cwrap(params);
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, VecRef_to_CVecRef) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  auto wparams = wrap(params);
  auto cwparams = cwrap(wparams);
  ASSERT_TRUE(std::is_const_v<decltype(cwparams)::value_type::type>);
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), wparams));
}

TEST(cwrap_arg, one_arguments) {
  int a = 1;
  auto wa = cwrap_arg(a);
  ASSERT_TRUE(std::is_const_v<decltype(wa)::value_type::type>);
  ASSERT_EQ(wa.size(), 1);
  ASSERT_EQ(a, wa.front());
  ASSERT_EQ(std::addressof(a), std::addressof(wa.front().get()));
}

TEST(cwrap_arg, two_arguments) {
  auto ab_vec = std::vector<int>{1, 2};
  auto& a = ab_vec[0];
  auto& b = ab_vec[1];
  auto ab = cwrap_arg(a, b);
  ASSERT_TRUE(std::is_const_v<decltype(ab)::value_type::type>);
  ASSERT_THAT(ab, Pointwise(Eq(), ab_vec));
  ASSERT_TRUE(CompareEqualVecRef{}(wrap(ab_vec), ab));
}

TEST(cwrap_arg, four_arguments) {
  auto vec = std::vector<int>{1, 2, 3, 4};
  auto& a = vec[0];
  auto& b = vec[1];
  auto& c = vec[2];
  auto& d = vec[3];
  auto w = cwrap_arg(a, b, c, d);
  ASSERT_TRUE(std::is_const_v<decltype(w)::value_type::type>);
  ASSERT_THAT(w, Pointwise(Eq(), vec));
  ASSERT_TRUE(CompareEqualVecRef{}(wrap(vec), w));
}

TEST(wrap_arg, four_arguments) {
  auto vec = std::vector<int>{1, 2, 3, 4};
  auto& a = vec[0];
  auto& b = vec[1];
  auto& c = vec[2];
  auto& d = vec[3];
  auto w = wrap_arg(a, b, c, d);
  ASSERT_THAT(w, Pointwise(Eq(), vec));
  ASSERT_TRUE(CompareEqualVecRef{}(wrap(vec), w));
}

TEST(cwrap, const_VecRef_to_CVecRef) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  const auto wparams = wrap(params);
  auto cwparams = cwrap(wparams);
  ASSERT_TRUE(std::is_const_v<decltype(cwparams)::value_type::type>);
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), params));
  ASSERT_THAT(cwparams, Pointwise(Eq(), wparams));
  ASSERT_TRUE(CompareEqualVecRef{}(wparams, cwparams));
}

TEST(cwrap, forward_iterator) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  VecRef<int> wparams = wrap<int>(begin(params), end(params));
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, const_forward_iterator) {
  const auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  CVecRef<int> wparams = cwrap<int>(begin(params), end(params));
  ASSERT_TRUE(std::is_const_v<decltype(wparams)::value_type::type>);
  ASSERT_EQ(wparams.size(), params.size());
  ASSERT_THAT(wparams, Pointwise(Eq(), params));
}

TEST(cwrap, forward_iterator_on_VecRef) {
  auto params = std::vector<int>{1, 2, 3, 4, 5, 6};
  VecRef<int> wparams = wrap<int>(begin(params), end(params));
  VecRef<int> wparams2 = wrap<int>(begin(wparams), end(wparams));
  ASSERT_EQ(wparams.size(), wparams2.size());
  ASSERT_THAT(wparams2, Pointwise(Eq(), wparams));
  ASSERT_TRUE(CompareEqualVecRef{}(wparams2, wparams));
}
