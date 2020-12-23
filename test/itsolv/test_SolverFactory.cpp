#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/SolverFactory.h>

using molpro::linalg::itsolv::Options;
using molpro::linalg::itsolv::options_map;
using R = std::vector<double>;
using Q = std::vector<double>;
using P = std::vector<double>;

TEST(SolverFactory, split_string) {
  using factory = molpro::linalg::itsolv::SolverFactory<R, Q, P>;
  std::string good{" key1=value1 , key2=value2,key3=  , key4=value4"};
  std::string bad{"keywithoutvalue"};
  auto sgood = factory::split_string(good);
  EXPECT_EQ(sgood["key1"], "value1");
  EXPECT_EQ(sgood["key2"], "value2");
  EXPECT_EQ(sgood["key3"], "");
  EXPECT_EQ(sgood["key4"], "value4");
  EXPECT_THROW(auto sbad = factory::split_string(bad), std::runtime_error);
}
TEST(SolverFactory, string_constructor) {
  {
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<R, Q, P>(
        "convergence_threshold=1e-3,max_size_qspace=73, n_roots=4");
    auto options = dynamic_cast<molpro::linalg::itsolv::Options*>(solver->get_options().get());
    EXPECT_TRUE(options->convergence_threshold.has_value());
    EXPECT_EQ(options->convergence_threshold.value(), 1e-3);
    {
      auto options = dynamic_cast<molpro::linalg::itsolv::LinearEigensystemOptions*>(solver->get_options().get());
      EXPECT_TRUE(options->norm_thresh.has_value());
      EXPECT_NE(options->norm_thresh.value_or(777), 777);
      EXPECT_EQ(options->max_size_qspace.value(), 73);
      EXPECT_EQ(options->n_roots.value(), 4);
    }
  }

  {
    auto solver = molpro::linalg::itsolv::create_LinearEquations<R, Q, P>(
        "convergence_threshold=1e-3,max_size_qspace=73, rubbish=trash");
    auto options = dynamic_cast<molpro::linalg::itsolv::LinearEquationsOptions*>(solver->get_options().get());
    EXPECT_TRUE(options->convergence_threshold.has_value());
    EXPECT_EQ(options->convergence_threshold.value(), 1e-3);
    EXPECT_TRUE(options->norm_thresh.has_value());
    EXPECT_NE(options->norm_thresh.value_or(777), 777);
    EXPECT_EQ(options->max_size_qspace.value(), 73);
  }

  for (const auto& method : std::vector<std::string>{"DIIS"}) {
    auto solver = molpro::linalg::itsolv::create_NonLinearEquations<R, Q, P>(
        "DIIS", "convergence_threshold=1e-3,max_size_qspace=73, rubbish=trash");
    auto options = dynamic_cast<molpro::linalg::itsolv::Options*>(solver->get_options().get());
    EXPECT_TRUE(options->convergence_threshold.has_value());
    EXPECT_EQ(options->convergence_threshold.value(), 1e-3);
    if (auto options =
            dynamic_cast<molpro::linalg::itsolv::NonLinearEquationsDIISOptions*>(solver->get_options().get())) {
      EXPECT_TRUE(options->norm_thresh.has_value());
      EXPECT_NE(options->norm_thresh.value_or(777), 777);
      EXPECT_EQ(options->max_size_qspace.value(), 73);
    }
  }

  for (const auto& method : std::vector<std::string>{"BFGS", "SD"}) {
    auto solver = molpro::linalg::itsolv::create_Optimize<R, Q, P>(
        method, std::string{"convergence_threshold=1e-3"} + (method == "BFGS" ? ",max_size_qspace=73" : ""));
    auto options = dynamic_cast<molpro::linalg::itsolv::Options*>(solver->get_options().get());
    EXPECT_TRUE(options->convergence_threshold.has_value());
    EXPECT_EQ(options->convergence_threshold.value(), 1e-3);
    if (method == "BFGS")
      if (auto options = dynamic_cast<molpro::linalg::itsolv::OptimizeBFGSOptions*>(solver->get_options().get())) {
        EXPECT_EQ(options->max_size_qspace.value(), 73);
      }
  }
}
