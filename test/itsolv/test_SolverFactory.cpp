#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/linalg/itsolv/SolverFactory.h>
#include "vector_types.h"

using molpro::linalg::itsolv::Options;
using molpro::linalg::itsolv::options_map;

TEST(SolverFactory, split_string) {
  using factory = molpro::linalg::itsolv::SolverFactory<Rvector, Qvector, Pvector>;
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
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>(
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
    auto solver = molpro::linalg::itsolv::create_LinearEquations<Rvector, Qvector, Pvector>(
        "convergence_threshold=1e-3,max_size_qspace=73, rubbish=trash");
    auto options = dynamic_cast<molpro::linalg::itsolv::LinearEquationsOptions*>(solver->get_options().get());
    EXPECT_TRUE(options->convergence_threshold.has_value());
    EXPECT_EQ(options->convergence_threshold.value(), 1e-3);
    EXPECT_TRUE(options->norm_thresh.has_value());
    EXPECT_NE(options->norm_thresh.value_or(777), 777);
    EXPECT_EQ(options->max_size_qspace.value(), 73);
  }

  for (const auto& method : std::vector<std::string>{"DIIS"}) {
    auto solver = molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector, Pvector>(
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
    auto solver = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector, Pvector>(
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
