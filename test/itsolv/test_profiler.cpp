#include <gtest/gtest.h>
#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/linalg/options.h>

#include <filesystem>
#include <molpro/linalg/itsolv/SolverFactory.h>

TEST(Profiler, default) {
  auto single = molpro::Profiler::single();
  const int initial_depth = 2;
  const std::string dotgraph_file{"test_profiler.gv"};
  single->set_max_depth(initial_depth);
  {
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>();
    EXPECT_EQ(single->get_max_depth(), 0);
  }
  EXPECT_EQ(single->get_max_depth(), initial_depth);
  molpro::linalg::set_options(molpro::Options("LINEARALGEBRA", std::string{"PROFILER_DEPTH=3, PROFILER_DOTGRAPH="} +
                                                                   dotgraph_file + ", PROFILER_THRESHOLD=0"));
  std::filesystem::remove(dotgraph_file);

  EXPECT_EQ(single->get_max_depth(), initial_depth);
  {
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<std::vector<double>>();
    EXPECT_EQ(single->get_max_depth(), 3);
    solver->profiler()->set_max_depth(27);
    EXPECT_EQ(single->get_max_depth(), 27);
    solver->profiler()->push("Some operation or other");
//    std::cout << *solver->profiler() << std::endl;
  }
//  std::cout << *single << std::endl;
  EXPECT_EQ(single->get_max_depth(), initial_depth);
  EXPECT_TRUE(std::filesystem::exists(dotgraph_file));
//  auto dotgraph = single->dotgraph("", 0);
//  std::ifstream t(dotgraph_file);
//  std::stringstream buffer;
//  buffer << t.rdbuf();
//  EXPECT_EQ(buffer.str(),dotgraph);
}