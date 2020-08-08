#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "parallel_util.h"
#include <molpro/linalg/array/util/DistrFlags.h>

using molpro::linalg::array::util::DistrFlags;
using molpro::linalg::test::mpi_comm;

TEST(DistrFlags, constructor_comm) { DistrFlags lock{mpi_comm}; }
TEST(DistrFlags, constructor_default) { DistrFlags lock{}; }
TEST(DistrFlags, constructor_copy) { DistrFlags lock{}; }
TEST(DistrFlags, constructor_move) { DistrFlags lock{}; }
TEST(DistrFlags, copy_assignment) { DistrFlags lock{}; }
TEST(DistrFlags, move_assignment) { DistrFlags lock{}; }


