#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
#include <vector>
using R = std::vector<double>;
using Q = std::vector<double>;
using P = std::vector<double>;

template class molpro::linalg::itsolv::SolverFactory<R,Q,P>;
