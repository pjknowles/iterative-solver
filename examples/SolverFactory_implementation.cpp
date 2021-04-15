#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
//#include <molpro/linalg/array/DistrArrayFile.h>
#include <vector>
#include <map>
template class molpro::linalg::itsolv::SolverFactory<std::vector<double>,std::vector<double>,std::map<size_t,double>>;
//template class molpro::linalg::itsolv::SolverFactory<std::vector<double>,molpro::linalg::array::DistrArrayFile,std::map<size_t,double>>;
