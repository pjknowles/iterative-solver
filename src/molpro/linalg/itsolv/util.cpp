#include "util.h"
namespace molpro::linalg::itsolv::util {

std::string StringFacet::toupper(std::string in) {
  facet.toupper(in.data(), in.data() + in.size());
  return in;
}

std::string StringFacet::tolower(std::string in) {
  facet.tolower(in.data(), in.data() + in.size());
  return in;
}

} // namespace molpro::linalg::itsolv::util
