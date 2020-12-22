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

void StringFacet::crop_space(std::string &path) {
  path.erase(path.begin(), std::find_if(path.begin(), path.end(), [](auto &el) { return !std::isspace(el); }));
  path.erase(std::find_if(path.rbegin(), path.rend(), [](auto &el) { return !std::isspace(el); }).base(), path.end());
}

} // namespace molpro::linalg::itsolv::util
