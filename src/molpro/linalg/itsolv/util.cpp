#include "util.h"
namespace molpro::linalg::itsolv::util {

std::string StringFacet::toupper(std::string in) const {
  facet.toupper(in.data(), in.data() + in.size());
  return in;
}

std::string StringFacet::tolower(std::string in) const {
  facet.tolower(in.data(), in.data() + in.size());
  return in;
}

bool StringFacet::tobool(const std::string &in) const {
  auto val = toupper(in);
  crop_space(val);
  bool result;
  if (val == "TRUE" || val == "T" || val == "1") {
    result = true;
  } else if (val == "FALSE" || val == "F" || val == "0") {
    result = false;
  } else {
    throw std::runtime_error("value =" + val + ", must be one of {TRUE, T, 1, FALSE, F, 0}");
  }
  return result;
}

void StringFacet::crop_space(std::string &path) {
  path.erase(path.begin(), std::find_if(path.begin(), path.end(), [](auto &el) { return !std::isspace(el); }));
  path.erase(std::find_if(path.rbegin(), path.rend(), [](auto &el) { return !std::isspace(el); }).base(), path.end());
}

static std::string trim_string(std::string s) {
  StringFacet::crop_space(s);
  return s;
}

std::map<std::string, std::string> StringFacet::parse_keyval_string(std::string s) {
  std::map<std::string, std::string> result;
  const std::string field_separators = ",;";
  const std::string keyval_separators = "=:";
  s += field_separators;
  crop_space(s);
  while (!s.empty()) {
    auto end = s.find_first_of(field_separators);
    auto entry = trim_string(s.substr(0, end));
    if (entry.empty())
      break;
    auto equal = entry.find_first_of(keyval_separators);
    if (equal == std::string::npos)
      throw std::runtime_error("String " + entry + " cannot be parsed as key" + field_separators.front() + "value");
    result[trim_string(entry.substr(0, equal))] = trim_string(entry.substr(equal + 1));
    s = trim_string(s.substr(end + 1));
  }
  return result;
}

} // namespace molpro::linalg::itsolv::util
