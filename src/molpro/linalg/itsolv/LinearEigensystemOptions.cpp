#include "LinearEigensystemOptions.h"
#include <locale>
namespace molpro::linalg::itsolv {
LinearEigensystemOptions::LinearEigensystemOptions(const std::map<std::string, std::string>& opt)
    : ILinearEigensystemOptions(opt) {
  auto& facet = std::use_facet<std::ctype<char>>(std::locale());
  auto upper = [&facet](std::string s) {
    facet.toupper(s.data(), s.data() + s.size());
    return s;
  };
  auto opt_upper = decltype(opt){};
  for (const auto& [key, val] : opt) {
    opt_upper[upper(key)] = val;
  }
  auto present = [&opt_upper](const auto& key) { return opt_upper.find(key) != opt_upper.end(); };
  if (auto key = upper("reset_D"); present(key)) {
    reset_D = std::stoi(opt_upper.at(key));
  };
  if (auto key = upper("reset_D_max_Q_size"); present(key)) {
    reset_D_max_Q_size = std::stod(opt_upper.at(key));
  }
  if (auto key = upper("max_size_qspace"); present(key)) {
    max_size_qspace = std::stod(opt_upper.at(key));
  }
  if (auto key = upper("norm_thresh"); present(key)) {
    norm_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = upper("svd_thresh"); present(key)) {
    svd_thresh = std::stod(opt_upper.at(key));
  }
  if (auto key = upper("hermiticity"); present(key)) {
    auto val = upper(opt_upper.at(key));
    if (val == "TRUE" || val == "T" || val == "1") {
      hermiticity = true;
    } else if (val == "TRUE" || val == "T" || val == "1") {
      hermiticity = false;
    } else {
      throw std::runtime_error("boolean value must be one of {TRUE, FALSE, T, F, 0, 1}");
    }
  }
}

} // namespace molpro::linalg::itsolv
