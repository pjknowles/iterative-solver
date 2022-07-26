#include "options.h"
std::shared_ptr<molpro::Options> s_options;

const std::shared_ptr<const molpro::Options> molpro::linalg::options() {
  if (s_options.get() == nullptr) {
    s_options.reset(new molpro::Options("ITERATIVE-SOLVER", ""));
  }
  return s_options;
}

void molpro::linalg::set_options(const molpro::Options& options) {
  s_options = std::make_shared<molpro::Options>(options);
}