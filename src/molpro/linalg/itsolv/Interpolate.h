#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#include <cmath>
#include <string>

/*!
 * @brief 4-parameter interpolation of a 1-dimensional function given two points for which function values and first
 * derivatives are known.
 */
namespace molpro::linalg::itsolv {

class Interpolate {
public:
  struct point {
    double x;
    double f;
    double f1;
    double f2 = std::nan("unset");
  };
  explicit Interpolate(point p0, point p1, std::string interpolant = "cubic");
  /*!
   * @brief Evaluate the interpolant and its derivative at a given point
   * @param x
   * @return
   */
  point operator()(double x) const;
  /*!
   * @brief Find the minimum of the interpolant within a range
   * @param xmin lower bound of range
   * @param xmax upper bound of range
   * @return
   */
  point minimize(double xmin, double xmax) const;

private:
  const point m_p0, m_p1;
  const std::string m_interpolant;
  double m_c0, m_c1, m_c2, m_c3; //< parameters defining the interpolant
};

inline bool operator==(const Interpolate::point& lhs, const Interpolate::point& rhs) {
  return lhs.x == rhs.x && lhs.f == rhs.f && lhs.f1 == rhs.f1;
}

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
