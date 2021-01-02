#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#include <cmath>
#include <string>
#include <vector>

/*!
 * @brief 4-parameter interpolation of a 1-dimensional function given two points for which function values and first
 * derivatives are known.
 */
namespace molpro::linalg::itsolv {

class Interpolate {
public:
  struct point {
    double x;                      //< abscissa
    double f = std::nan("unset");  //< function value at x
    double f1 = std::nan("unset"); //< function first gradient at x
    double f2 = std::nan("unset"); //< function second gradient at x
  };
  /*!
   * @brief Construct the interpolant
   * @param p0 Defining point
   * @param p1 Defining point
   * @param interpolant The interpolation method. An exception is thrown if it is not one of the implemented values.
   */
  explicit Interpolate(point p0, point p1, std::string interpolant = "cubic");
  /*!
   * @brief Evaluate the interpolant and its derivative at a given point
   * @param x
   * @return
   */
  point operator()(double x) const;
  /*!
   * @brief Find the minimum of the interpolant within a range
   * @param xa first bound of range
   * @param xb second bound of range
   * @param bracket_grid number of intervals in |xa-xb|, to be considered in initial bracketing of the minimum. Large
   * values result in many function evaluations, but if set too small in cases of multiple minima, the global minimum
   * may not be found.
   * @return The minimum point. The result may be one of Interpolate(xa), Interpolate(xb) with non-zero first
   * derivative, if no other minimum was found in the interval
   */
  Interpolate::point minimize(double xa, double xb, size_t bracket_grid = 100, size_t max_bracket_grid = 100000) const;
  static std::vector<std::string> interpolants();

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
