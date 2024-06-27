#ifndef LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#define LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
#include <cmath>
#include <string>
#include <vector>
#include <optional>

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
    friend std::ostream& operator<<(std::ostream& os, const point& p);
  };
  /*!
   * @brief Construct the interpolant
   * @param p0 Defining point
   * @param p1 Defining point
   * @param interpolant The interpolation method. An exception is thrown if it is not one of the implemented values.
   * @param verbosity Values greater than zero show information on constructing and using the interpolant.
   */
  explicit Interpolate(point p0, point p1, std::optional<std::string> interpolant = "cubic", int verbosity = 0);
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
  Interpolate::point minimize(double xa, double xb, size_t bracket_grid = 100, size_t max_bracket_grid = 100000, bool analytic=true) const;
  Interpolate::point minimize_cubic() const;
  static std::vector<std::string> interpolants();
  const std::string& interpolant() const;

  const std::vector<double>& parameters() const { return m_parameters; }
  friend std::ostream& operator<<(std::ostream& os, const Interpolate& interpolant);

private:
  const point m_p0, m_p1;
  const std::string m_interpolant;
  double m_c0, m_c1, m_c2, m_c3; //< parameters defining the interpolant
  std::vector<double> m_parameters;
};

inline bool operator==(const Interpolate::point& lhs, const Interpolate::point& rhs) {
  return lhs.x == rhs.x && lhs.f == rhs.f && lhs.f1 == rhs.f1;
}
std::ostream& operator<<(std::ostream& os, const Interpolate& interpolant);
std::ostream& operator<<(std::ostream& os, const Interpolate::point& p);

} // namespace molpro::linalg::itsolv
#endif // LINEARALGEBRA_SRC_MOLPRO_LINALG_ITSOLV_INTERPOLATE_H_
