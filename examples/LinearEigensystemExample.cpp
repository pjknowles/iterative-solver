#include <iomanip>
#include <molpro/linalg/IterativeSolver.h>
#include <molpro/linalg/array/ArrayHandlerSparse.h>

using molpro::linalg::array::ArrayHandler;
using molpro::linalg::array::ArrayHandlerIterable;
using molpro::linalg::array::ArrayHandlerIterableSparse;
using molpro::linalg::array::ArrayHandlerSparse;
using molpro::linalg::itsolv::ArrayHandlers;

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  for (const auto& s : v)
    o << " " << s;
  return o;
}

// Find lowest eigensolutions of M(i,j) = alpha*(i+1)*delta(i,j) + i + j
// Storage of vectors in-memory via class pv
using scalar = double;
size_t n;                     // dimension of problem
constexpr scalar alpha = 100; // separation of diagonal elements
std::vector<double> hmat;

class pv {
public:
  std::vector<scalar> buffer;
  using scalar_type = scalar;
  using value_type = scalar;

  pv() = default;
  explicit pv(size_t length, int option = 0) : buffer(length) {}

  // explicit pv(const pv& source, int option = 0) { *this = source; }
  pv(const pv& source, int option = 0) { *this = source; }

  std::vector<scalar>::const_iterator begin() const { return buffer.begin(); }

  std::vector<scalar>::iterator begin() { return buffer.begin(); }

  std::vector<scalar>::const_iterator end() const { return buffer.end(); }

  std::vector<scalar>::iterator end() { return buffer.end(); }

  size_t size() const { return buffer.size(); }

  void push_back(const scalar& elem) { buffer.push_back(elem); }

  void axpy(scalar a, const pv& other) {
    for (size_t i = 0; i < buffer.size(); i++)
      buffer[i] += a * other.buffer[i];
  }

  void axpy(scalar a, const std::map<size_t, scalar>& other) {
    for (const auto& o : other)
      (*this)[o.first] += a * o.second;
  }

  void scal(scalar a) {
    for (auto& b : buffer)
      b *= a;
  }

  void zero() {
    for (auto& b : buffer)
      b = 0;
  }

  scalar dot(const pv& other) const {
    scalar result = 0;
    for (size_t i = 0; i < buffer.size(); i++)
      result += buffer[i] * other.buffer[i];
    return result;
  }

  scalar dot(const std::map<size_t, scalar>& other) const {
    scalar result = 0;
    for (const auto& o : other)
      result += o.second * (*this)[o.first];
    return result;
  }

  std::tuple<std::vector<size_t>, std::vector<scalar>> select(const pv& measure, const size_t maximumNumber = 1000,
                                                              const scalar threshold = 0) const {
    return std::make_tuple(std::vector<size_t>(0), std::vector<scalar>(0)); // null implementation
  };

  const scalar& operator[](const size_t i) const { return buffer[i]; }

  scalar& operator[](size_t i) { return buffer[i]; }
};

scalar matrix(const size_t i, const size_t j) { return (i == j ? alpha * (i + 1) : 0) + i + j; }

void action(const std::vector<pv>& psx, std::vector<pv>& outputs) {
  for (size_t k = 0; k < psx.size(); k++) {
    for (size_t i = 0; i < n; i++) {
      outputs[k][i] = 0;
      for (size_t j = 0; j < n; j++)
        outputs[k][i] += matrix(i, j) * psx[k][j];
    }
  }
}

void update(std::vector<pv>& psc, const std::vector<pv>& psg, std::vector<scalar> shift = std::vector<scalar>()) {
  //  for (size_t k = 0; k < psc.size(); k++) {
  //   std::cout << "update root "<<k<<", shift="<<shift[k]<<": ";
  //    for (size_t i = 0; i < n; i++)
  //      std::cout <<" "<< -psg[k][i] / (1e-12-shift[k] + matrix(i,i));
  //    std::cout<<std::endl;
  //  }
  for (size_t k = 0; k < psc.size() and k < shift.size(); k++)
    for (size_t i = 0; i < n; i++)
      psc[k][i] -= psg[k][i] / (1e-12 - shift[k] + matrix(i, i));
}

int main(int argc, char* argv[]) {
  for (const auto& nn : std::vector<int>{1, 2, 4, 10, 100}) {
    n = nn;
    for (const auto& nroot : std::vector<int>{1, 2, 4}) {
      if (nroot > n)
        break;
      std::cout << "\nMatrix dimension=" << n << ", looking for " << nroot << " roots" << std::endl;
      //      std::cout << "diagonal elements";
      //      for (auto i = 0; i < n; i++)
      //        std::cout << " " << matrix(i, i);
      //      for (auto i = 0; i < n; i++) {
      //        for (auto j = 0; j < n; j++)
      //          std::cout << " " << matrix(i, j);
      //        std::cout << std::endl;
      //      }
      std::vector<double> diagonals(n);
      for (auto i = 0; i < n; i++)
        diagonals.push_back(matrix(i, i));
      std::cout << std::endl;

      auto handlers = std::make_shared<ArrayHandlers<pv, pv, std::map<size_t, double>>>();
      molpro::linalg::LinearEigensystem<pv> solver{handlers};
      solver.m_verbosity = 1;
      solver.m_roots = nroot;
      solver.m_thresh = 1e-9;
      std::vector<pv> g;
      std::vector<pv> x;
      for (size_t root = 0; root < solver.m_roots; root++) {
        x.emplace_back(n);
        g.emplace_back(n);
        x.back().zero();
        x.back()[std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin()] = 1; // initial guess
        *std::min_element(diagonals.begin(), diagonals.end()) = 1e99;
      }
      std::vector<std::vector<scalar>> Pcoeff(solver.m_roots);
      int nwork = solver.m_roots;
      for (auto iter = 0; iter < 10; iter++) {
        action(x, g);
        nwork = solver.addVector(x, g, Pcoeff);
        solver.report();
        if (nwork == 0)
          break;
        update(x, g, solver.working_set_eigenvalues());
      }
      std::cout << "Error={ ";
      for (const auto& e : solver.errors())
        std::cout << e << " ";
      std::cout << "} after " << solver.statistics().iterations << " iterations" << std::endl;
      for (size_t root = 0; root < solver.m_roots; root++) {
        std::cout << "Eigenvalue " << std::fixed << std::setprecision(9) << solver.eigenvalues()[root] << std::endl;
        solver.solution(root, x.front(), g.front(), Pcoeff.front());
        std::cout << "Eigenvector: (norm=" << std::sqrt(x[0].dot(x[0])) << "): ";
        for (size_t k = 0; k < n; k++)
          std::cout << " " << (x[0])[k];
        std::cout << std::endl;
      }
    }
  }
}
