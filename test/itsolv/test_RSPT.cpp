#include "test.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <molpro/iostream.h>
#include <numeric>
#include <regex>
#include <vector>

#include "vector_types.h"
#include <molpro/linalg/itsolv/CastOptions.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <molpro/linalg/itsolv/helper.h>

// Find lowest eigensolution of a matrix obtained from an external file
// Storage of vectors in-memory via class Rvector

using molpro::linalg::array::Span;
using molpro::linalg::itsolv::CastOptions;
using molpro::linalg::itsolv::CVecRef;
using molpro::linalg::itsolv::cwrap;
using molpro::linalg::itsolv::VecRef;
using molpro::linalg::itsolv::wrap;

struct RSPT : ::testing::Test {
    using MatrixXdr = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MatrixXdc = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using vectorP = std::vector<scalar>;

    size_t n = 0;
    size_t verbosity = 0;
    MatrixXdc hmat;
    Eigen::VectorXd h0;

    void load_matrix(const std::string &file, double degeneracy_split = 0) {
        std::ifstream f(std::string{"./"} + file + ".hamiltonian");
        f >> n;
        hmat.resize(n, n);
        for (auto i = 0; i < n; i++)
            for (auto j = 0; j < n; j++)
                f >> hmat(i, j);
        //    std::cout << "hmat\n" << hmat << std::endl;
        // split degeneracies
        for (auto i = 0; i < n; i++)
            hmat(i, i) += degeneracy_split * i;
        //    std::cout << "hmat\n" << hmat << std::endl;
        h0.resize(n);
        std::ifstream f0(std::string{"./"} + file + ".h0");
        for (auto i = 0; i < n; i++)
            f0 >> h0(i);
        //    std::cout << "h0 " << h0.adjoint() << std::endl;
    }

    void action(const Rvector &psx, Rvector &outputs) {
        auto x = Eigen::Map<const Eigen::VectorXd>(psx.data(), psx.size());
        auto r = Eigen::Map<Eigen::VectorXd>(outputs.data(), psx.size());
        r = hmat * x;
    }

    template<class scalar>
    scalar dot(const std::vector<scalar> &a, const std::vector<scalar> &b) {
        return std::inner_product(std::begin(a), std::end(a), std::begin(b), scalar(0));
    }

    void update(Rvector &psg) {
        double e0 = 1e50;
        for (size_t i = 0; i < n; i++)
            if (h0(i) < e0)
                e0 = h0(i);
        for (size_t i = 0; i < n; i++) {
            psg[i] = psg[i] / (1e-12 - e0 + h0(i));
        }
    }

    auto initial_guess(Rvector &x) {
        std::vector<double> diagonals;
        diagonals.reserve(n);
        for (auto i = 0; i < n; i++)
            diagonals.push_back(h0(i));
        x[(std::min_element(diagonals.begin(), diagonals.end()) - diagonals.begin())] = 1; // initial guess
    }

    void test_eigen(const std::string &title = "", const int n_working_vectors_max = 0) {
        {
            molpro::cout << "\n\n*** " << title << " by perturbation theory" << std::endl;
            auto solver = molpro::linalg::itsolv::create_LinearEigensystem<Rvector, Qvector, Pvector>("RSPT");
            int nwork = 1;
            Rvector x(n, 0);
            Rvector g(n, 0);
            initial_guess(x);
            for (auto iter = 1; iter < 10; iter++) {
                action(x, g);
//        std::cout << "x " << x << std::endl;
//        std::cout << "g " << g << std::endl;
                nwork = solver->add_vector(x, g);
//        std::cout << "g after add_vector " << g << std::endl;
                update(g);
//        std::cout << "g after update " << g << std::endl;
//        std::cout << "x after update " << x << std::endl;
                nwork = solver->end_iteration(x, g);
//        std::cout << "g after end_iteration " << g << std::endl;
//        std::cout << "x after end_iteration " << x << std::endl;
                solver->report();
            }
        }
    }

    enum method_e {
        BFGS, DIIS
    };

    template<method_e method = BFGS>
    double
    test_Hylleraas(const std::string &title = "", const int n_working_vectors_max = 0, bool precondition = true,
                   double expected = double(0)) {
//        molpro::cout << "\n\n*** " << title << " Hylleraas minimum by " << (method == BFGS ? "BFGS" : "DIIS")
//                     << std::endl;
        auto optimizer = molpro::linalg::itsolv::create_Optimize<Rvector, Qvector, Pvector>("BFGS");
        auto nonlinear_solver = molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector, Pvector>("DIIS");
        molpro::linalg::itsolv::IterativeSolver<Rvector, Qvector, Pvector> *solver;
        if (method == BFGS)
            solver = dynamic_cast<molpro::linalg::itsolv::IterativeSolver<Rvector, Qvector, Pvector> *>(optimizer.get());
        if (method == DIIS)
            solver = dynamic_cast<molpro::linalg::itsolv::IterativeSolver<Rvector, Qvector, Pvector> *>(nonlinear_solver.get());
//            std::cout << "h\n" << hmat << std::endl;
//            std::cout << "h0\n" << h0 << std::endl;
        auto ham0 = hmat;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++)ham0(i, j) = 0;
            ham0(i, i) = h0(i);
        }
        auto ham1 = hmat - ham0;
//            std::cout << "ham0\n" << ham0 << std::endl;
//            std::cout << "ham1\n" << ham1 << std::endl;
        int nwork = 1;
        Rvector x(n, 0);
        Rvector x0(n, 0);
        Rvector g(n, 0);
        auto xv = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
        auto rv = Eigen::Map<Eigen::VectorXd>(g.data(), g.size());
        initial_guess(x0);
        auto xv0 = Eigen::Map<const Eigen::VectorXd>(x0.data(), x0.size());
        auto e0 = h0.dot(xv0);
        auto e1 = xv0.dot(ham1 * xv0);
//            std::cout << "x0=" << xv0.transpose() << std::endl;
//            std::cout << "e0=" << e0 << std::endl;
//            std::cout << "e1=" << e1 << std::endl;
        double e2;
        for (auto iter = 1; iter < 10; iter++) {
            e2 = 2 * xv0.dot(ham1 * xv - e1 * xv) + xv.dot(ham0 * xv - e0 * xv);
            rv = ham1 * xv0 - e1 * xv0 + ham0 * xv - e0 * xv;
//                std::cout << "x " << x << std::endl;
//                std::cout << "g " << g << std::endl;
//                std::cout << "e2 " << e2 << std::endl;
            if constexpr(method == BFGS)
                nwork = optimizer->add_value(x, e2 / 2, g);
            else
                nwork = nonlinear_solver->add_vector(x, g);
//        std::cout << "g after add_vector " << g << std::endl;
            if (precondition) update(g);
//        std::cout << "g after update " << g << std::endl;
//        std::cout << "x after update " << x << std::endl;
            nwork = solver->end_iteration(x, g);
//        std::cout << "g after end_iteration " << g << std::endl;
//        std::cout << "x after end_iteration " << x << std::endl;
//            solver->report();
            if (nwork < 1)break;
        }
        if (expected != 0)
            EXPECT_NEAR(e2, expected, 1e-11);
        return e2;
    }

};

TEST_F(RSPT, file_eigen) {
    for (const auto &file : std::vector<std::string>{"he", "bh", "hf"}) {
        load_matrix(file, file == "phenol" ? 0 : 1e-8);
        test_eigen(file);
    }
}

TEST_F(RSPT, file_Hylleraas_BFGS) {
    for (const auto &file : std::vector<std::string>{"he", "bh", "hf"}) {
        load_matrix(file, file == "phenol" ? 0 : 1e-8);
        auto expected = test_Hylleraas(file, 0, true, 0);
        test_Hylleraas<BFGS>(file, 0, false, expected);
        test_Hylleraas<DIIS>(file, 0, true, expected);
//        test_Hylleraas_BFGS<DIIS>(file,0,false, expected);

    }
}
