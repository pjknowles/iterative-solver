import unittest
import iterative_solver
import numpy as np


class TestCase(unittest.TestCase):
    class RayleighQuotient(iterative_solver.Problem):
        def __init__(self, n, rho=0.1, **kwargs):
            self.size = n
            self.rho = rho

        @property
        def matrix(self):
            # M_ij = (i+1) delta_ij + rho
            m = self.rho * np.ones((self.size, self.size))
            for i in range(self.size):
                m[i, i] += i + 1
            return m

        def residual(self, parameters, gradient):
            # f = (sum_i (i+1)*x[i]^2 + rho*sum_ij x[i]*x[j])/(sum_i x[i]**2)
            # if len(parameters.shape) == 2:
            #     return self.residual(parameters[0, :], gradient[0, :])
            assert len(parameters.shape) == 1
            self.action(parameters.reshape([1, parameters.size]), gradient.reshape([1, parameters.size]))
            f = np.dot(parameters, gradient) / np.dot(parameters, parameters)
            for i in range(gradient.size):
                gradient[i] = 2 * (gradient[i] - f * parameters[i]) / np.dot(parameters, parameters)
            return f

        def action(self, parameters, residual):
            assert len(parameters.shape) == 2
            np.matmul(parameters, self.matrix, out=residual)

        def diagonals(self, diagonals):
            for i in range(self.size):
                diagonals[i] = self.matrix[i, i] - 0  # TODO underlying bug in that 1 is added to diagonals
            return True

    def test_problem(self):
        problem = TestCase.RayleighQuotient(4, 0.01)
        parameters = np.ones(problem.size) * 77
        residual = np.zeros(problem.size)
        step = 1e-5
        parameters[0] += step
        f1 = problem.residual(parameters, residual)
        parameters[0] -= 2 * step
        fm1 = problem.residual(parameters, residual)
        parameters[0] += step
        f0 = problem.residual(parameters, residual)
        # print(f0,f1, fm1)
        # print((f1-fm1)/(2*step),residual[0])
        self.assertAlmostEqual((f1 - fm1) / (2 * step), residual[0])

    def test_optimize(self):
        problem = TestCase.RayleighQuotient(4, 0.01)
        parameters = np.zeros(problem.size);
        parameters[0] = 1
        residual = np.zeros(problem.size)

        solver = iterative_solver.Optimize(problem.size, verbosity=0)
        self.assertEqual(type(solver), iterative_solver.Optimize)
        # print('initial parameters', parameters)
        solver.solve(parameters, residual, problem)
        answer = solver.solution([0], parameters, residual)
        # print('final parameters', parameters)
        # print('final function value', answer)
        # print('final residual', residual)
        for i in range(problem.size):
            self.assertAlmostEqual(residual[i], 0.0)

    def test_diagonalize(self):
        problem = TestCase.RayleighQuotient(8, 0.1)
        # print(problem.matrix)
        nroot = 2
        parameters = np.zeros([nroot, problem.size])
        for i in range(nroot):
            parameters[i, i] = 1
        residual = np.ndarray([nroot, problem.size])

        solver = iterative_solver.LinearEigensystem(problem.size, nroot, verbosity=0)
        self.assertEqual(type(solver), iterative_solver.LinearEigensystem)
        solver.solve(parameters, residual, problem)
        solver.solution([i for i in range(nroot)], parameters, residual)
        # print('final parameters', parameters)
        # print('final residual', residual)
        # print('final eigenvalues', solver.eigenvalues)
        # print('final errors', solver.errors)
        self.assertEqual(solver.errors.size, nroot)
        for e in solver.errors:
            self.assertAlmostEqual(e, 0.0)
        for root in range(nroot):
            for i in range(problem.size):
                self.assertAlmostEqual(residual[root, i], 0.0)
            # print('residual before residual function', residual[root])
            value = problem.residual(parameters[root], residual[root])
            self.assertAlmostEqual(value, solver.eigenvalues[root])
            # print('residual from residual function', residual[root])
            for i in range(problem.size):
                self.assertAlmostEqual(residual[root, i], 0.0)


if __name__ == '__main__':
    unittest.main()
