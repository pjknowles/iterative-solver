import unittest
import iterative_solver
import numpy as np


class TestCase(unittest.TestCase):
    class OptimisationProblem(iterative_solver.Problem):
        def __init__(self, n, **kwargs):
            self.size = n

        def residual(self, parameters, residual):
            # f = sum_i (i+1)*(x[i]-1)^2
            f = 0.0
            for i in range(self.size):
                f += (i + 1) * (parameters[i] - 1) ** 2
                residual[i] = 2 * (i + 1) * (parameters[i] - 1)
            return f

        def diagonals(self, diagonals):
            for i in range(self.size):
                diagonals[i] = 2 * (i + 1)
            return True

    def test_optimize(self):
        problem = TestCase.OptimisationProblem(1000)
        parameters = np.zeros(problem.size)
        residual = np.zeros(problem.size)

        solver = iterative_solver.Optimize(problem.size, verbosity=0)
        self.assertEqual(type(solver), iterative_solver.Optimize)  # add assertion here
        solver.solve(parameters, residual, problem)
        self.assertAlmostEqual(solver.solution([0], parameters, residual), 0.0)
        # print(parameters)
        # print(residual)
        for i in range(problem.size):
            self.assertAlmostEqual(residual[i], 0.0)
            self.assertAlmostEqual(parameters[i], 1.0)


if __name__ == '__main__':
    unittest.main()
