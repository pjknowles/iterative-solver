import unittest
import iterative_solver
import numpy as np


class MyTestCase(unittest.TestCase):
    class TestProblem(iterative_solver.Problem):
        def residual(self, parameters, residual):
            # f = sum_i (i+1)*(x[i]-1)^2
            f = 0.0
            for i in range(parameters.size):
                f += (i + 1) * (parameters[i]-1) ** 2
                residual[i] = 2 * (i + 1) * (parameters[i]-1)
            return f

        def diagonals(self, diagonals):
            for i in range(diagonals.size):
                diagonals[i] = 2 * (i + 1)
            return True

    def test_optimize(self):
        problem = MyTestCase.TestProblem()
        parameters = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        residual = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

        solver = iterative_solver.Optimize(parameters.size, verbosity=0)
        self.assertEqual(type(solver), iterative_solver.Optimize)  # add assertion here
        solver.solve(parameters, residual, problem)
        self.assertAlmostEqual(solver.solution([0], parameters, residual),0.0)
        print(parameters)
        print(residual)
        for i in range(parameters.size):
            self.assertAlmostEqual(residual[i],0.0)
            self.assertAlmostEqual(parameters[i],1.0)


if __name__ == '__main__':
    unittest.main()
