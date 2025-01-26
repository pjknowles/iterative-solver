import unittest
import iterative_solver
from iterative_solver.example_problems import Rosenbrock
import numpy as np
import inspect


class TestCase(unittest.TestCase):

    @property
    def verbosity(self):
        frame = inspect.currentframe()
        while frame:
            self = frame.f_locals.get('self')
            if isinstance(self, unittest.TestProgram):
                return self.verbosity
            frame = frame.f_back
        return 0

    def test_difficult(self):
        x_initial = np.array([-3.67701376, -2.68322017])
        x_initial = np.array([-3.02080173, 3.60275987])
        x_initial = np.array([-2.20716128, -0.1547407])
        x_initial = [-0.63524643, -3.19916027]
        for x_initial in [
            [-3.67701376, -2.68322017],
            [-3.02080173, 3.60275987],
            [-2.20716128, -0.1547407],
            [-0.63524643, -3.19916027],
        ]:
            x = np.zeros_like(x_initial)
            g = np.zeros_like(x_initial)
            print(x_initial)
            problem = Rosenbrock()
            natural_coordinate_problem = iterative_solver.NaturalCoordinateProblem(problem, x_initial)
            natural_solver = iterative_solver.Solve_Optimization(x, g, natural_coordinate_problem, verbosity=5,
                                                                 options="max_size_qspace=1",
                                                                 max_iter=100)
            self.assertTrue(natural_solver.converged)


if __name__ == '__main__':
    unittest.main()
