from iterative_solver import Problem
import numpy as np
import math


class Box(Problem):
    r"""
    Sum of exponentials

    M. J. Box, A comparison of several current optimization methods, and the use of transformations in constrained problems, Comput. J., v. 9, 1966, pp. 67-77.
    """

    def __init__(self):
        super().__init__()

    def residual(self, x, g):
        g[:] = 0.0
        value = 0.0
        for i in range(1, 11):
            t = i * 0.1
            d = math.exp(-x[0] * t) - math.exp(-x[1] * t) - math.exp(-t) + math.exp(-10 * t)
            value += d * d
            g[0] += 2 * d * (-t) * math.exp(-x[0] * t)
            g[1] -= 2 * d * (-t) * math.exp(-x[1] * t)
        return value

    def diagonals(self, diag):
        diag[:] = 1
        return True

    def hessian(self, x):
        # print('hessian',x)
        x0 = min(max(x[0], -10), 10)
        x1 = min(max(x[1], -10), 10)
        h = np.zeros([2, 2])
        for i in range(1, 11):
            t = i * 0.1
            d = math.exp(-x0 * t) - math.exp(-x1 * t) - math.exp(-t) + math.exp(-10 * t)
            g[0] -= 2 * d * t * math.exp(-x0 * t)
            g[1] += 2 * d * t * math.exp(-x1 * t)
            h[0, 0] += 2 * d * t * t * math.exp(-x0 * t) + 2 * t * t * math.exp(-2 * x0 * t)
            h[1, 1] += -2 * d * t * t * math.exp(-x1 * t) + 2 * t * t * math.exp(-2 * x1 * t)
            h[1, 1] += -2 * t * t * math.exp(-x0 * t - x1 * t)
        return h

    def test_parameters(self, instance, parameters):
        parameters[:] = .1
        parameters[instance % parameters.shape[0]] = .1 * (1 + instance // parameters.shape[0])
        return True

class Rosenbrock(Problem):
    r"""
    Generalised Rosenbrock function, as defined in scipy.optimize.

    Rosenbrock, H.H. (1960). An automatic method for finding the greatest or least value of a function. The Computer Journal. 3 (3): 175–184. doi:10.1093/comjnl/3.3.175
    """

    def __init__(self, a=1, b=100):
        super().__init__()
        self.a = a
        self.b = b

    def residual(self, x, g):
        value = 0.0
        g[:] = 0.0
        for i in range(x.shape[0] - 1):
            t1 = (self.a - x[i])
            t2 = (x[i + 1] - x[i] * x[i]);
            value += t1 * t1 + self.b * t2 * t2
            g[i] += -2 * t1 - 4 * self.b * t2 * x[i]
            g[i + 1] += 2 * self.b * t2
        return value

    def diagonals(self, diag):
        diag[:] = 200
        return True

    def hessian(self, x):
        h = np.zeros([x.shape[0], x.shape[0]])
        for i in range(x.shape[0] - 1):
            x0 = min(max(x[i], -10), 10)
            x1 = min(max(x[i+1], -10), 10)
            h[i, i] += 1200 * x0 * x0 - 400 * x1 + 1
            h[i, i + 1] += -400 * x0
            h[i + 1, i] += - 400 * x0
            h[i + 1, i + 1] += 200
        return h

    def test_parameters(self, instance, parameters):
        parameters[:] = .1
        parameters[instance%parameters.shape[0]] = .1 * (1+instance//parameters.shape[0])
        return True


