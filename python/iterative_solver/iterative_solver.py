from .problem import Problem
from .iterative_solver_extension import Optimize, LinearEigensystem, LinearEquations, NonLinearEquations
from .natural_coordinate import NaturalCoordinateProblem
import numpy as np


def Solve_Linear_Eigensystem(parameters, actions, problem: Problem, nroot=1, generate_initial_guess=True,
                             max_iter=1000000, max_p=0, **kwargs):
    solver = LinearEigensystem(parameters.shape[-1], nroot, **kwargs)
    solver.solve(parameters, actions, problem, generate_initial_guess=generate_initial_guess, max_iter=max_iter,
                 max_p=max_p)
    solver.solution(range(min(nroot, parameters.shape[0])), parameters, actions)
    return solver


def Solve_Linear_Equations(parameters, actions, problem: Problem, generate_initial_guess=True, max_iter=1000000,
                           max_p=0, **kwargs):
    solver = LinearEquations(parameters.shape[-1], **kwargs)
    solver.solve(parameters, actions, problem, generate_initial_guess=generate_initial_guess, max_iter=max_iter,
                 max_p=max_p)
    return solver


def Solve_NonLinear_Equations(parameters, actions, problem: Problem, generate_initial_guess=True, max_iter=1000000,
                              **kwargs):
    solver = NonLinearEquations(parameters.shape[-1], **kwargs)
    solver.solve(parameters, actions, problem, generate_initial_guess=generate_initial_guess, max_iter=max_iter)
    return solver


def Solve_Optimization(parameters, actions, problem: Problem, generate_initial_guess=True, max_iter=1000000, **kwargs):
    solver = Optimize(parameters.shape[-1], **kwargs)
    solver.solve(parameters, actions, problem, generate_initial_guess=generate_initial_guess, max_iter=max_iter)
    return solver


def test_problem_gradient(problem: Problem, dimension=1, step=1e-3, rtol=1e-7, atol=1e-13, instance=0):
    x = np.zeros(dimension)
    g_problem = np.zeros_like(x)
    g_numerical = np.zeros_like(x)
    if not hasattr(problem, 'test_parameters') or not problem.test_parameters(instance=instance, parameters=x):
        return False
    v0 = problem.residual(x, g_problem)
    for i in range(dimension):
        g = np.zeros_like(x)
        x[i] -= 2 * step
        vmm = problem.residual(x, g)
        x[i] += step
        vm = problem.residual(x, g)
        x[i] += 2 * step
        vp = problem.residual(x, g)
        x[i] += step
        vpp = problem.residual(x, g)
        x[i] -= 2 * step
        g_numerical[i] = (vmm - 8 * vm + 8 * vp - vpp) / (12 * step)
    np.testing.assert_allclose(g_problem, g_numerical, rtol=rtol, atol=atol)
