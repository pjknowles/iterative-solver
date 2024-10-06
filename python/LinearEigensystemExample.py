import iterative_solver
import numpy as np

n = 300
nroot = 3
m = np.array([[1 if i != j else 3 * (n - i) for i in range(n)] for j in range(n)])


def solve(m, nP, max_p, nroot):
    problem = iterative_solver.MatrixProblem()
    problem.attach(m)
    problem.p_space.add_simple(range(min(nP, m.shape[0])))  # the first nP components, so not the best
    c = np.empty((nroot, m.shape[0]), dtype=float)
    g = np.empty((nroot, m.shape[0]), dtype=float)
    iterative_solver.Solve_Linear_Eigensystem(c, g, problem, nroot, verbosity=2, thresh=1e-8, hermitian=True,
                                              max_p=max_p)


for nP in range(0, 51, 50):
    for max_p in range(0, 51, 10):
        if max_p > 0 and nP > 0: continue
        print('Explicit P-space =', nP, ', auto P-Space =', max_p, ', dimension =', m.shape[0], ', roots =', nroot)
        solve(m, nP, max_p, nroot)
