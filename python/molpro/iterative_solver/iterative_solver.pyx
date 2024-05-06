import numpy as np
cimport numpy as np
import cython
from libcpp.vector cimport vector

# distutils : language = c++
m_mpicomm_compute = None
class IterativeSolver:
    # cdef CppIterativeSolver* thisptr # thisptr is a pointer that will hold to the instance of the C++ classptr
    # def __cinit__(self):  # defines the python wrapper class' init function
        # cdef double[::1] InputArrayC = InputArray # defines a memoryview containnig a 1D numpy array - this can be passed as a C-like array by providing a pointer to the 1st element and the length
        # print('in python constructor, mapped:',mapped)
        # self.thisptr = new CppTestClass(InputArray.size, &InputArrayC[0], 1 if mapped else 0) # creates an instance of the C++ class and puts allocates the pointer to this
        # self.nroot = 1
        # self.nq = n
        # print('IterativeSolver.__cinit__')
    def __init__(self, n, nroot=1):
        self.n = n
        self.nroot = nroot

    def mpicomm_compute(self):
        global m_mpicomm_compute
        if m_mpicomm_compute is None:
            m_mpicomm_compute = IterativeSolver_mpicomm_global()
        return m_mpicomm_compute

    def solution(self,roots, parameters, residual, sync=True):
        # from cpython cimport array
        # cimport array
        cdef double[::1] parameters_ = parameters
        cdef double[::1] residual_ = residual
        cdef vector[int] roots_ = roots
        IterativeSolverSolution(len(roots),&roots_[0],&parameters_[0], &residual_[0], 1 if sync else 0)
        return self.value

    def add_value(self, value, parameters, action, sync=True):
        cdef double[::1] parameters_ = parameters
        cdef double[::1] action_ = action
        result = IterativeSolverAddValue(value, &parameters_[0], &action_[0], 1 if sync else 0)
        self.value = value
        return result

    def add_vector(self, buffer_size, parameters, action, sync=True):
        cdef double[::1] parameters_ = parameters
        cdef double[::1] action_ = action
        cdef size_t buffer_size_ = buffer_size
        result = IterativeSolverAddVector(buffer_size_, &parameters_[0], &action_[0], 1 if sync else 0)
        return result

    def end_iteration(self, solution, residual, sync=True):
        cdef double[::1] solution_ = solution
        cdef double[::1] residual_ = residual
        nbuffer = solution.shape[1] if len(solution.shape) > 1 else 1
        cdef size_t nbuffer_ = nbuffer
        result = IterativeSolverEndIteration(nbuffer_, &solution_[0], &residual_[0], 1 if sync else 0)
        return int(result)

    def solve(self, parameters, actions, problem, generate_initial_guess=False, max_iter=None):
        cdef double[::1] parameters_ = parameters
        cdef double[::1] actions_ = actions
        nbuffer = parameters.shape[1] if len(parameters.shape) > 1 else 1
        ev = np.array([self.nroot],dtype=float)
        cdef double[::1] ev_ = ev
        errors = np.array([self.nroot],dtype=float)
        cdef double[::1] errors_ = errors
        verbosity = IterativeSolverVerbosity()
        use_diagonals = problem.diagonals(actions)
        cdef int max_iter_
        if max_iter is not None:
            max_iter_ = max_iter
            IterativeSolverSetMaxIter(max_iter_)
        if use_diagonals:
            IterativeSolverSetDiagonals(&actions_[0])
        if generate_initial_guess:
            if not use_diagonals:
                raise ValueError('Default initial guess requested, but diagonal elements are not available')
            raise NotImplementedError('Default initial guess not yet implemented')

        nwork = nbuffer
        for iter in range(IterativeSolverMaxIter()):
            if IterativeSolverNonLinear() > 0:
                value = problem.residual(parameters, actions)
                nwork = self.add_value(value, parameters, actions)
            else:
                problem.action(parameters,  actions)
                nwork = self.add_vector(nbuffer, parameters, actions)
            if nwork > 0:
                IterativeSolverWorkingSetEigenvalues(&ev_[0])
                if use_diagonals:
                    IterativeSolverDiagonals(&parameters_[0])
                    problem.precondition(actions, ev, parameters)
                else:
                    problem.precondition(actions, ev)
            nwork = self.end_iteration(parameters,actions)
            if nwork <= 0: verbosity += 1
            IterativeSolverErrors(&errors_[0])
            self.value = IterativeSolverValue()
            if IterativeSolverHasValues() != 0:
                reported = problem.report(iter, verbosity, errors, value=value)
            elif IterativeSolverHasEigenvalues() != 0:
                IterativeSolverEigenvalues(&ev_[0])
                reported = problem.report(iter, verbosity, errors, eigenvalues=ev)
            else:
                reported = problem.report(iter, verbosity, errors)
            if not reported and verbosity >=2:
                print('Iteration', iter,'log10(|residual|)=', np.log10(errors))
                if IterativeSolverHasValues() != 0:
                    print('Objective function value', value)
            if nwork < 1: break

class Optimize(IterativeSolver):
    def __init__(self, n, range=None, thresh=1e-10, thresh_value=1e50, verbosity=0, minimize=True,
                  pname='', mpicomm=None, algorithm='', options=''):
        super().__init__(n)
        cdef size_t n_ = n
        cdef size_t range_[2]
        if range is None:
            range_ = [0, 0]
        else:
            range_ = range
        cdef size_t * rb = &range_[0]
        cdef size_t * re = &range_[1]
        cdef double thresh_ = thresh
        cdef double thresh_value_ = thresh_value
        cdef int verbosity_ = verbosity
        cdef bytes pname__ = pname.encode()
        cdef char * pname_ = pname__
        cdef int mpicomm_ = mpicomm if mpicomm is not None else self.mpicomm_compute()
        cdef bytes algorithm__ = algorithm.encode()
        cdef char * algorithm_ = algorithm__
        cdef bytes options__ = options.encode()
        cdef char * options_ = options__
        IterativeSolverOptimizeInitialize(n_, rb, re, thresh_, thresh_value_, verbosity_,
                                          1 if minimize else 0, pname_, mpicomm_, algorithm_,
                                          options_)
        if range is not None:
            range[0] = rb[0]
            range[1] = re[0]

    # def __init__(self):
    #     print('Optimize.__init__()')
