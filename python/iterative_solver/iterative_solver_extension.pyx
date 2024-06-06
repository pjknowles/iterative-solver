import sys

import numpy as np
cimport numpy as np
import cython
from libcpp.vector cimport vector

# distutils : language = c++
m_mpicomm_compute = None
class IterativeSolver:
    '''
    Base class inherited by Optimize, NonLinearEquations, LinearEquations, LinearEigensystem
    '''
    def __init__(self, n, nroot=1):
        '''

        :param n: The dimension of the parameter space
        :type n: int
        :param nroot: For the case of linear methods, the number of solutions that will be sought
        :type nroot: int
        '''
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
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*len(roots)])
        cdef double[::1] residual_ = residual.reshape([parameters.shape[-1]*len(roots)])
        cdef vector[int] roots_ = roots
        IterativeSolverSolution(len(roots),&roots_[0],&parameters_[0], &residual_[0], 1 if sync else 0)
        return self.value

    def add_value(self, value, parameters, action, sync=True):
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] action_ = action.reshape([parameters.shape[-1]*nbuffer])
        result = IterativeSolverAddValue(value, &parameters_[0], &action_[0], 1 if sync else 0)
        self.value = value
        return result

    def add_vector(self, parameters, action, sync=True):
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] action_ = action.reshape([parameters.shape[-1]*nbuffer])
        cdef size_t nbuffer_ = nbuffer
        result = IterativeSolverAddVector(nbuffer_, &parameters_[0], &action_[0], 1 if sync else 0)
        return result

    def end_iteration(self, parameters, residual, sync=True):
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] residual_ = residual.reshape([parameters.shape[-1]*nbuffer])
        cdef size_t nbuffer_ = nbuffer
        result = IterativeSolverEndIteration(nbuffer_, &parameters_[0], &residual_[0], 1 if sync else 0)
        return int(result)
    @property
    def errors(self):
        e = np.ndarray(self.nroot)
        cdef double[::1] e_ = e
        IterativeSolverErrors(&e_[0])
        return e
        pass

    def solve(self, parameters, actions, problem, generate_initial_guess=False, max_iter=None):
        '''
        Driver for iterating towards problem solution.
        :param parameters:
        :type parameters:
        :param actions:
        :type actions:
        :param problem:
        :type problem:
        :param generate_initial_guess:
        :type generate_initial_guess:
        :param max_iter:
        :type max_iter:
        :return:
        :rtype:
        '''
        # print ('solve:',parameters.shape, actions.shape, self.nroot, self.n)
        if len(parameters.shape) <2 or len(actions.shape) < 2:
            parameters_reshape = parameters.reshape([self.nroot, self.n])
            actions_reshape = actions.reshape([self.nroot, self.n])
            return self.solve(parameters_reshape, actions_reshape, problem, generate_initial_guess, max_iter)
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] actions_ = actions.reshape([parameters.shape[-1]*nbuffer])
        ev = np.zeros([self.nroot],dtype=float)
        cdef double[::1] ev_ = ev
        errors = np.array([self.nroot],dtype=float)
        cdef double[::1] errors_ = errors
        verbosity = IterativeSolverVerbosity()
        use_diagonals = problem.diagonals(actions.reshape([actions.size]))
        cdef int max_iter_
        if max_iter is not None:
            max_iter_ = max_iter
            IterativeSolverSetMaxIter(max_iter_)
        if use_diagonals:
            IterativeSolverSetDiagonals(&actions_[0])
        if generate_initial_guess:
            parameters[:,:] = 0
            if type(self) == LinearEigensystem:
                if not use_diagonals:
                    raise ValueError('Default initial guess requested, but diagonal elements are not available')
                IterativeSolverDiagonals(&actions_[0])
                for i in range(self.nroot):
                    argmin = np.argmin(actions[0,:])
                    actions[0,argmin] = sys.float_info.max
                    parameters[i,argmin] = 1.0
            elif type(self) == LinearEquations:
                for i in range(self.nroot):
                    parameters[i,i]=1


        nwork = nbuffer
        for iter in range(IterativeSolverMaxIter()):
            # print('start of iteration')
            # print('parameters',parameters)
            if IterativeSolverNonLinear() > 0:
                value = problem.residual(parameters.reshape([parameters.shape[-1]]), actions.reshape([parameters.shape[-1]]))
                nwork = self.add_value(value, parameters, actions)
            else:
                problem.action(parameters,  actions)
                nwork = self.add_vector(parameters[:nwork,:], actions[:nwork,:])
                # print('actions',actions)
            if nwork > 0:
                IterativeSolverWorkingSetEigenvalues(&ev_[0])
                if use_diagonals:
                    IterativeSolverDiagonals(&parameters_[0])
                    problem.precondition(actions[:nwork,:], shift=ev[:nwork], diagonals=parameters.reshape([parameters.size])[:parameters.shape[-1]])
                else:
                    problem.precondition(actions[:nwork,:], shift=ev[:nwork])
            # print('before end_iteration')
            # print('parameters',parameters, parameters.shape)
            # print('actions',actions, actions.shape)
            nwork = self.end_iteration(parameters,actions)
            # print('after end_iteration', nwork)
            # print('parameters',parameters, parameters.shape)
            # print('actions',actions, actions.shape)
            if nwork <= 0: verbosity += 1
            IterativeSolverErrors(&errors_[0])
            # print('after errors',nwork,errors)
            self.value = IterativeSolverValue()
            if IterativeSolverHasValues() != 0:
                reported = problem.report(iter, verbosity, errors, value=value)
            elif IterativeSolverHasEigenvalues() != 0:
                IterativeSolverEigenvalues(&ev_[0])
                reported = problem.report(iter, verbosity, errors, eigenvalues=ev[:self.nroot])
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

class NonLinearEquations(IterativeSolver):
    def __init__(self, n, range=None, thresh=1e-10, verbosity=0,
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
        cdef int verbosity_ = verbosity
        cdef bytes pname__ = pname.encode()
        cdef char * pname_ = pname__
        cdef int mpicomm_ = mpicomm if mpicomm is not None else self.mpicomm_compute()
        cdef bytes algorithm__ = algorithm.encode()
        cdef char * algorithm_ = algorithm__
        cdef bytes options__ = options.encode()
        cdef char * options_ = options__
        IterativeSolverNonLinearEquationsInitialize(n_, rb, re, thresh_, verbosity_,
                                          pname_, mpicomm_, algorithm_,
                                          options_)
        if range is not None:
            range[0] = rb[0]
            range[1] = re[0]

class LinearEquations(IterativeSolver):
    def __init__(self, rhs, range=None, aughes=0.0, thresh=1e-10, thresh_value=1e50, hermitian=False, verbosity=0,
                 pname='', mpicomm=None, algorithm='', options=''):
        n = rhs.shape[-1]
        nroot = rhs.shape[0] if len(rhs.shape)>1 else 1
        super().__init__(n, nroot)
        cdef size_t n_ = n
        cdef size_t nroot_ = nroot
        cdef double[::1] rhs_ = rhs.reshape([n*nroot])
        cdef size_t range_[2]
        if range is None:
            range_ = [0, 0]
        else:
            range_ = range
        cdef size_t * rb = &range_[0]
        cdef size_t * re = &range_[1]
        cdef double aughes_ = aughes
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
        IterativeSolverLinearEquationsInitialize(n_, nroot_, rb, re, &rhs_[0], aughes_, thresh_, thresh_value_,
                                                   1 if hermitian else 0, verbosity_,
                                                   pname_, mpicomm_, algorithm_,
                                                   options_)
        if range is not None:
            range[0] = rb[0]
            range[1] = re[0]

class LinearEigensystem(IterativeSolver):
    def __init__(self, n, nroot, range=None, thresh=1e-10, thresh_value=1e50, hermitian=False, verbosity=0,
                 pname='', mpicomm=None, algorithm='', options=''):
        super().__init__(n, nroot)
        cdef size_t n_ = n
        cdef size_t nroot_ = nroot
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
        IterativeSolverLinearEigensystemInitialize(n_, nroot_, rb, re, thresh_, thresh_value_,
                                                 1 if hermitian else 0, verbosity_,
                                                 pname_, mpicomm_, algorithm_,
                                                 options_)
        if range is not None:
            range[0] = rb[0]
            range[1] = re[0]
    @property
    def eigenvalues(self):
        e = np.ndarray(self.nroot)
        cdef double[::1] e_ = e
        IterativeSolverEigenvalues(&e_[0])
        return e
