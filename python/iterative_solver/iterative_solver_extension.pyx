import sys

import numpy as np
cimport numpy as np


# distutils : language = c++
m_mpicomm_compute = None

cdef extern from "../src/molpro/linalg/IterativeSolverC.h":
    ctypedef void (*apply_on_p_t)(const double*, double*, const size_t, const size_t*)
    size_t IterativeSolverAddP(size_t buffer_size, size_t nP, const size_t* offsets, const size_t* indices,
                             const double* coefficients, const double* pp, double* parameters, double* action,
                             int sync, apply_on_p_t apply_on_p)



current_problem = None

from libcpp.vector cimport vector
cdef void apply_on_p(const double* coefficients, double* action, const size_t size, const size_t* ranges) noexcept:
    cdef size_t size_ = size
    coefficients_ = np.asarray(<np.float64_t[:current_problem.p_space.size*size_]> coefficients).reshape((size_,current_problem.p_space.size))
    action_ = np.asarray(<np.float64_t[:current_problem.size_*size_]> action).reshape((size_,current_problem.size_))
    current_problem.p_action(coefficients_, action_, [ranges[0],ranges[1]])

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

    def add_P(self, pp, parameters, action):
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] action_ = action.reshape([parameters.shape[-1]*nbuffer])
        cdef size_t nbuffer_ = nbuffer
        cdef size_t nP_ = self.problem.p_space.size
        cdef double[::1] pp_ = pp.reshape([nP_*nP_])
        cdef size_t[::1] offsets_ = np.array(self.problem.p_space.offsets,dtype=np.uintp)
        cdef size_t[::1] indices_ = np.array(self.problem.p_space.indices,dtype=np.uintp)
        cdef double[::1] coefficients_ = np.array(self.problem.p_space.coefficients)
        global current_problem
        current_problem = self.problem
        current_problem.size_ = parameters.shape[-1]
        IterativeSolverAddP(nbuffer_, nP_, &offsets_[0], &indices_[0], &coefficients_[0], &pp_[0], &parameters_[0], &action_[0], 1, apply_on_p)
        return nbuffer


    def end_iteration(self, parameters, residual, sync=True):
        nbuffer = parameters.shape[0] if len(parameters.shape) > 1 else 1
        cdef double[::1] parameters_ = parameters.reshape([parameters.shape[-1]*nbuffer])
        cdef double[::1] residual_ = residual.reshape([parameters.shape[-1]*nbuffer])
        cdef size_t nbuffer_ = nbuffer
        result = IterativeSolverEndIteration(nbuffer_, &parameters_[0], &residual_[0], 1 if sync else 0)
        return int(result)

    @property
    def end_iteration_needed(self):
        result = IterativeSolverEndIterationNeeded()
        return result !=0

    @property
    def errors(self):
        e = np.ndarray(self.nroot)
        cdef double[::1] e_ = e
        IterativeSolverErrors(&e_[0])
        return e

    @property
    def converged(self):
        return IterativeSolverConverged() != 0
        pass

    def solve(self, parameters, actions, problem, generate_initial_guess=False, max_iter=None, max_p=0):
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
        if len(parameters.shape) <2 or len(actions.shape) < 2:
            parameters_reshape = parameters.reshape([self.nroot, self.n])
            actions_reshape = actions.reshape([self.nroot, self.n])
            return self.solve(parameters_reshape, actions_reshape, problem, generate_initial_guess, max_iter)
        self.iterations = 0
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
        if IterativeSolverNonLinear() == 0 and max_p>0 and problem.p_space.size == 0:
            if max_p >= self.nroot:
                for i in range(min(max_p,self.n)):
                    argmin = np.argmin(actions[0,:])
                    actions[0,argmin] = sys.float_info.max
                    problem.p_space.add_simple([argmin])
        if generate_initial_guess and problem.p_space.size == 0:
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
        for self.iterations in range(IterativeSolverMaxIter()):
            if IterativeSolverNonLinear() > 0:
                value = problem.residual(parameters.reshape([parameters.shape[-1]]), actions.reshape([parameters.shape[-1]]))
                if type(self) == Optimize:
                    nwork = self.add_value(value, parameters, actions)
                else:
                    nwork = self.add_vector(parameters[0,:], actions[0,:])
            elif self.iterations == 0 and problem.p_space.size> 0:
               self.problem = problem # TODO ? needed
               pp=problem.pp_action_matrix()
               nwork = self.add_P(pp,parameters,actions)
            else:
                problem.action(parameters,  actions)
                nwork = self.add_vector(parameters, actions)
            while self.end_iteration_needed:
                if nwork > 0:
                    IterativeSolverWorkingSetEigenvalues(&ev_[0])
                    if use_diagonals:
                        IterativeSolverDiagonals(&parameters_[0])
                        if type(self) == LinearEigensystem:
                            problem.precondition(actions[:nwork,:], shift=ev[:nwork], diagonals=parameters.reshape([parameters.size])[:parameters.shape[-1]])
                        else:
                            problem.precondition(actions[:nwork,:], diagonals=parameters.reshape([parameters.size])[:parameters.shape[-1]])
                    else:
                        problem.precondition(actions[:nwork,:], shift=ev[:nwork])
                nwork = self.end_iteration(parameters,actions)
            IterativeSolverErrors(&errors_[0])
            self.value = IterativeSolverValue()
            if IterativeSolverHasValues() != 0:
                reported = problem.report(self.iterations+1 if nwork > 0 else 0, verbosity, errors, value=value)
            elif IterativeSolverHasEigenvalues() != 0:
                IterativeSolverEigenvalues(&ev_[0])
                reported = problem.report(self.iterations+1 if nwork > 0 else 0, verbosity, errors, eigenvalues=ev[:self.nroot])
            else:
                reported = problem.report(self.iterations+1 if nwork > 0 else 0, verbosity, errors)
            if not reported and verbosity >=2:
                print('Iteration', self.iterations,'log10(|residual|)=', np.log10(errors))
                if IterativeSolverHasValues() != 0:
                    print('Objective function value', value)
            if nwork < 1: break
        return self.converged
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
