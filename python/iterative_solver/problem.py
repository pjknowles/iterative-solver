import numpy as np
import sys

from .pspace import PSpace


class Problem:
    '''
    Semi-abstract base class for specifying the problem to be solved by iterative_solver.
    '''

    def __init__(self):
        self.p_space = PSpace()
        self.dimension = None

    def residual(self, parameters, residual):
        """
        Calculate the residual vector. Used by non-linear solvers (NonLinearEquations,  Optimize) only.
        :return: In the case where the residual is an exact differential, the corresponding function value. Used by Optimize but not NonLinearEquations.
        :rtype: float
        :param parameters: The trial solution for which the residual is to be calculated
        :type parameters: np.ndarray(dtype=float)
        :param residual:  The residual vector
        :type residual: np.ndarray(dtype=float)
        """
        raise NotImplementedError

    def action(self, parameters, action):
        """
        Calculate the action of the kernel matrix on a set of parameters. Used by
        linear solvers, but not by the non-linear solvers (NonLinearEquations, Optimize).

        :param parameters: The trial solutions for which the action is to be calculated
        :type parameters: np.ndarray(dtype=float)
        :param action: The action vectors
        :type action: np.ndarray(dtype=float)
        """
        raise NotImplementedError

    def diagonals(self, diagonals):
        """
        Optionally provide the diagonal elements of the underlying kernel. If
        implemented and returning true, the provided diagonals will be used by
        IterativeSolver for preconditioning (and therefore the precondition() function does
        not need to be implemented), and, in the case of linear problems, for selection of
        the P space. Otherwise, preconditioning will be done with precondition(), and any P
        space has to be provided manually.

        :param diagonals: On exit, contains the diagonal elements
        :type diagonals: np.ndarray(dtype=bool, ndim=1)
        :return: Whether diagonals have been provided.
        :rtype: bool
        """
        return False

    def precondition(self, residual, shift=None, diagonals=None):
        """
        Apply preconditioning to a residual vector in order to predict a step towards
        the solution
   */

        :param residual: On entry, assumed to be the residual. On exit, the negative of the predicted step.
        :type residual: np.ndarray(dtype=float)
        :param shift: When called from LinearEigensystem, contains the corresponding current
           eigenvalue estimates for each of the parameter vectors in the set. All other solvers
           pass a vector of zeroes.
        :type shift: float
        :param diagonals: Diagonal elements of kernel
        :type diagonals: np.ndarray(dtype=float, ndim=1)
        """
        small = 1e-14
        if len(residual.shape) > 1:
            for i in range(residual.shape[0]):
                if shift is not None:
                    assert type(shift) == np.ndarray
                self.precondition(residual[i, :], float(shift[i]) if shift is not None else None, diagonals)
            return
        if diagonals is not None:
            if shift is not None:
                for j in range(residual.size):
                    residual[j] = residual[j] / (diagonals[j] - shift + small)
            else:
                for j in range(residual.size):
                    residual[j] = residual[j] / (diagonals[j] + small)
        else:
            raise NotImplementedError
    def pp_action_matrix(self):
        """
        Calculate the representation of the kernel matrix in the P space. Implementation required only for linear hermitian problems for which P-space acceleration is wanted.
        """
        if self.p_space.size > 0:
            raise NotImplementedError('P-space unavailable: unimplemented pp_action_matrix() in Problem class')
        return np.array([], dtype=np.double)

    def p_action(self, p_coefficients, actions, ranges):
        """
        Calculate the action of the kernel matrix on a set of vectors in the P space. Implementation required only for linear hermitian problems for which P-space acceleration is wanted.

        :param p_coefficients The projection of the vectors onto to the P space
        :param actions On exit, the computed action has been added to the original contents
        :param range The range of the full space for which actions should be computed.
        """
        if self.p_space.size > 0:
            raise NotImplementedError('P-space unavailable: unimplemented p_action() in Problem class')

    def test_parameters(self, instance, parameters):
        return False

    def report(self, iteration, verbosity, errors, value=None, eigenvalues=None):
        if (iteration <= 0 and verbosity >= 1) or verbosity >= 2:
            if iteration > 0 and verbosity >= 2:
                print('Iteration', iteration, 'log10(|residual|)=', np.log10(errors))
            elif iteration == 0:
                print('Converged', 'log10(|residual|)=', np.log10(errors + sys.float_info.min))
            else:
                print('Unconverged', 'log10(|residual|)=', np.log10(errors + sys.float_info.min))
            if value is not None:
                print('Objective function value', value)
            if eigenvalues is not None:
                print('Eigenvalues', eigenvalues)
            return True
        else:
            return False
