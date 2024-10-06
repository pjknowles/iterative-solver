import numpy as np
from iterative_solver import Problem
import copy
import sys


class MatrixProblem(Problem):
    '''
    A specialisation of the Problem class for linear problems in which the kernel matrix is stored in full in an existing array
    '''

    def __init__(self):
        super().__init__()

    def attach(self, matrix=None, RHS=None):
        if matrix is not None:
            self.matrix = matrix
        if RHS is not None:
            self.m_RHS = RHS

    def RHS(self, vector, instance):
        if instance < 0 or instance >= self.m_RHS.shape[1]:
            return False
        vector = copy.deepcopy(self.m_RHS[instance, :])  # check index order
        return True

    def action(self, parameters, action):
        """
        Calculate the action of the kernel matrix on a set of parameters. Used by
        linear solvers, but not by the non-linear solvers (NonLinearEquations, Optimize).

        :param parameters: The trial solutions for which the action is to be calculated
        :type parameters: np.ndarray(dtype=float)
        :param action: The action vectors
        :type action: np.ndarray(dtype=float)
        """
        np.copyto(action, np.matmul(parameters, self.matrix))

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
        np.copyto(diagonals[:self.matrix.shape[0]], self.matrix.diagonal())
        return True

    def pp_action_matrix(self):
        """
        Calculate the representation of the kernel matrix in the P space. Implementation required only for linear hermitian problems for which P-space acceleration is wanted.
        """
        matrix = np.zeros((self.p_space.size, self.p_space.size), dtype=np.double)
        for i in range(self.p_space.size):
            for j in range(self.p_space.size):
                matrix[i, j] = 0.0
                for ic in range(self.p_space.offsets[i], self.p_space.offsets[i + 1]):
                    for jc in range(self.p_space.offsets[j], self.p_space.offsets[j + 1]):
                        matrix[i, j] = matrix[i, j] + self.matrix[self.p_space.indices[ic], self.p_space.indices[jc]] * \
                                       self.p_space.coefficients[ic] * self.p_space.coefficients[jc]
        return matrix

    def p_action(self, p_coefficients, actions, ranges):
        """
        Calculate the action of the kernel matrix on a set of vectors in the P space. Implementation required only for linear hermitian problems for which P-space acceleration is wanted.

        :param p_coefficients The projection of the vectors onto to the P space
        :param actions On exit, the computed action has been added to the original contents
        :param range The range of the full space for which actions should be computed.
        """
        for i in range(actions.shape[0]):
            for k in range(self.p_space.size):
                for kc in range(self.p_space.offsets[k], self.p_space.offsets[k + 1]):
                    j = range(ranges[0], ranges[1])
                    actions[i, j] = actions[i, j] + self.matrix[self.p_space.indices[kc], j] * \
                                    self.p_space.coefficients[kc] * p_coefficients[i, k]
