import numpy as np
import sys
from . import iterative_solver


class NaturalCoordinateProblem(iterative_solver.Problem):
    r"""
    Wrapper problem class that implements optimization in natural coordinates.
    """

    def __init__(self, problem: iterative_solver.Problem, initial_parameters: np.ndarray, hessian_threshold=1e-8,
                 coordinate_refinements=1,
                 base_coordinate_inverse_algorithm=True, history=False):
        r"""

        :param problem: The problem class for the optimization in base coordinates. It should contain an additional method `hessian`, which returns an exact or approximate hessian matrix expressed as a numpy array
        :type problem: iterative_solver.Problem
        :param initial_parameters:
        :type initial_parameters: np.ndarray
        :param hessian_threshold: Hessian eigenvalues less than this value are assumed to represent natural coordinates that are redundant, and are explicitly excluded to avoid numerical noise
        :type hessian_threshold: float
        :param coordinate_refinements: How many refinements in the generation of base coordinates from natural coordinates
        :type coordinate_refinements: int
        :param base_coordinate_inverse_algorithm: Whether to use the inverse-average algorithm in the generation of base coordinates
        :type base_coordinate_inverse_algorithm: bool
        :param history: Whether to store the history of the optimization. If selected, the history is available as instance variable `history`
        :type history: bool
        """
        super().__init__()
        self.coordinate_refinements = coordinate_refinements
        self.base_coordinate_inverse_algorithm = base_coordinate_inverse_algorithm
        self.hessian_threshold = hessian_threshold
        if history:
            self.history = []
        self.reset(problem, initial_parameters)

    def reset(self, problem: iterative_solver.Problem, initial_parameters: np.ndarray):
        self.first = True
        self.problem = problem
        self.last_base_coordinates = np.array(initial_parameters, dtype=np.float64)
        if hasattr(self, 'history'):
            self.history = []

    def residual(self, natural_coordinates, residual):
        if self.first:
            base_coordinates_ = self.last_base_coordinates
        else:
            base_coordinates_ = self.base_coordinates(natural_coordinates, refinements=self.coordinate_refinements)
        hessian = self.problem.hessian(base_coordinates_)
        hessian_eigensolution = aligned_eigensolution(hessian,
                                                      reference_eigenvectors=None if self.first else self.last_hessian_eigenvectors)
        self.first = False

        base_coordinate_gradient = np.zeros_like(base_coordinates_)
        value = self.problem.residual(base_coordinates_, base_coordinate_gradient)
        self.hessian_diagonal = np.copy(hessian_eigensolution[0])
        for j in range(len(self.hessian_diagonal)):
            self.hessian_diagonal[j] = abs(self.hessian_diagonal[j])
            if self.hessian_diagonal[j] < self.hessian_threshold:
                self.hessian_diagonal[j] = sys.float_info.max
        hessian_eigenvectors = np.copy(hessian_eigensolution[1])
        residual[:] = hessian_eigenvectors.transpose() @ base_coordinate_gradient
        if hasattr(self, 'history'):
            self.history.append(
                {'natural_coordinates': np.copy(natural_coordinates),
                 'base_coordinates': np.copy(base_coordinates_),
                 'value': value,
                 'base_coordinate_gradient': np.copy(base_coordinate_gradient),
                 'hessian': np.copy(hessian),
                 'hessian_eigenvectors': np.copy(hessian_eigensolution[1]),
                 'hessian_eigenvalues': np.copy(hessian_eigensolution[0]),
                 'natural_coordinate_gradient': np.copy(residual),
                 }
            )
        if not hasattr(self, 'last_value') or value < self.last_value:
            self.last_value = value
            self.last_natural_coordinates = np.copy(natural_coordinates)
            self.last_base_coordinates = np.copy(base_coordinates_)
            self.last_hessian_eigenvectors = np.copy(hessian_eigenvectors)
        return value

    def base_coordinates(self, natural_coordinates, refinements=0):
        base_coordinates_ = self.last_base_coordinates + self.last_hessian_eigenvectors @ (
                natural_coordinates - self.last_natural_coordinates)
        for iter in range(refinements):
            hessian_eigensolution = aligned_eigensolution(self.problem.hessian(base_coordinates_),
                                                          reference_eigenvectors=self.last_hessian_eigenvectors)
            if self.base_coordinate_inverse_algorithm:
                base_coordinates_ = self.last_base_coordinates + 2 * np.linalg.inv(
                    (self.last_hessian_eigenvectors + hessian_eigensolution[1]).transpose()) @ (
                                            natural_coordinates - self.last_natural_coordinates)
            else:
                base_coordinates_ = self.last_base_coordinates + 0.5 * (
                        self.last_hessian_eigenvectors + hessian_eigensolution[1]) @ (
                                            natural_coordinates - self.last_natural_coordinates)
        return base_coordinates_

    def precondition(self, residual, shift=None, diagonals=None):
        if len(residual.shape) > 1:
            for i in range(residual.shape[0]):
                if shift is not None:
                    assert type(shift) == np.ndarray
                self.precondition(residual[i, :], float(shift[i]) if shift is not None else None, diagonals)
            return
        residual[:] = residual[:] / self.hessian_diagonal[:]


def aligned_eigensolution(matrix, reference_eigenvectors=None):
    r"""
    Find the eigensolutions of a symmetric matrix, then order and phase-align them with a set of reference vectors

    :param matrix:
    :type matrix: np.ndarray
    :param reference_eigenvectors:
    :type reference_eigenvectors: np.ndarray
    :return: eigenvalues, eigenvectors
    :rtype: (np.ndarray, np.ndarray)
    """
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if True:
        reconstructed_matrix = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.transpose()
        assert (np.isclose(reconstructed_matrix, matrix, rtol=1e-14)).all()
        temp = np.zeros_like(reconstructed_matrix)
        for i in range(eigenvalues.shape[0]):
            for j in range(eigenvalues.shape[0]):
                temp[i, j] = eigenvectors[i, j] * eigenvalues[j]
        reconstructed_matrix[:, :] = 0
        for i in range(eigenvalues.shape[0]):
            for j in range(eigenvalues.shape[0]):
                for k in range(eigenvalues.shape[0]):
                    reconstructed_matrix[i, j] += temp[i, k] * eigenvectors[j, k]
        assert (np.isclose(reconstructed_matrix, matrix, rtol=1e-14)).all()
        # print('eigenvectors after they are made\n', eigenvectors)
    if reference_eigenvectors is None:
        return (eigenvalues, eigenvectors)
    overlap = eigenvectors.transpose() @ reference_eigenvectors
    partners = []
    for i in range(overlap.shape[0]):
        partner_overlap = 0
        for j in range(overlap.shape[1]):
            if j not in partners and abs(overlap[i, j]) > abs(partner_overlap):
                partner = j
                partner_overlap = overlap[i, j]
        partners.append(partner)
    new_eigenvalues = np.array(eigenvalues)
    new_eigenvectors = np.array(eigenvectors)
    for i in range(overlap.shape[0]):
        new_eigenvalues[partners[i]] = eigenvalues[i]
        new_eigenvectors[:, partners[i]] = eigenvectors[:, i] if overlap[i, partners[i]] > 0 else -eigenvectors[:, i]
    return (new_eigenvalues, new_eigenvectors)
