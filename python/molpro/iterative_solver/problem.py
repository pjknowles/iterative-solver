import numpy as np
class Problem:
    def __init__(self):
        pass

    def residual(self, parameters, residual):
        return None

    def action(self, parameters, action):
        return None

    def diagonals(self, diagonals):
        return False

    def precondition(self, residual, shift=None, diagonals=None):
        # print('precondition',residual.shape,residual,shift,diagonals)
        small = 1e-14
        if len(residual.shape) > 1:
            for i in range(residual.shape[1]):
                self.precondition(residual[i, :], shift[i] if shift is not None else None, diagonals)
            return
        if diagonals is not None:
            # nbuffer = residual.shape[1] if len(residual.shape) > 1 else 1
            if shift is not None:
                shift_ = shift if type(shift) is float else shift[0]
                for j in range(residual.size):
                    residual[j] = residual[j] / (diagonals[j] + shift_ + small)
            else:
                for j in range(residual.size):
                    residual[j] = residual[j] / (diagonals[j] + small)

    def pp_action_matrix(self, pparams):
        return np.array([], dtype=np.double)

    def p_action(self, p_coefficients, pparams, actions):
        raise NotImplementedError('P-space unavailable: unimplemented p_action() in Problem class')
        pass

    def test_parameters(self, instance, parameters):
        return False

    def report(self, iteration, verbosity, errors, value=None, eigenvalues=None):
        if (iteration <=0 and verbosity >=1) or verbosity>=2:
            if iteration > 0 and verbosity >= 2:
                print('Iteration', iteration,'log10(|residual|)=', np.log10(errors))
            elif iteration == 0:
                print('Converged', 'log10(|residual|)=', np.log10(errors))
            else:
                print('Unconverged', 'log10(|residual|)=', np.log10(errors))
            if value is not None:
                print('Objective function value', value)
            if eigenvalues is not None:
                print('Eigenvalues', eigenvalues)
            return True
        else:
            return False
