import numpy as np


class PSpace:
    '''
    Container for the definition of the P space in Iterative_Solver
    '''

    def __init__(self):
        self.size = 0
        self.offsets = [0]
        self.indices = []
        # self.coefficients = np.empty(shape=(0), dtype=float)
        self.coefficients = []
        self.simple = True
        pass

    def add_complex(self, indices, coefficients):
        assert(len(indices) == len(coefficients))
        # print('add_complex', indices, coefficients)
        self.simple = self.simple and len(indices) <= 1
        self.offsets.append(len(self.indices) + len(indices))
        self.indices = self.indices + indices
        self.coefficients = self.coefficients + coefficients
        # np.append(self.coefficients, coefficients)
        self.size = self.offsets[-1]
        # print(self.indices, self.coefficients,self.offsets)

    def add_simple(self, indices):
        # print ('PSpace.add_simple', indices)
        for i in indices:
            self.add_complex([i], [1.0])
