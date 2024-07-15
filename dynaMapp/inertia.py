import numpy as np
import random

class Inertia:
    
    def __init__(self, XXi, XYi, XZi, YYi, YZi, ZZi) -> None:
        self.XX = XXi
        self.XY = XYi
        self.XZ = XZi
        self.YY = YYi
        self.YZ = YZi
        self.ZZ = ZZi

    def __str__(self):
        return (f"Inertia:\n"
                f"XX: {self.XX}, XY: {self.XY}, XZ: {self.XZ}\n"
                f"YY: {self.YY}, YZ: {self.YZ}, ZZ: {self.ZZ}\n")

    def _str_items(self, items):
        row_format = '\t{0:^6} : {1}\n'
        str_format = ""
        for idx, item in enumerate(items):
            str_format += row_format.format(str(idx), str(item))
        return str_format

    def init_random(self):
        self.XX = random.random()
        self.XY = random.random()
        self.XZ = random.random()
        self.YY = random.random()
        self.YZ = random.random()
        self.ZZ = random.random()

    @property
    def get_matrix(self):
        return np.array([
            [self.XX, self.XY, self.XZ],
            [self.XY, self.YY, self.YZ],
            [self.XZ, self.YZ, self.ZZ]
        ])

    @property
    def consistance(self):
        I = self.get_matrix
        return np.array_equal(I, I.T) and (
            self.XX + self.YY > self.ZZ and
            self.XX + self.ZZ > self.YY and
            self.YY + self.ZZ > self.XX
        )

    def rescale(self):
        if self.consistance:
            return self.get_matrix
        else:
            I = self.get_matrix
            eigvals, eigvecs = np.linalg.eigh(I)
            I = np.zeros_like(I)
            eigvals = np.maximum(np.abs(eigvals), 1e-2)
            for i in range(len(eigvals)):
                I += eigvals[i] * np.outer(eigvecs[:, i], eigvecs[:, i])
            return I

    def dh_to_com(self, m, cx, cy, cz):
        """ Convert inertia matrix from joint dh frame to com based frame
        using the parallel axis formula
        """
        I = self.get_matrix
        I[0, 0] += m * (cy**2 + cz**2)
        I[1, 1] += m * (cx**2 + cz**2)
        I[2, 2] += m * (cx**2 + cy**2)
        I[0, 1] += m * cx * cy
        I[0, 2] += m * cx * cz
        I[1, 2] += m * cy * cz
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        return I

    def com_to_dh(self, m, cx, cy, cz):
        """ Convert inertia matrix from joint com frame to dh based frame
        using the parallel axis formula
        """
        I = self.get_matrix
        I[0, 0] -= m * (cy**2 + cz**2)
        I[1, 1] -= m * (cx**2 + cz**2)
        I[2, 2] -= m * (cx**2 + cy**2)
        I[0, 1] -= m * cx * cy
        I[0, 2] -= m * cx * cz
        I[1, 2] -= m * cy * cz
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        return I

    def __eq__(self, other):
        """ Check if """
        return (self.XX == other.XX and self.XY == other.XY and self.XZ == other.XZ and
                self.YY == other.YY and self.YZ == other.YZ and self.ZZ == other.ZZ)

    def __add__(self, other):
        """ Add 2 instances of an inertia matrices"""
        return Inertia(self.XX + other.XX, self.XY + other.XY, self.XZ + other.XZ,
                       self.YY + other.YY, self.YZ + other.YZ, self.ZZ + other.ZZ)

    def matrix_to_vector(self):
        return np.array([self.XX, self.XY, self.XZ, self.YY, self.YZ, self.ZZ])
    
    def principal_moments(self):
        """Compute and return the principal moments of inertia."""
        I = self.get_matrix
        eigvals, _ = np.linalg.eigh(I)
        return eigvals
    
    def rotate(self, rotation_matrix):
        """Rotate the inertia tensor by a given rotation matrix."""
        I = self.get_matrix
        R = np.array(rotation_matrix)
        return R @ I @ R.T
    
    def axis_inertia_moment(self, axis):
        """Calculate the moment of inertia for a given axis."""
        I = self.get_matrix
        axis = np.array(axis)
        return axis.T @ I @ axis