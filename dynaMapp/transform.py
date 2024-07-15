import numpy as np
from urdfpy import matrix_to_rpy, rpy_to_matrix

class Transformation:
    
    def __init__(self, alpha, theta, d, a) -> None:
        self.alpha = alpha
        self.theta = theta
        self.d = d
        self.a = a

    def __str__(self):
        return (f"Transformation:\n"
                f"alpha: {self.alpha}, theta: {self.theta}, d: {self.d}, a: {self.a}\n")

    def _init_default(self) -> None:
        self.alpha = 0.0
        self.theta = 0.0
        self.a = 0.0
        self.d = 0.0

    def matrix(self):
        H = np.array([
            [np.cos(self.theta), -np.sin(self.theta) * np.cos(self.alpha), np.sin(self.theta) * np.sin(self.alpha), self.a * np.cos(self.theta)],
            [np.sin(self.theta), np.cos(self.theta) * np.cos(self.alpha), -np.cos(self.theta) * np.sin(self.alpha), self.a * np.sin(self.theta)],
            [0, np.sin(self.alpha), np.cos(self.alpha), self.d],
            [0, 0, 0, 1]
        ])
        return H

    def rotation(self):
        R = np.array([
            [np.cos(self.theta), -np.sin(self.theta) * np.cos(self.alpha), np.sin(self.theta) * np.sin(self.alpha)],
            [np.sin(self.theta), np.cos(self.theta) * np.cos(self.alpha), -np.cos(self.theta) * np.sin(self.alpha)],
            [0, np.sin(self.alpha), np.cos(self.alpha)]
        ])
        return R

    def translation(self):
        H = self.matrix()
        return H[:3, 3]

    def inv(self):
        H = self.matrix()
        return np.linalg.inv(H)

    def inv_rot(self):
        R = self.rotation()
        return np.linalg.inv(R)

    def rpy(self, solution=1):
        R = self.rotation()
        return matrix_to_rpy(R, solution)

    def to_xyz_rpy(self):
        H = self.matrix()
        xyz = H[:3, 3]
        rpy = matrix_to_rpy(H[:3, :3])
        return np.hstack((xyz, rpy))

    def compose(self, other):
        H1 = self.matrix()
        H2 = other.matrix()
        H_combined = H1 @ H2
        return Transformation.from_matrix(H_combined)

    def is_valid(self):
        H = self.matrix()
        return np.allclose(H @ np.linalg.inv(H), np.eye(4))

    def scale(self, factor):
        return Transformation(self.alpha * factor, self.theta * factor, self.d * factor, self.a * factor)

    def to_dict(self):
        return {
            "alpha": self.alpha,
            "theta": self.theta,
            "d": self.d,
            "a": self.a
        }

    @classmethod
    def from_dict(cls, data):
        return cls(data["alpha"], data["theta"], data["d"], data["a"])

    @classmethod
    def from_matrix(cls, H):
        a = H[0, 3]
        d = H[2, 3]
        theta = np.arctan2(H[1, 0], H[0, 0])
        alpha = np.arctan2(H[2, 1], H[2, 2])
        return cls(alpha, theta, d, a)

    def to_matrix(self):
        return self.matrix()

    def to_rotation_matrix(self):
        return self.rotation()


s = Transformation(0.1,0,1,1)
p= s.rpy(1.0)
print(p)