import numpy as np

class Friction:
    
    def __init__(self, coulomb, stribeck, velocity, tspan, tstep=0.01) -> None:
        """
        Initialize the friction model with given parameters.

        Args:
            coulomb (float): Coulomb friction coefficient.
            stribeck (float): Stribeck friction coefficient.
            velocity (float): Joint velocity.
            tspan (float): Final simulation time.
            tstep (float): Step time for simulation. Default is 0.01.
        """
        self.coulomb = coulomb
        self.stribeck = stribeck
        self.velocity = velocity
        self.vs = 0.01
        self.tstep = tstep
        self.tspan = tspan
        self.sigma0 = 0.12
        self.sigma1 = 0.25
        self.sigma2 = 0.21
        self.z = 0.0
        self.tinit = 0.0

    @property
    def coulomb_friction(self):
        """Calculate Coulomb friction."""
        return self.coulomb * np.sign(self.velocity)
    
    @property
    def stribeck_friction(self):
        """Calculate Stribeck friction."""
        return (self.coulomb + (self.stribeck - self.coulomb) * 
                np.exp(-(self.velocity / self.vs) ** 2))
    
    def initialize_state(self, z0=0.0, tinit=0.0):
        """Initialize the internal state variables."""
        self.z = z0
        self.tinit = tinit

    def update_coefficients(self, sigma0, sigma1, sigma2):
        """Update the model's coefficients."""
        self.sigma0 = sigma0
        self.sigma1 = sigma1
        self.sigma2 = sigma2

    def lugre(self):
        """Compute the LuGre friction model."""
        gv = (self.coulomb + (self.stribeck - self.coulomb) * 
              np.exp(-(self.velocity / self.vs) ** 2)) / self.sigma0
        z_dot = self.velocity - abs(self.velocity) * self.z / gv
        self.z += z_dot * self.tstep
        F = self.sigma0 * self.z + self.sigma1 * z_dot + self.sigma2 * self.velocity
        return F

    def lugre_friction(self):
        """
        Compute LuGre Friction Model over a time span.

        Returns:
            tuple: (t, F, Fss, error) where
                t (numpy.ndarray): Simulation time vector.
                F (numpy.ndarray): Friction force.
                Fss (float): Steady-state friction force.
                error (numpy.ndarray): Difference between friction force and steady-state force.
        """
        self.initialize_state()
        t = np.arange(self.tinit, self.tspan + self.tstep, self.tstep)
        Fss = self.coulomb * np.sign(self.velocity) + (self.stribeck - self.coulomb) * \
              np.exp(-(self.velocity / self.vs) ** 2) * np.sign(self.velocity) + self.sigma2 * self.velocity
        F = np.empty_like(t)
        error = np.empty_like(t)

        for j in range(len(t)):
            F[j] = self.lugre()
            error[j] = abs(F[j] - Fss)

        return t, F, Fss, error

    def hysteresis_force(self, Fb, W, X, K):
        """
        Compute the hysteresis force for the Leuven model.

        Args:
            Fb (float): Base hysteresis force.
            W (list): List of weights.
            X (list): List of displacement values.
            K (list): List of stiffness coefficients.

        Returns:
            float: The computed hysteresis force.
        """
        assert len(W) == len(X) == len(K), "All input lists must have the same length"
        N = len(W)
        Fh = 0
        for i in range(N):
            if abs(self.z - X[i]) < W[i] / K[i]:
                Fh += K[i] * (self.z - X[i])
            else:
                Fh += np.sign(self.z - X[i]) * W[i]
        Fh += Fb
        return Fh

    def leuven(self, Fb, W, X, K, n):
        """
        Compute the Leuven friction model.

        Args:
            Fb (float): Base hysteresis force.
            W (list): List of weights.
            X (list): List of displacement values.
            K (list): List of stiffness coefficients.
            n (float): Transition curve shape coefficient.

        Returns:
            float: The computed Leuven friction force.
        """
        Fh = self.hysteresis_force(Fb, W, X, K)
        sv = self.stribeck_friction
        z_dot = self.velocity * ((1 - np.sign(Fh / sv)) * abs(Fh / sv) ** n)
        self.z += z_dot * self.tstep
        F = Fh + self.sigma1 * z_dot + self.sigma2 * self.velocity
        return F

    def leuven_friction(self, Fb, W, X, K, n):
        """
        Compute the Leuven Friction Model over a time span.

        Args:
            Fb (float): Base hysteresis force.
            W (list): List of weights.
            X (list): List of displacement values.
            K (list): List of stiffness coefficients.
            n (float): Transition curve shape coefficient.

        Returns:
            tuple: (t, F, Fss, error) where
                t (numpy.ndarray): Simulation time vector.
                F (numpy.ndarray): Friction force.
                Fss (float): Steady-state friction force.
                error (numpy.ndarray): Difference between friction force and steady-state force.
        """
        self.initialize_state()
        t = np.arange(self.tinit, self.tspan + self.tstep, self.tstep)
        F = np.empty_like(t)
        Fss = self.coulomb * np.sign(self.velocity) + (self.stribeck - self.coulomb) * \
              np.exp(-(self.velocity / self.vs) ** 2) * np.sign(self.velocity) + self.sigma2 * self.velocity
        error = np.empty_like(t)
        for j in range(len(t)):
            F[j] = self.leuven(Fb, W, X, K, n)
            error[j] = abs(F[j] - Fss)

        return t, F, Fss, error

    def dahl_friction(self, K, disp):
        """
        Compute the Dahl Friction Model.

        Args:
            K (float): Stiffness coefficient.
            disp (list): List of incremental displacements.

        Returns:
            list: List of computed friction forces.
        """
        F = 0     
        self.initialize_state()
        friction_force = []
        for dx in disp:
            if dx != 0:
                dz = dx - (dx / abs(dx)) * (F / K)
            else:
                dz = 0
            self.z += dz
            F = self.coulomb * (1 - (2 * abs(self.z)) / (abs(self.z) + 1))
            friction_force.append(F)

        return friction_force
