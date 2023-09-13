import numpy as np
from scipy.integrate import quad
from scipy.constants import mu_0, pi

class BSIntegrator:
    
    def __init__(self, current_profile, r_range, z_range):
        
        assert len(current_profile.shape) == 2, "Current profile grid should be 2 dimensional"
        self.current_profile = current_profile
        
        # Get minimum and maximum coordinate values
        
        self.r_min = r_range[0]
        self.r_max = r_range[1]
        
        self.z_min = z_range[0]
        self.z_max = z_range[1]
        
        # Get r and z step sizes
        
        self.dr = (self.r_max - self.r_min) / current_profile.shape[1]
        self.dz = (self.z_max - self.z_min) / current_profile.shape[0]
        
        # Create coordinate grid
        
        self.r_grid, self.z_grid = np.meshgrid(np.linspace(self.r_min, self.r_max, current_profile.shape[1]), np.linspace(self.z_min, self.z_max, current_profile.shape[0]))
        
        # Center coordinates on the grid centers
        
        self.r_grid += self.dr / 2
        self.z_grid += self.dz / 2
        
    @staticmethod
    def diff_B_r(phi, r_prime, z_prime, r, z):
        
        # Calculate r component of magnetic field along a current loop
        
        d_squared = np.square(r) + np.square(r_prime) + np.square(z - z_prime) - (2 * r * r_prime * np.cos(phi))
        
        return np.cos(phi) / np.power(d_squared, 1.5)
    
    @staticmethod
    def diff_B_z(phi, r_prime, z_prime, r, z):
        
        # Calculate z component of magnetic field along a current loop
        
        d_squared = np.square(r) + np.square(r_prime) + np.square(z - z_prime) - (2 * r * r_prime * np.cos(phi))
        
        return (r_prime - r * np.cos(phi)) / np.power(d_squared, 1.5)
    
    def CurrentLoop_r(self, j, r_prime, z_prime, r, z):
 
        # Integrate diff_B_r to get the r component magnetic field created by one current loop
        
        # If the current is zero, the entire integral is zero
        if j == 0:
            return 0
        
        return j * r_prime * (z - z_prime) * quad(self.diff_B_r, -pi, pi, args = (r_prime, z_prime, r, z))[0]
    
    def CurrentLoop_z(self, j, r_prime, z_prime, r, z):
        
        # Integrate diff_B_z to get the z component magnetic field created by one current loop
        
        if j == 0:
            return 0
        
        return j * r_prime * quad(self.diff_B_z, -pi, pi, args = (r_prime, z_prime, r, z))[0]
    
    def get_B_r(self, r, z):
        
        # Calculate magnetic field contributions in each grid and sum them up
        
        # Vectorize CurrentLoop_r() to work with coordinate grids
        func = np.vectorize(self.CurrentLoop_r, excluded = ['r', 'z'])
        
        values = self.dr * self.dz * func(self.current_profile, self.r_grid, self.z_grid, r, z)
        
        return mu_0 / (4 * pi) * values.sum()
    
    def get_B_z(self, r, z):
        
        # Same for the z component of the field
        
        func = np.vectorize(self.CurrentLoop_z, excluded = ['r', 'z'])
        
        values = self.dr * self.dz * func(self.current_profile, self.r_grid, self.z_grid, r, z)
        
        return mu_0 / (4 * pi) * values.sum()
        