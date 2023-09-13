from biot_savart import BSIntegrator
import matplotlib.pyplot as plt
from scipy.constants import mu_0
import numpy as np

# Test the integrator against the analytical Biot-Savart solution on the z axis

def analytical(z):
    
    return mu_0 * current / 2 * 1 / np.power(1 + z ** 2, 1.5)

current_profile = np.zeros((101, 101), dtype = float)

dr = float(2 / 101)
dz = float(2 / 101)

current = 10000
current_profile[50, 50] = current / dr / dz

integrator = BSIntegrator(current_profile, (0, 2), (-1, 1))

z_range = np.linspace(-5, 5, 100)

B_z = [integrator.get_B_z(0, z) for z in z_range]
anal_B_z = [analytical(z) for z in z_range]

plt.plot(z_range, B_z, color = "blue")
plt.plot(z_range, anal_B_z, color = "red")
