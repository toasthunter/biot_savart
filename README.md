# What does it do?

Let's assume a cylindrical coordinate system $(r, \phi, z)$ is given. Let's take a current density $\mathbf j (r, z)$ in the $rz$ plane (i.e. the cross section of this coordinate system), and assume that the current is flowing in the $\phi$ direction(around the z axis circularly). This current density can be taken as rotationally symmetric, i.e. it doesn't depend on the $\phi$ coordinate. We want to solve for the magnetic field $\mathbf B(r, z)$. Thus we use the Biot-Savart law:

$$\mathbf B(\mathbf r) = \frac{\mu_0}{4 \pi}\iiint_V \frac{\mathbf j(\mathbf r') \times (\mathbf r -  \mathbf r')}{|\mathbf r - \mathbf r'|^3} d^3 r'$$

...after some parametrization, we get the following formulas for the $r$ and $z$ components of the magnetic field:

$$B_r(r, z) = \frac{\mu_0}{4 \pi} \int dr' \int d \phi' \int dz' \: r' \frac{j(r', z')(z - z') \cos(\phi')}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi))^{\frac{3}{2}}}$$
$$B_z(r, z) = \frac{\mu_0}{4 \pi} \int dr' \int d \phi' \int dz' \: r' \frac{j(r', z')(r' - r \cos(\phi'))}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi))^{\frac{3}{2}}}$$

And so, the program needs to solve for these integrals. 

# How does it work?

We discretize $j(r, z)$ into a 2-dimensional array as an input for our integrator. Each little "square" in the current grid forms a circular conductor in the coordinate system, and we use numpy's `quad` scheme to numerically integrate on each circular conductor:

$$dB_r = jr'(z - z') \int_{-\pi}^{\pi} d \phi' \frac{\cos(\phi')}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi'))^{\frac{3}{2}}} $$

$$dB_z = jr' \int_{-\pi}^{\pi} d \phi' \frac{r' - r \cos(\phi')}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi'))^{\frac{3}{2}}}$$

This is done by `CurrentLoop_r()` and `CurrentLoop_z()` respectively.

Second, we need to evaluate all of these loop integrals on the entire grid. We use a discrete sum over the current grid for that. So, the program evaluates the following:

$$B_r(r, z) = \frac{\mu_0}{4 \pi} \sum_{r', z'} \Delta r' \Delta z'\cdot jr'(z - z') \int_{-\pi}^{\pi} d \phi' \frac{\cos(\phi')}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi'))^{\frac{3}{2}}} $$
$$B_z(r, z) =  \frac{\mu_0}{4 \pi} \sum_{r', z'} \Delta r' \Delta z'\cdot jr' \int_{-\pi}^{\pi} d \phi' \frac{r' - r \cos(\phi')}{(r^2 + r'^2 + (z - z')^2 - 2 r r' \cos(\phi'))^{\frac{3}{2}}}$$

The integral is done by scipy's quad function, and the discrete sum is performed with numpy arrays. `get_B_r()` and `get_B_z()` perform these calculations and return the final value for the magnetic field components.

