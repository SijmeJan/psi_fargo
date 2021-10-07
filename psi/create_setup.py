import numpy as np

from fargo_setup import ShearingBox, Output, FargoSetup
import single_mode
from polydust import Polydust, SizeDistribution

# For continuum: min and max Stokes numbers
# For discrete multifluid: list all Stokes numbers
#stokes_range = [0.1]
stokes_range = [1.0e-2, 0.1]

pd = Polydust(n_dust = 4,
              stokes_range = stokes_range,
              dust_density = 3.0,
              gas_density = 1.0,
              #size_distribution = [1.0, 0.5],
              size_distribution = SizeDistribution(stokes_range),
              gauss_legendre=True,
              discrete_equilibrium=False)

# Add single mode perturbation
Kx = 30
Kz = 30
#mode = None
#mode = single_mode.GasEpicycle(30, 30, 20, 0.01)
#mode = single_mode.Linear3(1.0e-5)
#mode = single_mode.RandomFixedK(pd.N, 1.0e-5, Kx, Kz)
mode = single_mode.PSI_pert(pd, 1.0e-5, Kx, Kz)

#print('Discrete growth rate:', pd.growth_rate(Kx, Kz))

Ly = 2*np.pi/Kx                  # 'radial' box size
Lz = 2*np.pi/Kz                  # vertical box size
Ny = 32                          # 'radial' number of grid points
Nz = 32                          # vertical number of grid points

shearing_box = ShearingBox(dims=[0, Ly, Lz], mesh_size=[1, Ny, Nz])

output = Output('/Users/sjp/Codes/psi_fargo/data/mu1_K100/test',
                dt=0.01, Ninterm=10, Ntot=5000)

try:
    setup = FargoSetup('psi')
except:
    raise

setup.create(pd, mode, shearing_box, output, cfl=0.44)
