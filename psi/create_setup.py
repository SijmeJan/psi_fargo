import os
import shutil
import numpy as np

import single_mode
import initial_conditions
import fargo_opt
import fargo_par
import fargo_boundary
from polydust import Polydust, SizeDistribution

class ShearingBox:
    def __init__(self, dims=None, mesh_size=None):
        self.dims = dims
        self.mesh_size = mesh_size

class FargoSetup:
    def __init__(self, setup_name, fargo_dir=None):
        self.setup_name = setup_name

        self.fargo_dir = fargo_dir

        # Installation directory not supplied: guess
        if fargo_dir is None:
            guess_dir = ['/Users/sjp/Codes/psi_fargo',
                         '/astro/sjp/Codes/psi_fargo']

            for direc in guess_dir:
                if os.path.isdir(direc):
                    self.fargo_dir = direc

        if self.fargo_dir is None:
            raise RuntimeError('Could not figure out FARGO installation directory!')
        print('Using Fargo directory ' + self.fargo_dir)

    def create(self, polydust, mode, shearing_box):
        created_files = [self.setup_name + '.opt',
                         self.setup_name + '.par',
                         'condinit.c']
        fargo_opt.write_opt_file(self.setup_name, polydust)
        fargo_par.write_par_file(self.setup_name, polydust, shearing_box)
        initial_conditions.write_condinit_file(polydust,
                                               perturbation=mode.to_string())
        created_files.extend(fargo_boundary.write_boundary_files(self.setup_name, polydust.N))

        # Create setup directory if not exists
        output_dir = self.fargo_dir + '/public/setups/' + self.setup_name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            print("Warning: setup directory exists; contents will be overwritten")

        # Move all created files to setup directory
        for f in created_files:
            os.replace(f, output_dir + '/' + f)

        shutil.copy('boundaries.txt', output_dir + '/' + 'boundaries.txt')

        print('Setup created: ' + output_dir)

# Number of dust fluids/ dust nodes
n_dust = 8

xi1 = np.log(0.001)
xi2 = np.log(0.1)
xi_edge = np.linspace(xi1, xi2, n_dust+1)

Ts_edge = np.exp(xi_edge)
Ts = xi_edge + 0.5*(xi_edge[1] - xi_edge[0])
Ts = np.exp(Ts)[0:-1]

q = np.power(Ts_edge, 0.5)

eps = 2*(np.roll(q, -1) - q)/(q[-1] - q[0])
eps = eps[0:-1]

# For continuum: min and max Stokes numbers
# For discrete multifluid: all Stokes numbers
#stokes_range = Ts
#stokes_range = [0.00316228, 0.03162278]  #[0.0425, 0.1]
stokes_range = [1.0e-3, 0.1]

# For continuum: SizeDistribution object
# For discrete multifluid: all dust densities
#size_distribution = eps
#size_distribution = [0.48050615, 1.51949385]
#size_distribution = [0.8403834887142569, 4.2359685932218305]
size_distribution = SizeDistribution(stokes_range)
dust_density = 2.0  #np.sum(size_distribution)
gas_density = 1.0

pd = Polydust(n_dust, stokes_range, dust_density,
              gas_density, size_distribution, gauss_legendre=False)

# Add single mode perturbation
Kx = 60
Kz = 60
#mode = single_mode.GasEpicycle(30, 30, 20, 0.01)
#mode = single_mode.Linear3(1.0e-5)
#mode = single_mode.RandomFixedK(n_dust, 1.0e-5, Kx, Kz)
mode = single_mode.PSI_pert(pd, 1.0e-5, Kx, Kz)

Ly = 2*np.pi/Kx            # 'radial' box size
Lz = 2*np.pi/Kz            # vertical box size
Ny = 32                          # 'radial' number of grid points
Nz = 32                          # vertical number of grid points

shearing_box = ShearingBox(dims=[0, Ly, Lz], mesh_size=[1, Ny, Nz])

try:
    setup = FargoSetup('psi_linearA')
except:
    raise

setup.create(pd, mode, shearing_box)
