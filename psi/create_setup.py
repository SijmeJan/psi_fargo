import os
import shutil
import numpy as np

import single_mode
import initial_conditions
import fargo_opt
import fargo_par
import fargo_boundary
from polydust import Polydust, SizeDistribution

# Number of dust fluids/ dust nodes
n_dust = 8

# For continuum: min and max Stokes numbers
# For discrete multifluid: all Stokes numbers
#stokes_range = [0.0425, 0.1]
stokes_range = [1.0e-8, 1.0]

# For continuum: SizeDistribution object
# For discrete multifluid: all dust densities
#size_distribution = [1.0, 0.5]
size_distribution = SizeDistribution(stokes_range)

pd = Polydust(n_dust, stokes_range, size_distribution)
dust_total_density = 10.0

# Add single mode perturbation
Kx = 10
Kz = 10
#mode = single_mode.Linear3(1.0e-4)
mode = single_mode.RandomFixedK(n_dust, 1.0e-4, Kx, Kz)

fargo_dir = '/Users/sjp/Codes/psi_fargo/public/setups'
setup_name = 'psi_mu3'         # Name of the PSI setup
Ly = 2*np.pi/mode.Kx            # 'radial' box size
Lz = 2*np.pi/mode.Kz            # vertical box size
Ny = 32                         # 'radial' number of grid points
Nz = 32                         # vertical number of grid points


# Create all files needed for FARGO setup
created_files = []
created_files.append(fargo_opt.write_opt_file(setup_name, pd))
created_files.append(fargo_par.write_par_file(setup_name, [Ly, Lz], [Ny, Nz], pd))
created_files.append(initial_conditions.write_condinit_file(pd, dust_total_density, perturbation=mode.to_string()))
created_files.extend(fargo_boundary.write_boundary_files(setup_name, n_dust))

# Create setup directory if not exists
output_dir = fargo_dir + '/' + setup_name
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
    print("Warning: setup directory exists; contents will be overwritten")

# Move all created files to setup directory
for f in created_files:
    os.replace(f, output_dir + '/' + f)

shutil.copy('boundaries.txt', output_dir + '/' + 'boundaries.txt')

print('Setup created: ' + output_dir)
