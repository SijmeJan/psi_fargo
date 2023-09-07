import numpy as np

from fargo_setup import ShearingBox, Output, FargoSetup
import single_mode
from polydust import Polydust, SizeDistribution

# Name of the setup. Any existing setup with the same name will be
# overwritten!
setup_name = 'psi'

# FARGO directory: /path/to/fargo/public
# If None, it is assumed that FARGO was installed as a submodule.
fargo_dir = None

##############################################################################
# Define how output: how much and when
#
# Follow the FARGO terminology: dt is the amount of time between lightweight
# output, and Ninterm is the amount of lightweight outputs between full data
# dumps. Therefore, the code will perform a dump every dt*Ninterm. The total
# simulation time is dt*Ntot.
##############################################################################

# Data output directory, to be put in FARGO par file
output_direc = './'

# Output object: how many dumps, how much total time
output = Output(output_direc, dt=0.1, Ninterm=10, Ntot=10000)

##############################################################################
# Gas and dust properties
##############################################################################

# Gas viscosity
viscous_alpha = 0.0

# Gas density (constant, since not stratified)
gas_density = 1.0

# Total dust density (constant, since not stratified)
dust_density = 3.0

##############################################################################
# Define the Stokes number distribution
#
# We can either have a discrete number of Stokes numbers, or a continuous
# distribution between a minimum and a maximum Stokes number.
##############################################################################

# Number of dust species. In the case of a continuous size distribution, this
# is the number of Gauss-Legendre nodes
n_dust = 1

# For discrete multifluid: list all Stokes numbers
stokes_range = [0.1]
# For continuum: min and max Stokes numbers
#stokes_range = [1.0e-2, 0.1]

# Discrete sizes: list of dust densities
size_dist = [3.0]
# Continuous size distribution
#size_dist = SizeDistribution(stokes_range),

# Create polydust object
pd = Polydust(n_dust=n_dust,
              stokes_range=stokes_range,
              dust_density=dust_density,
              gas_density=gas_density,
              size_distribution=size_dist,
              gauss_legendre=True,
              discrete_equilibrium=False)

# Add single mode perturbation
Kx = 30
Kz = 30
#mode = None
mode = single_mode.LinearA(1.0e-5)
#mode = single_mode.GasEpicycle(30, 30, 20, 0.01)
#mode = single_mode.Linear3(1.0e-5)
#mode = single_mode.RandomFixedK(pd.N, 1.0e-5, Kx, Kz)
#mode = single_mode.PSI_pert(pd, 1.0e-5, Kx, Kz,
#                            viscous_alpha=viscous_alpha)

##############################################################################
# Define the shearing box
##############################################################################

Ly = 2*np.pi/Kx                  # 'radial' box size
Lz = 2*np.pi/Kz                  # vertical box size
Ny = 16                          # 'radial' number of grid points
Nz = 16                          # vertical number of grid points

shearing_box = ShearingBox(dims=[0, Ly, Lz], mesh_size=[1, Ny, Nz])

##############################################################################
# Write the FARGO setup files
##############################################################################

try:
    setup = FargoSetup(setup_name, fargo_dir=fargo_dir)
except:
    raise

setup.create(pd, mode, shearing_box, output,
             cfl=0.44, viscous_alpha=viscous_alpha)
