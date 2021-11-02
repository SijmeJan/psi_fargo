import sys
import numpy as np

from fargo_setup import ShearingBox, Output, FargoSetup
from polydust import Polydust

n_arg = len(sys.argv)
if n_arg != 2:
    print('Usage: python create_setup_from_output.py output_direc')
    exit()

output_direc = sys.argv[1]
print('Recreating setup from directory ' + output_direc)

# Get number of dust species from directory
pd = Polydust.from_output_direc(output_direc)

mode = None

# Get shearing box domain setup from domain_xyz.dat
shearing_box = ShearingBox.from_output_direc(output_direc)

output = Output(output_direc, dt=0.1, Ninterm=10, Ntot=10000)

try:
    setup = FargoSetup('psi')
except:
    raise

setup.create(pd, mode, shearing_box, output)
