import shutil
import numpy as np
import os

import single_mode
import initial_conditions
import fargo_opt
import fargo_par
import fargo_boundary
from polydust import Polydust, SizeDistribution

class Output:
    def __init__(self, output_dir, dt, Ninterm, Ntot):
        self.output_dir = output_dir
        self.dt = dt
        self.Ninterm = Ninterm
        self.Ntot = Ntot

class ShearingBox:
    def __init__(self, dims=None, mesh_size=None):
        self.dims = dims
        self.mesh_size = mesh_size

class FargoSetup:
    def __init__(self, setup_name, fargo_dir=None):
        self.setup_name = setup_name

        # Path of current file, should be /path/to/psi_fargo/psi
        self.psi_dir = os.path.dirname(os.path.abspath(__file__))

    def create(self, polydust, mode, shearing_box, output, cfl=0.44):
        # Going to create a .opt file, a .par file and initial conditions
        created_files = [self.setup_name + '.opt',
                         self.setup_name + '.par',
                         'condinit.c']
        fargo_opt.write_opt_file(self.setup_name, polydust)
        fargo_par.write_par_file(self.setup_name,
                                 polydust,
                                 shearing_box,
                                 output,
                                 cfl=cfl)
        initial_conditions.write_condinit_file(polydust,
                                               perturbation=mode.to_string())
        created_files.extend(fargo_boundary.write_boundary_files(self.setup_name, polydust.N))

        # Create setup directory if not exists
        output_dir = self.psi_dir + '/../public/setups/' + self.setup_name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            print("Warning: setup directory exists; contents will be overwritten")

        # Move all created files to setup directory
        for f in created_files:
            os.replace(f, output_dir + '/' + f)

        shutil.copy(self.psi_dir + '/boundaries.txt', output_dir + '/' + 'boundaries.txt')

        print('Setup created: ' + output_dir)
