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
    '''Class governing FARGO output

    Args:
        output_dir: Output directory
        dt: Time between fine-grain outputs
        Ninterm: Number of fine-grain outputs per full output
        Ntot: Total number of full outputs
    '''
    def __init__(self, output_dir, dt, Ninterm, Ntot):
        self.output_dir = output_dir
        self.dt = dt
        self.Ninterm = Ninterm
        self.Ntot = Ntot

class ShearingBox:
    '''Class governing shearing box parameters

    Args:
        dims: Dimensions of the box, should have 3 elements.
        mesh_size: grid points for each dimension, should have 3 elements.
    '''
    def __init__(self, dims, mesh_size):
        self.dims = dims
        self.mesh_size = mesh_size

class FargoSetup:
    '''Class for creating a new PSI FARGO setup

    Args:
        setup_name: Name of setup to be created
        fargo_dir (optional): Path to public FARGO directory when not using vanila submodule.
    '''
    def __init__(self, setup_name, fargo_dir=None):
        self.setup_name = setup_name

        # Path of current file, should be /path/to/psi_fargo/psi
        self.psi_dir = os.path.dirname(os.path.abspath(__file__))

        self.fargo_dir = fargo_dir

    def create(self, polydust, mode, shearing_box, output, cfl=0.44):
        '''Create FARGO setup

        Args:
            polydust: PolyDust object
            mode: SingleMode object, or None.
            shearing_box: ShearingBox object
            output: Output object
            cfl (optional): Courant number to use. Defaults to 0.44.
        '''
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

        perturbation = None
        if mode is not None:
            perturbation=mode.to_string()

        initial_conditions.write_condinit_file(polydust,
                                               perturbation=perturbation)
        created_files.extend(fargo_boundary.write_boundary_files(self.setup_name, polydust.N))

        # Setup directory
        output_dir = self.psi_dir + '/../public/setups/' + self.setup_name
        # If not using submodule FARGO
        if self.fargo_dir is not None:
            output_dir = self.fargo_dir + '/setups' + self.setup_name

        # Create setup directory if not exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            print("Warning: setup directory exists; contents will be overwritten")

        # Move all created files to setup directory
        for f in created_files:
            os.replace(f, output_dir + '/' + f)

        shutil.copy(self.psi_dir + '/boundaries.txt', output_dir + '/' + 'boundaries.txt')

        print('Setup created: ' + output_dir)
