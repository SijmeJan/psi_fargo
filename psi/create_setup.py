import os
import shutil
import numpy as np

def str_perturbation(amp, Kx, Kz, eig, stagger=None):
    k = '{}*Ymed(j) + {}*Zmed(k)'.format(Kx, Kz)
    # Stagger in y
    if stagger == 'y':
        k = '{}*Ymin(j) + {}*Zmed(k)'.format(Kx, Kz)
    if stagger == 'z':
        k = '{}*Ymed(j) + {}*Zmin(k)'.format(Kx, Kz)

    return '{}*({}*cos('.format(amp, np.real(eig)) + k + ') - {}*sin('.format(np.imag(eig)) + k + '))'

def read_file(file_name):
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()

    return lines

def write_file(file_name, lines):
    f = open(file_name, "w")
    f.writelines(lines)
    f.close()

def add_line(lines, str_line):
    '''Add line to end of list'''
    n = len(lines)
    lines[n:n] = [str_line]

def add_option(lines, option):
    '''Add option to opt file content'''
    add_line(lines, 'FARGO_OPT +=  -D'+ option + '\n')

def write_opt_file(setup_name, n_dust):
    lines = ['# PSI FARGO2D setup using ' + str(n_dust) + ' dust fluids\n\n']

    line = 'FLUIDS := 0'
    for i in range(1, n_dust + 1):
        line = line + ' ' + str(i)
    add_line(lines, line +'\n')
    add_line(lines, 'NFLUIDS = ' + str(n_dust + 1) + '\n')
    add_line(lines, 'FARGO_OPT += -DNFLUIDS=${NFLUIDS}\n\n')

    add_option(lines, 'X')
    add_option(lines, 'Y')
    add_option(lines, 'Z')

    add_option(lines, 'CARTESIAN')
    add_option(lines, 'SHEARINGBOX')
    add_option(lines, 'SHEARINGBC')
    add_option(lines, 'ISOTHERMAL')

    add_option(lines, 'DRAGFORCE')
    add_option(lines, 'CONSTANTSTOKESNUMBER')

    add_line(lines, 'MONITOR_SCALAR = MOM_X | MOM_Y | MOM_Z\n')

    add_line(lines, '\n#Cuda blocks\n')
    add_line(lines, 'ifeq (${GPU}, 1)\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_X=16\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_Y=16\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_Z=1\n')
    add_line(lines, 'endif\n')

    fname = setup_name + '.opt'
    write_file(fname, lines)

    return fname

def write_par_file(setup_name, stokes, dust_to_gas_ratio, box_size, grid):
    n_dust = len(stokes)

    lines = ['# PSI FARGO2D setup using ' + str(n_dust) + ' dust fluids\n\n']

    add_line(lines, 'Setup          ' + setup_name + '\n')

    add_line(lines, '\n### Disk parameters\n\n')
    add_line(lines, 'AspectRatio        1.0            Make eta source term 2\n')
    add_line(lines, 'Nu                 0.0            Uniform kinematic viscosity\n')
    add_line(lines, 'OmegaFrame         1.0      Rotation rate\n')
    add_line(lines, 'ShearParam         1.5      Shear rate\n')

    add_line(lines, '\n### Dust parameters\n\n')

    for i in range(1, n_dust+1):
        add_line(lines, 'Invstokes' + str(i) + '        ' + str(1/stokes[i-1]) + '   Inverse of the Stokes number for dust' + str(i) + '\n')

    add_line(lines, '\nEpsilon                 ' + str(dust_to_gas_ratio) + '    Dust-to-gas mass ratio\n')

    add_line(lines, '\n### Mesh parameters\n\n')

    add_line(lines, 'Ny         {}      Number of \'radial\' zones\n'.format(grid[0]))
    add_line(lines, 'Nz         {}      Number of vertical zones\n'.format(grid[1]))

    add_line(lines, 'Ymin              ' + str(-0.5*box_size[0]) + '\n')
    add_line(lines, 'Ymax              ' + str(0.5*box_size[0]) + '\n')
    add_line(lines, 'Zmin              ' + str(-0.5*box_size[1]) + '\n')
    add_line(lines, 'Zmax              ' + str(0.5*box_size[1]) + '\n')

    add_line(lines, 'PeriodicY          YES\n')
    add_line(lines, 'PeriodicZ          YES\n')

    add_line(lines, '\n### Output control parameters\n\n')

    add_line(lines, 'DT         0.1     Time step length\n')
    add_line(lines, 'Ninterm    10      Time steps between outputs\n')
    add_line(lines, 'Ntot       2000    Total number of time steps\n')

    add_line(lines, 'OutputDir      @outputs/' + setup_name + '\n')

    fname = setup_name + '.par'
    write_file(fname, lines)

    return fname

def write_condinit_file(stokes, sigma):
    AN = np.sum(sigma*stokes/(1 + stokes*stokes))
    BN = 1.0 + np.sum(sigma/(1 + stokes*stokes))

    vgx = 2*AN/(AN*AN + BN*BN)
    vgy = -BN/(AN*AN + BN*BN)     # Excluding background shear
    vgz = 0.0

    vdx = (vgx + 2*stokes*vgy)/(1 + stokes*stokes)
    vdy = (vgy - 0.5*stokes*vgx)/(1 + stokes*stokes)
    vdz = 0.0*stokes

    # Eigenvector LinA
    drhog = 0.0000074637 + 1j*0.0000070677
    dvgx = -0.0563787907 + 1j*0.0120535455
    dvgy =  0.0445570113 + 1j*0.0197224299
    dvgz =  0.0563784989 - 1j*0.0120536242
    drhod = [1.0]
    dvdx =  [-0.0466198076 + 1j*0.0124333223]
    dvdy =  [0.0435211557 + 1j*0.0213517453]
    dvdz =  [0.0546507401 - 1j*0.0077776652]
    Kx = 30.0
    Kz = 30.0
    amp = 1.0e-4


    # Convert to FARGO standard where y=x and x=y....
    vgx, vgy = vgy, vgx
    vdx, vdy = vdy, vdx
    dvgx, dvgy = dvgy, dvgx
    dvdx, dvdy = dvdy, dvdx

    lines = \
      ['#include "fargo3d.h"\n\n',
       'void CondInit() {\n',
       '  int id_gas = 0;\n',
       '  int feedback = YES;\n',
       '  //We first create the gaseous fluid and store it in the array Fluids\n',
       '  Fluids[id_gas] = CreateFluid("gas",GAS);\n\n',
       '  //We now select the fluid\n',
       '  SelectFluid(id_gas);\n',
       '  //and fill its fields\n',
       '  int i,j,k;\n',
       '  real *rho  = Density->field_cpu;\n',
       '  real *cs   = Energy->field_cpu;\n',
       '  real *vx   = Vx->field_cpu;\n',
       '  real *vy   = Vy->field_cpu;\n',
       '  real *vz   = Vz->field_cpu;\n',
       '  for (k=0; k<Nz+2*NGHZ; k++) {\n',
       '    for (j=0; j<Ny+2*NGHY; j++) {\n',
       '      for (i=0; i<Nx+2*NGHX; i++) {\n',
       '        rho[l] = 1.0 + 0.0*exp(-((Ymed(j)-0.05)*(Ymed(j)-0.05) + (Zmed(k)-0.05)*(Zmed(k)-0.05))/(0.025*0.025)) + '+ str_perturbation(amp, Kx, Kz, drhog) + ';\n',
       '        vx[l]  = {} + '.format(vgx) + str_perturbation(amp, Kx, Kz, dvgx) + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
       '        vy[l]  = {} + '.format(vgy) + str_perturbation(amp, Kx, Kz, dvgy, stagger='y') + ';\n',
       '        vz[l]  = {} + '.format(vgz) + str_perturbation(amp, Kx, Kz, dvgz, stagger='z') + ';\n',
       '        cs[l]  = 20.0;      // H*Omega*Omega/eta\n',
       '      }\n',
       '    }\n',
       '  }\n',
       '  \n  //We repeat the process for the dust fluids\n',
       '  char dust_name[MAXNAMELENGTH];\n']

    for n in range(0, len(stokes)):
        lines.extend(\
          ['  sprintf(dust_name,"dust%d",{}); //We assign different names to the dust fluids\n'.format(n+1),
           '  Fluids[{}]  = CreateFluid(dust_name, DUST);\n'.format(n+1),
           '  SelectFluid({});\n'.format(n+1),
           '  rho  = Density->field_cpu;\n',
           '  cs   = Energy->field_cpu;\n',
           '  vx   = Vx->field_cpu;\n',
           '  vy   = Vy->field_cpu;\n',
           '  vz   = Vz->field_cpu;\n',
           '  for (k=0; k<Nz+2*NGHZ; k++) {\n',
           '    for (j=0; j<Ny+2*NGHY; j++) {\n',
           '      for (i=0; i<Nx+2*NGHX; i++) {\n',
           '        rho[l] = {} + '.format(sigma[n]) + str_perturbation(amp, Kx, Kz, drhod[n]) + ';\n',
           '        vx[l]  = {} + '.format(vdx[n]) + str_perturbation(amp, Kx, Kz, dvdx[n]) + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
           '        vy[l]  = {} + '.format(vdy[n]) + str_perturbation(amp, Kx, Kz, dvdy[n], stagger='y') + ';\n',
           '        vz[l]  = {} + '.format(vdz[n]) + str_perturbation(amp, Kx, Kz, dvdz[n], stagger='z') + ';\n',
           '        cs[l]  = 0.0;\n',
           '      }\n',
           '    }\n',
           '  }\n\n'])

    lines.extend(\
      ['  /*We now fill the collision matrix (Feedback from dust included)\n',
       '  Note: ColRate() moves the collision matrix to the device.\n',
       '  If feedback=NO, gas does not feel the drag force.*/\n'])


    for i in range(1, len(stokes) + 1):
     add_line(lines, '  ColRate(INVSTOKES'+str(i)+', id_gas, '+str(i)+', feedback);\n')

    add_line(lines, '}\n')

    fname = 'condinit.c'
    write_file(fname, lines)

    return fname

def write_boundary_files(setup_name, n_dust):
    fnames = []

    for i in range(0, n_dust+1):
        lines = \
          ['#Boundaries file for fluid number '+str(i)+'\n',
           '#----------------------------------\n\n',
           'Density:\n',
           '    Ymin: NOBOUNDARY\n',
           '    Ymax: NOBOUNDARY\n',
           '    Zmin: NOBOUNDARY\n',
           '    Zmax: NOBOUNDARY\n',
           'Vx:\n',
           '    Ymin: NOBOUNDARY\n',
           '    Ymax: NOBOUNDARY\n',
           '    Zmin: NOBOUNDARY\n',
           '    Zmax: NOBOUNDARY\n',
           'Vy:\n',
           '    Ymin: NOBOUNDARY\n',
           '    Ymax: NOBOUNDARY\n',
           '    Zmin: NOBOUNDARY\n',
           '    Zmax: NOBOUNDARY\n',
           'Vz:\n',
           '    Ymin: NOBOUNDARY\n',
           '    Ymax: NOBOUNDARY\n',
           '    Zmin: NOBOUNDARY\n',
           '    Zmax: NOBOUNDARY\n']

        lines = ['']
        file_name = setup_name + '.bound.'+str(i)
        write_file(file_name, lines)

        fnames = fnames + [file_name]

    return fnames

fargo_dir = '/Users/sjp/Codes/fargo3d/public/setups'
setup_name = 'psi_linearA'      # Name of the PSI setup
stokes = np.asarray([0.1])      # List of Stokes numbers
dust_to_gas_ratio = 3.0         # Dust to gas ratio
Ly = 2*np.pi/30.0               # 'radial' box size
Lz = 2*np.pi/30.0               # vertical box size
Ny = 64                         # 'radial' number of grid points
Nz = 64                         # vertical number of grid points


n_dust = len(stokes)            # number of dust fluids
sigma = np.asarray([dust_to_gas_ratio])     # for single size!
created_files = []

created_files.append(write_opt_file(setup_name, n_dust))
created_files.append(write_par_file(setup_name, stokes, dust_to_gas_ratio, [Ly, Lz], [Ny, Nz]))
created_files.append(write_condinit_file(stokes, sigma))
created_files.extend(write_boundary_files(setup_name, n_dust))

output_dir = fargo_dir + '/' + setup_name
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for f in created_files:
    os.replace(f, output_dir + '/' + f)

shutil.copy('boundaries.txt', output_dir + '/' + 'boundaries.txt')

print(created_files)
