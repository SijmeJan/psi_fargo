from common import *

def write_par_file(setup_name, pd, shearing_box, output, cfl=0.44):
    n_dust = pd.N
    stokes = pd.dust_nodes()

    box_size = shearing_box.dims
    grid = shearing_box.mesh_size

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

    #add_line(lines, '\nEpsilon                 ' + str(dust_to_gas_ratio) + '    Dust-to-gas mass ratio\n')

    add_line(lines, '\n### Mesh parameters\n\n')

    add_line(lines, 'Ny         {}      Number of \'radial\' zones\n'.format(grid[1]))
    add_line(lines, 'Nz         {}      Number of vertical zones\n'.format(grid[2]))

    add_line(lines, 'Ymin              ' + str(-0.5*box_size[1]) + '\n')
    add_line(lines, 'Ymax              ' + str(0.5*box_size[1]) + '\n')
    add_line(lines, 'Zmin              ' + str(-0.5*box_size[2]) + '\n')
    add_line(lines, 'Zmax              ' + str(0.5*box_size[2]) + '\n')

    add_line(lines, 'PeriodicY          YES\n')
    add_line(lines, 'PeriodicZ          YES\n')

    add_line(lines, '\n### Output control parameters\n\n')

    add_line(lines, 'CFL        {}     Courant number\n'.format(cfl))
    add_line(lines, 'DT         {}     Time step length\n'.format(output.dt))
    add_line(lines, 'Ninterm    {}     Time steps between outputs\n'.format(output.Ninterm))
    add_line(lines, 'Ntot       {}    Total number of time steps\n'.format(output.Ntot))

    #add_line(lines, 'OutputDir      @outputs/' + setup_name + '\n')
    add_line(lines, 'OutputDir      {}\n'.format(output.output_dir))

    fname = setup_name + '.par'
    write_file(fname, lines)

    return fname
