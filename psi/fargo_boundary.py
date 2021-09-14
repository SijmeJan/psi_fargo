from common import *

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
