from common import *

def add_option(lines, option):
    '''Add option to opt file content'''
    add_line(lines, 'FARGO_OPT +=  -D'+ option + '\n')

def write_opt_file(setup_name, pd):
    n_dust = pd.N
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
    #add_option(lines, 'STANDARD')
    add_option(lines, 'NOSUBSTEP2')

    add_option(lines, 'DRAGFORCE')
    #add_option(lines, 'EXPLICIT_DRAG')
    add_option(lines, 'CONSTANTSTOKESNUMBER')
    #add_option(lines, 'COLLISIONPREDICTOR')

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
