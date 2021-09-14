import numpy as np

from common import *

def veq_num(eps, T):
    # Numerical equilibrium (vy's without background shear)
    #2*vgy + 2*chi*R - eps*(vgx - vdx)/T = 0
    #2*vgx + eps*(vgy - vdy)/T = 0
    #2*vdy - (vdx - vgx)/T = 0
    #2*vdx + (vdy - vgy)/T = 0

    A = np.array([[-eps/T, 2, eps/T, 0], [2, eps/T, 0, -eps/T],
                  [1/T, 0, -1/T, 2], [0, -1/T, 2, 1/T]])
    b = np.array([-2, 0, 0, 0])

    return np.linalg.solve(A, b)

def write_condinit_file(pd, dust_total_density, perturbation=None):
    stokes, sigma = pd.dust_densities(dust_total_density)

    # Numerical equilibrium
    #veq = veq_num(sigma[0], stokes[0])
    #vgx = veq[0]
    #vgy = veq[1]
    #vgz = 0.0
    #vdx = np.array([veq[2]])
    #vdy = np.array([veq[3]])
    #vdz = 0.0*stokes

    AN = np.sum(sigma*stokes/(1 + stokes*stokes))
    BN = 1.0 + np.sum(sigma/(1 + stokes*stokes))

    vgx = 2*AN/(AN*AN + BN*BN)
    vgy = -BN/(AN*AN + BN*BN)     # Excluding background shear
    vgz = 0.0

    vdx = (vgx + 2*stokes*vgy)/(1 + stokes*stokes)
    vdy = (vgy - 0.5*stokes*vgx)/(1 + stokes*stokes)
    vdz = 0.0*stokes

    # Convert to FARGO standard where y=x and x=y....
    vgx, vgy = vgy, vgx
    vdx, vdy = vdy, vdx

    pert = perturbation
    if perturbation is None:
        pert=[]
        for i in range(0, 4*n_dust+4):
            pert.append('')

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
#       '        rho[l] = 1.0 + '+ str_perturbation(amp, Kx, Kz, drhog) + ';\n',
#       '        vx[l]  = {} + '.format(vgx) + str_perturbation(amp, Kx, Kz, dvgx) + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
#       '        vy[l]  = {} + '.format(vgy) + str_perturbation(amp, Kx, Kz, dvgy, stagger='y') + ';\n',
#       '        vz[l]  = {} + '.format(vgz) + str_perturbation(amp, Kx, Kz, dvgz, stagger='z') + ';\n',
       '        rho[l] = 1.0 '+ pert[0] + ';\n',
       '        vx[l]  = {} '.format(vgx) + pert[1] + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
       '        vy[l]  = {} '.format(vgy) + pert[2] + ';\n',
       '        vz[l]  = {} '.format(vgz) + pert[3] + ';\n',
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
#           '        rho[l] = {} + '.format(sigma[n]) + str_perturbation(amp, Kx, Kz, drhod[n]) + ';\n',
#           '        vx[l]  = {} + '.format(vdx[n]) + str_perturbation(amp, Kx, Kz, dvdx[n]) + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
#           '        vy[l]  = {} + '.format(vdy[n]) + str_perturbation(amp, Kx, Kz, dvdy[n], stagger='y') + ';\n',
#           '        vz[l]  = {} + '.format(vdz[n]) + str_perturbation(amp, Kx, Kz, dvdz[n], stagger='z') + ';\n',
           '        rho[l] = {} '.format(sigma[n]) + pert[4*n+4] + ';\n',
           '        vx[l]  = {} '.format(vdx[n]) + pert[4*n+5] + ' - SHEARPARAM*OMEGAFRAME*OMEGAFRAME*Ymed(j);\n',
           '        vy[l]  = {} '.format(vdy[n]) + pert[4*n+6] + ';\n',
           '        vz[l]  = {} '.format(vdz[n]) + pert[4*n+7] + ';\n',
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
