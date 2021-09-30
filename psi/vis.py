import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy import linalg
from os import listdir
from scipy.special import roots_legendre

from single_mode import PSI_eigen

class Scalar:
    def __init__(self, filename):
        data = np.loadtxt(filename)
        self.t = data[:,0]
        self.val = data[:,1]

class Coordinates:
    def __init__(self, direc):
        # Coordinates of cell *edges*
        #x = np.loadtxt(direc + 'domain_x.dat')[3:-4]
        x = np.loadtxt(direc + 'domain_y.dat')[3:-4]
        z = np.loadtxt(direc + 'domain_z.dat')[3:-4]

        dx = 0.0
        if len(x) > 1:
            dx = x[1] - x[0]
        dz = 0.0
        if len(z) > 1:
            dz = z[1] - z[0]

        # Coordinates of cell *centres*
        self.x = x + 0.5*dx
        self.z = z + 0.5*dz

        self.dx = dx
        self.dz = dz

class Fluid:
    def __init__(self, direc, number=0):
        coord = Coordinates(direc)

        self.ny = 1
        self.nx = len(coord.x)
        self.nz = len(coord.z)

        self.basename = 'gas'
        if number > 0:
            self.basename = 'dust' + str(number)

        self.basename = direc + self.basename

    def read(self, number):
        fname = self.basename
        self.dens = np.fromfile(fname + 'dens' + str(number) + '.dat')
        self.vely = np.fromfile(fname + 'vx' + str(number) + '.dat')
        self.velx = np.fromfile(fname + 'vy' + str(number) + '.dat')
        self.velz = np.fromfile(fname + 'vz' + str(number) + '.dat')

        self.dens = self.dens.reshape(self.nz, self.ny, self.nx)
        self.velx = self.velx.reshape(self.nz, self.ny, self.nx)
        self.vely = self.vely.reshape(self.nz, self.ny, self.nx)
        self.velz = self.velz.reshape(self.nz, self.ny, self.nx)

class PolyFluid:
    def __init__(self, direc):

        data = np.genfromtxt(direc + '/variables.par',dtype='str')

        tau = []
        for d in data:
            if d[0].find('INVSTOKES') != -1:
                tau.append(float(d[1]))
        tau = np.sort(1.0/np.asarray(tau))

        self.stopping_times = tau
        self.n_dust = len(tau)

        # Gas fluid
        self.Fluids = [Fluid(direc, number=0)]

        # Dust fluids
        for n in range(1, self.n_dust + 1):
            self.Fluids.append(Fluid(direc, number=n))

    def read(self, number):
        for fluid in self.Fluids:
            fluid.read(number)

def time_stamps(direc, length):
    data = np.genfromtxt(direc + '/variables.par',dtype='str')
    dt = 0.0
    n = 1
    for d in data:
        if d[0] == 'DT':
            dt = float(d[1])
        if d[0] == 'NINTERM':
            n = int(d[1])

    return dt*n*np.arange(length)

def max_save(direc):
    n = 0
    for file in listdir(direc):
        if 'gasdens' in file:
            n += 1

    return n

def Fourier(direc, time_stamps, Kx, Kz):
    coord = Coordinates(direc)

    x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')

    xmin = x - 0.5*coord.dx
    zmin = z - 0.5*coord.dz

    expmed = np.exp(-1j*(Kx*x + Kz*z))
    expminx = np.exp(-1j*(Kx*xmin + Kz*z))
    expminz = np.exp(-1j*(Kx*x + Kz*zmin))

    pf = PolyFluid(direc)

    ret = np.empty((len(time_stamps), 4*(pf.n_dust+1)), dtype=np.complex128)

    for t, n in enumerate(time_stamps):
        pf.read(n)

        for indx, fluid in enumerate(pf.Fluids):
            ret[t,4*indx+0] = np.mean(fluid.dens[:,0,:]*expmed)
            ret[t,4*indx+1] = np.mean(fluid.velx[:,0,:]*expminx)
            ret[t,4*indx+2] = np.mean(fluid.vely[:,0,:]*expmed)
            ret[t,4*indx+3] = np.mean(fluid.velz[:,0,:]*expminz)

    if (Kx == 0 and Kz == 0):
        return ret
    return 2*ret

def FourierPlot(direcs, Kx, Kz):
    for direc in direcs:
        ret = Fourier(direc, range(0, max_save(direc)), Kx, Kz)

        t = time_stamps(direc, len(ret[:,0]))

        #for i in range(1, int(len(ret[0,:])/4)):
        #    plt.plot(t, np.abs(ret[:,4*i]))
        plt.plot(t, np.abs(ret[:,1]))
        #plt.plot(t, np.real(ret[:,2]))


    #plt.plot(t, 1.0e-7*np.exp(0.4190091323*t))
    #plt.plot(t, 1.0e-5*np.exp(0.3027262829*t))
    #plt.plot(t, 1.0e-5*np.exp(0.0980250*t))
    plt.plot(t, 1.0e-5*np.exp(0.17*t))
    #plt.plot(t, 1.0e-5*np.exp(0.42185357935887124*t))
    #plt.plot(t, 1.0e-6*np.exp(0.0154839*t))

    plt.yscale('log')
    plt.show()

def EigenVectorPlot(direcs, Kx, Kz,
                    dust_to_gas_ratio = 2.0,
                    stokes_range=[0.001,0.1]):
    rhog, vg, sigma, u = PSI_eigen(dust_to_gas_ratio=dust_to_gas_ratio,
                                   stokes_range=stokes_range,
                                   wave_number_x=Kx,
                                   wave_number_z=Kz)

    tau = np.logspace(np.log10(stokes_range[0]),
                      np.log10(stokes_range[1]), 1000)

    plt.subplot(2,2,1)
    plt.plot(tau, np.real(np.pi*sigma(tau)))
    plt.plot(tau, np.imag(np.pi*sigma(tau)))

    plt.subplot(2,2,2)
    plt.plot(tau, np.real(u[0](tau)))
    plt.plot(tau, np.imag(u[0](tau)))

    plt.subplot(2,2,3)
    plt.plot(tau, np.real(u[1](tau)))
    plt.plot(tau, np.imag(u[1](tau)))

    plt.subplot(2,2,4)
    plt.plot(tau, np.real(u[2](tau)))
    plt.plot(tau, np.imag(u[2](tau)))

    for direc in direcs:
        ret = Fourier(direc, [max_save(direc) - 1], Kx, Kz)

        normfac = 1/ret[-1,1]
        if Kx == 0 and Kz == 0:
            normfac = 1

        pf = PolyFluid(direc)
        tau = pf.stopping_times

        if (np.abs(np.log(tau[1]*tau[1]/tau[0]/tau[2])) < 1.0e-10):
            print(direc, ': equidistant nodes')

            # Equidistant
            dlog = np.log(tau[1]/tau[0])
            sigma = ret[-1,4::4]*normfac/tau/dlog
        else:
            print(direc, ': Gauss-Legendre nodes')
            # Gauss-Legendre
            xi, w = roots_legendre(len(tau))
            sigma = 2*ret[-1,4::4]*normfac/tau/w/np.log(stokes_range[1]/stokes_range[0])


        ux = ret[-1,5::4]*normfac
        uy = ret[-1,6::4]*normfac
        uz = ret[-1,7::4]*normfac

        rhog = ret[-1,0]*normfac
        vgx = ret[-1,1]*normfac
        vgy = ret[-1,2]*normfac
        vgz = ret[-1,3]*normfac

        plt.subplot(2,2,1)
        plt.plot(tau, np.real(sigma))
        plt.plot(tau, np.imag(sigma))
        plt.plot([np.min(tau)], [np.real(rhog)], marker='o')
        plt.plot([np.min(tau)], [np.imag(rhog)], marker='o')
        plt.xscale('log')
        plt.ylabel(r'$\hat\varsigma$')

        plt.subplot(2,2,2)
        plt.plot(tau, np.real(ux))
        plt.plot(tau, np.imag(ux))
        plt.plot([np.min(tau)], [np.real(vgx)], marker='o')
        plt.plot([np.min(tau)], [np.imag(vgx)], marker='o')
        plt.xscale('log')
        plt.ylabel(r'$\hat u_x$')
        print('Maximum norm ux error:', np.max(np.abs(ux - u[0](tau))))

        plt.subplot(2,2,3)
        plt.plot(tau, np.real(uy))
        plt.plot(tau, np.imag(uy))
        plt.plot([np.min(tau)], [np.real(vgy)], marker='o')
        plt.plot([np.min(tau)], [np.imag(vgy)], marker='o')
        plt.xscale('log')
        plt.ylabel(r'$\hat u_y$')

        plt.subplot(2,2,4)
        plt.plot(tau, np.real(uz))
        plt.plot(tau, np.imag(uz))
        plt.plot([np.min(tau)], [np.real(vgz)], marker='o')
        plt.plot([np.min(tau)], [np.imag(vgz)], marker='o')
        plt.xscale('log')
        plt.ylabel(r'$\hat u_z$')

    plt.tight_layout()

    plt.xscale('log')
    plt.show()

def ErrorPlot(direcs, Kx, Kz,
              dust_to_gas_ratio = 2.0,
              stokes_range=[0.001,0.1]):
    #rhog, vg, sigma, u = \
    #      PSI_eigen(dust_to_gas_ratio=dust_to_gas_ratio,
    #                stokes_range=stokes_range,
    #                wave_number_x=Kx,
    #                wave_number_z=Kz)
    #w = PSI_eigen.eigenvalue

    w = -1j*(-0.3027262829 + 0.3242790653j)

    for direc in direcs:
        coord = Coordinates(direc)
        x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
        xmin = x - 0.5*coord.dx
        zmin = z - 0.5*coord.dz

        t = time_stamps(direc, max_save(direc))

        pf = PolyFluid(direc)
        pf.read(0)

        # Initial background
        state0 = Fourier(direc, [0], 0, 0)[0]
        # Initial perturbation
        state1 = Fourier(direc, [0], Kx, Kz)[0]

        ret = np.empty((len(t), 4*(pf.n_dust+1)), dtype=float)

        for n in range(0, max_save(direc)):
            pf.read(n)

            expmed = np.exp(1j*(Kx*x + Kz*z - w*t[n]))
            expminx = np.exp(1j*(Kx*xmin + Kz*z - w*t[n]))
            expminz = np.exp(1j*(Kx*x + Kz*zmin - w*t[n]))

            for i, fluid in enumerate(pf.Fluids):
                dens_ana = np.real(state0[4*i] + state1[4*i]*expmed)
                velx_ana = np.real(state0[4*i+1] + state1[4*i+1]*expminx)
                vely_ana = np.real(state0[4*i+2] + state1[4*i+2]*expmed) - 1.5*x
                velz_ana = np.real(state0[4*i+3] + state1[4*i+3]*expminz)

                ret[n, 4*i] = np.max(np.abs(fluid.dens[:,0,:] - dens_ana))
                ret[n, 4*i+1] = np.max(np.abs(fluid.velx[:,0,:] - velx_ana))
                ret[n, 4*i+2] = np.max(np.abs(fluid.vely[:,0,:] - vely_ana))
                ret[n, 4*i+3] = np.max(np.abs(fluid.velz[:,0,:] - velz_ana))

                #if i == 0 and n==0:
                #    plt.subplot(121)
                #    plt.contourf(fluid.velz[:,0,:])
                #    plt.colorbar()
                #    plt.subplot(122)
                #    plt.contourf(velz_ana)
                #    plt.colorbar()

                #    plt.show()

        #plt.plot(t[1:], ret[1:, 1])
        plt.plot(t[1:], ret[1:, 2])
        #plt.plot(t[1:], ret[1:, 3])


        print(np.max([ret[:,1], ret[:,2]]))

    plt.yscale('log')
    plt.show()

direcs = [
          #'/Users/sjp/Codes/psi_fargo/data/psi_mu2/N32_ND32_gauss/',
          #'/Users/sjp/Codes/psi_fargo/data/psi_mu2/N32_ND16_gauss/',
          '/Users/sjp/Codes/psi_fargo/data/test/N64_ND8_gauss/',
          #'/Users/sjp/Codes/psi_fargo/data/test/N8_ND4_gauss_cfl044/',
          #'/Users/sjp/Codes/psi_fargo/data/test/N8_ND4_gauss_cfl0044/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N8/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N16/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N32/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N64/',
          #'/Users/sjp/Codes/psi_fargo/public/outputs/psi_linearA/'
          ]

#direcs = [
          #'/Users/sjp/Codes/psi_fargo/data/psi_mu2/N32/',
#          '/Users/sjp/Codes/psi_fargo/data/psi_mu2_discrete/N32_ND8/',
#          '/Users/sjp/Codes/psi_fargo/data/psi_mu2_discrete/N32_ND16/',
#          '/Users/sjp/Codes/psi_fargo/data/psi_mu2_discrete/N32_ND32/',
#          '/Users/sjp/Codes/psi_fargo/data/psi_mu2_discrete/N32_ND64/'
          #'/Users/sjp/Codes/psi_fargo/data/psi_mu2/N32_ND16/',
#          '/Users/sjp/Codes/psi_fargo/data/psi_mu2/N32_ND32/'
#         ]

#FourierPlot(direcs, 30, 30)
EigenVectorPlot(direcs, 30, 30,
                dust_to_gas_ratio = 3,
                stokes_range=[0.01,0.1])
#ErrorPlot(direcs, 50, 50,
#          dust_to_gas_ratio = 1.5,
#          stokes_range=[0.01,0.1])
exit()

#exit()

pf = PolyFluid(direcs[0])
pf.read(0)

for direc in direcs:
    pf = PolyFluid(direc)

    data = np.genfromtxt(direc + '/variables.par',dtype='str')
    dt = 0.0
    n = 1
    for d in data:
        if d[0] == 'DT':
            dt = float(d[1])
        if d[0] == 'NINTERM':
            n = int(d[1])

    t = dt*n*np.arange(max_save(direc))

    gasdens = []
    dustdens = []
    for n in range(0, max_save(direc)):
        pf.read(n)

        gasdens.append(pf.Fluids[0].dens[0, 0, 0] - 1.0)
        dustdens.append(pf.Fluids[1].dens[0, 0, 0] - 3.0)

    gasdens = np.asarray(gasdens/np.max(gasdens))
    dustdens = np.asarray(dustdens/np.max(dustdens))

    plt.plot(t, gasdens)
    plt.plot(t, dustdens)

plt.show()
    #ret = Fourier(direc, range(0, max_save(direc)), 60, 60)

#    plt.plot(pf.stopping_times, np.real(ret[-1,4::4]))

exit()



#plt.xscale('log')
#plt.show()

#direcs = ['/Users/sjp/Codes/psi_fargo/data/linearA/N8/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N16/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N32/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N64/']

#direc = '/Users/sjp/Codes/psi_fargo/public/outputs/psi_mu2_discrete/'

#coord = Coordinates(direc)

#gas = Fluid(direc)
#gas.read(max_save(direc)-1)
#dust = Fluid(direc, number=1)
#dust.read(max_save(direc)-1)

#plt.contourf(coord.x,coord.z, dust.velx[:,0,:], levels=100)
#plt.colorbar()

#plt.show()
