import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy import linalg
from os import listdir
from scipy.special import roots_legendre
from scipy.interpolate import BarycentricInterpolator, KroghInterpolator, lagrange

from single_mode import PSI_eigen
from polydust import Polydust

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

    def average_stopping_time(self):
        num = self.Fluids[1].dens*self.stopping_times[0]
        denom = self.Fluids[1].dens

        for n in range(2, self.n_dust + 1):
            num = num + self.Fluids[n].dens*self.stopping_times[n-1]
            denom = denom + self.Fluids[n].dens

        return num/denom

    def dust_density(self):
        num = self.Fluids[1].dens

        for n in range(2, self.n_dust + 1):
            num = num + self.Fluids[n].dens

        return num

    def size_distribution(self, i, j, k):
        # How to reconstruct the size distribution from nodal values?
        xi, weights = roots_legendre(self.n_dust)
        xi = np.asarray(xi)
        weights = np.asarray(weights)

        # Roundabout way to calc 0.5*log(taumax/taumin)
        logfac = \
          np.log(self.stopping_times[1]/self.stopping_times[0])/(xi[1]-xi[0])

        x = xi #self.stopping_times
        y = [self.Fluids[1].dens[i,j,k]]

        for n in range(2, self.n_dust +1):
            y.append(self.Fluids[n].dens[i,j,k])

        # Convert to sigma
        y = np.asarray(y)/weights/self.stopping_times/logfac

        xi = np.linspace(-1, 1, 100)
        res = BarycentricInterpolator(x, y)(xi)

        return x, y, xi, res


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
        plt.plot(t, np.abs(ret[:,2]))
        plt.plot(t, np.abs(ret[:,3]))


    #plt.plot(t, 1.0e-7*np.exp(0.4190091323*t))
    #plt.plot(t, 1.0e-5*np.exp(0.3027262829*t))
    #plt.plot(t, 1.0e-5*np.exp(0.0980250*t))
    #plt.plot(t, 1.0e-5*np.exp(0.17*t))
    #plt.plot(t, 1.0e-5*np.exp(0.42185357935887124*t))
    plt.plot(t, 2.0e-5*np.exp(0.14117601*t))

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
    plt.xscale('log')

    plt.subplot(2,2,2)
    plt.plot(tau, np.real(u[0](tau)))
    plt.plot(tau, np.imag(u[0](tau)))
    plt.xscale('log')

    plt.subplot(2,2,3)
    plt.plot(tau, np.real(u[1](tau)))
    plt.plot(tau, np.imag(u[1](tau)))
    plt.xscale('log')

    plt.subplot(2,2,4)
    plt.plot(tau, np.real(u[2](tau)))
    plt.plot(tau, np.imag(u[2](tau)))
    plt.xscale('log')

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
    rhog, vg, sigma, u = \
          PSI_eigen(dust_to_gas_ratio=dust_to_gas_ratio,
                    stokes_range=stokes_range,
                    wave_number_x=Kx,
                    wave_number_z=Kz)
    w = PSI_eigen.eigenvalue

    #w = -1j*(-0.3027262829 + 0.3242790653j)
    #w = 0.3733047861642763+0.0021788291419847284j

    for direc in direcs:
        coord = Coordinates(direc)
        x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
        xmin = x - 0.5*coord.dx
        zmin = z - 0.5*coord.dz

        t = time_stamps(direc, max_save(direc))

        pf = PolyFluid(direc)
        pf.read(0)

        #pd = Polydust(pf.n_dust, stokes_range,
        #              dust_to_gas_ratio, 1.0)

        #v_eq = pd.equilibrium_velocities(pf.stopping_times)

        # Initial background
        state0 = Fourier(direc, [0], 0, 0)[0]
        # Initial perturbation
        state1 = Fourier(direc, [0], Kx, Kz)[0]

        #state0[1] = v_eq[0]
        #state0[2] = v_eq[1]

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

        plt.plot(t[1:], ret[1:, 1])
        plt.plot(t[1:], ret[1:, 2])
        plt.plot(t[1:], ret[1:, 3])


        print(np.max([ret[:,1], ret[:,2]]))

    plt.yscale('log')
    plt.show()

def contour_dust_stop(direc, number):
    coord = Coordinates(direc)

    pf = PolyFluid(direc)
    pf.read(number)

    fig, axs = plt.subplots(1,2)

    p = axs[0].contourf(coord.x,coord.z,
                        pf.dust_density()[:,0,:], levels=100)
    plt.colorbar(p, ax=axs[0])

    q = axs[1].contourf(coord.x,coord.z,
                        pf.average_stopping_time()[:,0,:], levels=100)
    plt.colorbar(q, ax=axs[1])

    plt.show()

def scatter_dust_stop(direc, number):
    coord = Coordinates(direc)

    pf = PolyFluid(direc)
    pf.read(number)

    plt.plot(pf.dust_density()[:,0,:],
             pf.average_stopping_time()[:,0,:], ls='None', marker='o')

    plt.show()

def stopping_time_distribution(direc, number, stokes_range):
    pf = PolyFluid(direc)

    for n in number:
        pf.read(n)

        # Look at maximum dust density
        d = pf.dust_density()[:,0,:]
        indx = np.argmax(d)
        indx = np.unravel_index(indx, np.shape(d))
        print(d[indx])

        x, y, xi, sigma = pf.size_distribution(indx[0],0,indx[1])

        tau = \
          stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(xi + 1))
        tau_points = \
          stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(x + 1))

        plt.plot(tau, sigma)
        plt.plot(tau_points, y, ls='None', marker='o')

        # Look at minimum dust density
        d = pf.dust_density()[:,0,:]
        indx = np.argmin(d)
        indx = np.unravel_index(indx, np.shape(d))


        x, y, xi, sigma = pf.size_distribution(indx[0],0,indx[1])

        tau = \
          stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(xi + 1))
        tau_points = \
          stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(x + 1))

    #plt.plot(tau, sigma)
    #plt.plot(tau_points, y, ls='None', marker='o')

    plt.plot(tau, 3/np.sqrt(tau)/2/(np.sqrt(stokes_range[1]) - np.sqrt(stokes_range[0])))

    plt.xlabel(r'$\tau_{\rm s}$')
    plt.ylabel(r'$\sigma$')

    plt.xscale('log')
    plt.yscale('log')
    plt.show()

direcs = [
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K100/N8_ND64/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N32_ND16/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N8_ND32/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30_long/N32_ND8/',
          '/Users/sjp/Codes/psi_fargo/data/mu1_K100/test/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30_long/N32_ND4/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30_long/N64_ND4/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N64_ND4_gauss/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N128_ND8_gauss/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N8/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N16/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N32/',
          #'/Users/sjp/Codes/psi_fargo/data/lin3/N64/',
          #'/Users/sjp/Codes/psi_fargo/public/outputs/psi_linearA/'
          ]

#stopping_time_distribution(direcs[0], [700,750,800,850,900,950,999], [0.01, 0.1])

FourierPlot(direcs, 100, 100)
#EigenVectorPlot(direcs, 100, 200,
#                dust_to_gas_ratio = 3,
#                stokes_range=[1.0e-8,0.1])
#ErrorPlot(direcs, 100, 100,
#          dust_to_gas_ratio = 1.0,
#          stokes_range=[1.0e-4, 0.1])
#exit()


#exit()
#contour_dust_stop(direcs[0], 500)
#contour_dust_stop(direcs[1], 900)

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
