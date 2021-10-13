import matplotlib.pyplot as plt
import numpy as np
#from scipy.optimize import fsolve
#from scipy import linalg
#from os import listdir
from scipy.special import roots_legendre
#from scipy.interpolate import BarycentricInterpolator, KroghInterpolator, lagrange

from single_mode import PSI_eigen
from polydust import Polydust
import vistools as vt

def FourierPlot(direcs, Kx, Kz):
    for direc in direcs:
        ret = vt.Fourier(direc, range(0, vt.max_save(direc)), Kx, Kz)

        t = vt.time_stamps(direc, len(ret[:,0]))

        plt.plot(t, np.abs(ret[:,1]))
        plt.plot(t, np.abs(ret[:,2]))
        plt.plot(t, np.abs(ret[:,3]))

    #plt.plot(t, 2.0e-5*np.exp(0.14117601*t))

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
        ret = vt.Fourier(direc, [vt.max_save(direc) - 1], Kx, Kz)

        normfac = 1/ret[-1,1]
        if Kx == 0 and Kz == 0:
            normfac = 1

        pf = vt.PolyFluid(direc)
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
        coord = vt.Coordinates(direc)
        x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
        xmin = x - 0.5*coord.dx
        zmin = z - 0.5*coord.dz

        t = vt.time_stamps(direc, vt.max_save(direc))

        pf = vt.PolyFluid(direc)
        pf.read(0)

        # Initial background
        state0 = vt.Fourier(direc, [0], 0, 0)[0]
        # Initial perturbation
        state1 = vt.Fourier(direc, [0], Kx, Kz)[0]

        ret = np.empty((len(t), 4*(pf.n_dust+1)), dtype=float)

        for n in range(0, vt.max_save(direc)):
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

        plt.plot(t[1:], ret[1:, 1])
        plt.plot(t[1:], ret[1:, 2])
        plt.plot(t[1:], ret[1:, 3])

        print(np.max([ret[:,1], ret[:,2]]))

    plt.yscale('log')
    plt.show()

def contour_dust_stop(direc, number):
    if number < 0:
        number = vt.max_save(direc) - 1

    coord = vt.Coordinates(direc)

    pf = vt.PolyFluid(direc)
    pf.read(number)

    fig, axs = plt.subplots(1,2)

    p = axs[0].contourf(coord.x,coord.z,
                        pf.dust_density()[:,0,:], levels=100)
    plt.colorbar(p, ax=axs[0])

    q = axs[1].contourf(coord.x,coord.z,
                        pf.average_stopping_time()[:,0,:], levels=100)
    plt.colorbar(q, ax=axs[1])

    plt.suptitle(r'$\Omega t =$ {}'.format(vt.time_stamps(direc, number)[-1]))

    plt.show()

def scatter_dust_stop(direc, number):
    coord = vt.Coordinates(direc)

    pf = vt.PolyFluid(direc)
    pf.read(number)

    plt.plot(pf.dust_density()[:,0,:],
             pf.average_stopping_time()[:,0,:], ls='None', marker='o')

    plt.show()

def stopping_time_distribution_max(direc, number, stokes_range,
                                   max_dens=10.0):
    pf = vt.PolyFluid(direc)

    for n in number:
        if n < 0:
            n = vt.max_save(direc) - 1
        pf.read(n)

        dens = []
        sigma = []

        total_dens = pf.dust_density()
        d = np.ravel(total_dens)

        indx = np.asarray(d > max_dens).nonzero()[0]

        # Look at all dust densities
        for i in indx:
            j = np.unravel_index(i, np.shape(total_dens))

            x, y, xi, s = pf.size_distribution(j[0], j[1], j[2])

            tau = stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(xi + 1))
            tau_points = stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(x + 1))

            #s0 = 3/np.sqrt(tau)/2/(np.sqrt(stokes_range[1]) - np.sqrt(stokes_range[0]))

            plt.plot(tau, s)
            plt.plot(tau_points, y, ls='None', marker='o')


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\tau_{\rm s}$')
    plt.ylabel(r'$\sigma$')

    plt.show()

def stopping_time_distribution(direc, number, stokes_range):
    pf = vt.PolyFluid(direc)

    for n in number:
        if n < 0:
            n = vt.max_save(direc) - 1
        pf.read(n)

        dens = []
        sigma = []

        total_dens = pf.dust_density()
        d = np.ravel(total_dens)

        index_array = np.argsort(d)

        # Look at all dust densities
        for i in index_array:
            dens.append(d[i])
            indx = np.unravel_index(i, np.shape(total_dens))

            x, y, xi, s = pf.size_distribution(indx[0], indx[1], indx[2])

            tau = stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(xi + 1))
            #tau_points = stokes_range[0]*np.power(stokes_range[1]/stokes_range[0], 0.5*(x + 1))

            s0 = 3/np.sqrt(tau)/2/(np.sqrt(stokes_range[1]) - np.sqrt(stokes_range[0]))

            sigma.append(s/s0)

            #plt.plot(tau, sigma)
            #plt.plot(tau_points, y, ls='None', marker='o')

        sigma = np.asarray(sigma)

        plt.contourf(tau, dens, sigma,levels=100)

        #for j in range(0, len(sigma[:,0])):
        #    plt.plot(tau, sigma[j,:])
        #print(np.min(sigma), np.max(sigma))

    plt.xlabel(r'$\tau_{\rm s}$')
    plt.ylabel(r'$\rho_{\rm d}$')

    plt.colorbar()
    plt.tight_layout()

    plt.show()


    #plt.plot(tau, 3/np.sqrt(tau)/2/(np.sqrt(stokes_range[1]) - np.sqrt(stokes_range[0])))

    #plt.xlabel(r'$\tau_{\rm s}$')
    #plt.ylabel(r'$\sigma$')

    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

def Reynolds_plot(direcs, cs=20.0):
    for direc in direcs:
        t = vt.time_stamps(direc, vt.max_save(direc))

        coord = vt.Coordinates(direc)
        x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
        vy0 = -1.5*x

        pf = vt.PolyFluid(direc)
        res = []

        for n in range(0, vt.max_save(direc)):
            pf.read(n)

            dens = pf.Fluids[0].dens
            vx = pf.Fluids[0].velx
            vy = pf.Fluids[0].vely # - vy0

            #print(np.min(vy), np.max(vy))

            res.append(np.sum(dens*vx*vy)/np.sum(dens)/cs/cs)

        print(res)
        plt.plot(t, res)

    plt.show()

direcs = [
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K100/N8_ND64/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N32_ND16/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K30/N8_ND32/',
          '/Users/sjp/Codes/psi_fargo/data/mu3_K30_long/N128_ND16/',
          #'/Users/sjp/Codes/psi_fargo/data/mu3_K10/test/',
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


#Reynolds_plot(direcs)
#stopping_time_distribution(direcs[0], [-1], [0.01, 0.1])
#stopping_time_distribution_max(direcs[0], [-1], [0.01, 0.1], max_dens=12.0)

#FourierPlot(direcs, 10, 1)
#EigenVectorPlot(direcs, 100, 200,
#                dust_to_gas_ratio = 3,
#                stokes_range=[1.0e-8,0.1])
#ErrorPlot(direcs, 100, 100,
#          dust_to_gas_ratio = 1.0,
#          stokes_range=[1.0e-4, 0.1])
#exit()


#exit()
contour_dust_stop(direcs[0], -1)
#contour_dust_stop(direcs[1], 900)

#plt.xscale('log')
#plt.show()

#direcs = ['/Users/sjp/Codes/psi_fargo/data/linearA/N8/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N16/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N32/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N64/']

#direc = '/Users/sjp/Codes/psi_fargo/public/outputs/psi_mu2_discrete/'

#coord = Coordinates(direcs[0])
#pf = PolyFluid(direcs[0])
#pf.read(0)

#x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
#vy0 = -1.5*x

#plt.contourf(coord.x, coord.z, pf.Fluids[0].vely[:,0,:] - vy0, levels=100)
#plt.colorbar()

#plt.show()
