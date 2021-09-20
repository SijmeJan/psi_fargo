import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from os import listdir

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

        #dx = x[1] - x[0]
        dx = 0.0
        if len(x) > 1:
            dx = x[1] - x[0]
        dz = 0.0
        if len(z) > 1:
            dz = z[1] - z[0]

        # Coordinates of cell *centres*
        #self.x = x + 0.5*dx
        self.x = x + 0.5*dx
        self.z = z + 0.5*dz

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

def max_save(direc):
    n = 0
    for file in listdir(direc):
        if 'gasdens' in file:
            n += 1

    return n

def Fourier(direc, time_stamps, Kx, Kz):
    coord = Coordinates(direc)

    dx = 0.0
    if len(coord.x) > 1:
        dx = coord.x[1] - coord.x[0]
    dz = 0.0
    if len(coord.z) > 1:
        dz = coord.z[1] - coord.z[0]

    x, z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')

    xmin = x - 0.5*dx
    zmin = z - 0.5*dz

    expmed = np.exp(-1j*(Kx*x + Kz*z))
    expminx = np.exp(-1j*(Kx*xmin + Kz*z))
    expminz = np.exp(-1j*(Kx*x + Kz*zmin))

    pf = PolyFluid(direc)

    ret = np.empty((len(time_stamps), 4*(pf.n_dust+1)), dtype=np.complex128)

    for t, n in enumerate(time_stamps):
        pf.read(n)

        #print(t, pf.Fluids[0].velx)

        for indx, fluid in enumerate(pf.Fluids):
            ret[t,4*indx+0] = np.mean(fluid.dens[:,0,:]*expmed)
            ret[t,4*indx+1] = np.mean(fluid.velx[:,0,:]*expminx)
            ret[t,4*indx+2] = np.mean(fluid.vely[:,0,:]*expmed)
            ret[t,4*indx+3] = np.mean(fluid.velz[:,0,:]*expminz)

    return ret

def FourierPlot(direcs, Kx, Kz):
    for direc in direcs:
        ret = Fourier(direc, range(0, max_save(direc)), Kx, Kz)

        data = np.genfromtxt(direc + '/variables.par',dtype='str')
        dt = 0.0
        n = 1
        for d in data:
            if d[0] == 'DT':
                dt = float(d[1])
            if d[0] == 'NINTERM':
                n = int(d[1])

        t = dt*n*np.arange(len(ret[:,4]))

        #for i in range(1, int(len(ret[0,:])/4)):
        #    plt.plot(t, np.abs(ret[:,4*i]))
        plt.plot(t, np.abs(ret[:,1]))
        #plt.plot(t, np.real(ret[:,2] - ret[0,2]))

        #plt.plot(pf.stopping_times, np.abs(ret[-1,5::4]))

        #plt.plot(pf.stopping_times, np.real(ret[0,5::4])) # vx as a function of size
        #plt.plot(pf.stopping_times, np.real(ret[-1,5::4]))
        #tau = pf.stopping_times
        #plt.plot(pf.stopping_times, -2*tau/(1 + tau*tau))

    #plt.plot(t, 1.0e-5*np.exp(0.4190091323*t))
    plt.plot(t, 1.0e-5*np.exp(0.3027262829*t))
    #plt.plot(t, 1.0e-5*np.exp(0.0980250*t))
    #plt.plot(t, 1.0e-5*np.exp(0.42185357935887124*t))

    plt.yscale('log')
    plt.show()

direcs = [
          '/Users/sjp/Codes/psi_fargo/public/outputs/psi_linearA/'
          #'/Users/sjp/Codes/psi_fargo/data/psi_linearA/'
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

FourierPlot(direcs, 60, 60)

#for direc in direcs:
#    pf = PolyFluid(direc)

#    ret = Fourier(direc, range(0, max_save(direc)), 60, 60)

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
