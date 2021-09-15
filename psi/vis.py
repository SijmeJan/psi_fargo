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
        dx = x[1] - x[0]
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

    dx = coord.x[1] - coord.x[0]
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

        for indx, fluid in enumerate(pf.Fluids):
            ret[t,4*indx+0] = np.mean(fluid.dens[:,0,:]*expmed)
            ret[t,4*indx+1] = np.mean(fluid.velx[:,0,:]*expminx)
            ret[t,4*indx+2] = np.mean(fluid.vely[:,0,:]*expmed)
            ret[t,4*indx+3] = np.mean(fluid.velz[:,0,:]*expminz)

    return ret

direcs = ['/Users/sjp/Codes/psi_fargo/public/outputs/psi_mu3/']

pf = PolyFluid(direcs[0])

#direcs = ['/Users/sjp/Codes/psi_fargo/data/psi_mu10/N64/']

#direcs = ['/Users/sjp/Codes/psi_fargo/data/linearA/N8/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N16/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N32/',
#          '/Users/sjp/Codes/psi_fargo/data/linearA/N64/']

for direc in direcs:
    ret = Fourier(direc, range(0, max_save(direc)), 0, 0)

    data = np.genfromtxt(direc + '/variables.par',dtype='str')
    dt = 0.0
    n = 1
    for d in data:
        if d[0] == 'DT':
            dt = float(d[1])
        if d[0] == 'NINTERM':
            n = int(d[1])

    print(dt, n)
    t = dt*n*np.arange(len(ret[:,4]))


    #for i in range(0, len(ret[0,:])):
    #    plt.plot(t, np.abs(ret[:,i]))


    #for i in range(0, np.int(len(ret[0,:])/4)):
    #    plt.plot(t, np.abs(ret[:,4*i+1]))
    plt.plot(pf.stopping_times, np.abs(ret[0,5::4])) # vx as a function of size
    plt.plot(pf.stopping_times, np.abs(ret[-1,5::4]))
#plt.plot(t, 1.0e-4*np.exp(0.3027*t))

plt.xscale('log')
plt.show()

exit()


#direc = '/Users/sjp/Codes/psi_fargo/data/linearA/N32/'
#vz = Scalar(direc + 'monitor/gas/momz.dat')
#plt.plot(vz.t, vz.val)

#plt.plot(vz.t, 1.0e-11*np.exp(0.84*vz.t))

#direc = '/Users/sjp/Codes/psi_fargo/data/linearA/N64/'
#vz = Scalar(direc + 'monitor/gas/momz.dat')
#plt.plot(vz.t, vz.val)

#direc = '/Users/sjp/Codes/psi_fargo/data/linearA/N128/'
#vz = Scalar(direc + 'monitor/gas/momz.dat')
#plt.plot(vz.t, vz.val)

#direc = '/Users/sjp/Codes/psi_fargo/data/linearA/N256/'
#vz = Scalar(direc + 'monitor/gas/momz.dat')
#plt.plot(vz.t, vz.val)

direc = '/Users/sjp/Codes/psi_fargo/public/outputs/psi_linearA/'
vz = Scalar(direc + 'monitor/gas/momz.dat')
plt.plot(vz.t, vz.val)
plt.plot(vz.t, 1.0e-11*np.exp(2*0.0154862262*vz.t))


#vz = Scalar(direc + 'monitor/gas/momz.dat')
#plt.plot(vz.t, vz.val)

plt.yscale('log')

AN = 3*0.1/(1 + 0.1*0.1)
BN = 1 + 3/(1 + 0.1*0.1)

plt.show()

#exit()

#direc = '/Users/sjp/Codes/psi_fargo/data/linearA/N64/'

coord = Coordinates(direc)

gas = Fluid(direc)
gas.read(0)
dust = Fluid(direc, number=1)
dust.read(0)
#print('Initial dvx:', np.mean(gas.vely-dust.vely))
#gas.read(200)
#dust.read(200)
#print('Final dvx:', np.mean(gas.vely-dust.vely))



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

def veq_ana(eps, T):
    A = eps*T/(1+T*T)
    B = 1 + eps/(1+T*T)

    vgx = 2*A/(A*A + B*B)
    vgy = -B/(A*A + B*B)

    vdx = (vgx + 2*T*vgy)/(1+T*T)
    vdy = (vgy - 0.5*T*vgx)/(1 + T*T)

    return np.array([vgx, vgy, vdx, vdy])

def eq(p):
    eps, T = p

    v_num = veq_num(eps, T)
    v_ana = veq_ana(eps, T)

    return (v_num[0] - v_ana[0], v_num[1] - v_ana[1])


#veq = veq_num(3, 0.1)
#print(veq[0] - veq[2])

#import numpy as np
#rng = np.random.default_rng(12345)

#min_y = 0

#for times in range(0, 100000):
#    eps = 10*rng.random()
#    T = rng.random()

#    ret, info, err, msg =  fsolve(eq, (eps, T), full_output=1)

#    x = ret[0]
#    y = ret[1]

#    if x > 0.1 and x < 100 and y > min_y and y < 1.0 and err==1:
#        print(x, y, veq_num(x,y), veq_ana(x, y))

#        min_y = y


#nx = 1
#ny = len(coord.y)
#nz = len(coord.z)

#rho = np.fromfile(direc + 'gasvy0.dat').reshape(nz,ny,nx)

#plt.contourf(coord.y, coord.z, dust.dens[:,:,0], levels=100)
plt.plot(coord.y, dust.dens[:,0,0])
plt.plot(coord.y + 2.0*np.pi/6.0, dust.dens[:,0,0])
#plt.colorbar()

plt.show()
