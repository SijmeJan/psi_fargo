import matplotlib.pyplot as plt
import numpy as np

class Scalar:
    def __init__(self, filename):
        data = np.loadtxt(filename)
        self.t = data[:,0]
        self.val = data[:,1]

class Coordinates:
    def __init__(self, direc):
        # Coordinates of cell *edges*
        #x = np.loadtxt(direc + 'domain_x.dat')[3:-4]
        y = np.loadtxt(direc + 'domain_y.dat')[3:-4]
        z = np.loadtxt(direc + 'domain_z.dat')[3:-4]

        #dx = x[1] - x[0]
        dy = y[1] - y[0]
        dz = z[1] - z[0]

        # Coordinates of cell *centres*
        #self.x = x + 0.5*dx
        self.y = y + 0.5*dy
        self.z = z + 0.5*dz

class Fluid:
    def __init__(self, direc, number=0):
        coord = Coordinates(direc)

        self.nx = 1
        self.ny = len(coord.y)
        self.nz = len(coord.z)

        self.basename = 'gas'
        if number > 0:
            self.basename = 'dust' + str(number)

        self.basename = direc + self.basename

    def read(self, number):
        fname = self.basename
        self.dens = np.fromfile(fname + 'dens' + str(number) + '.dat')
        self.velx = np.fromfile(fname + 'vx' + str(number) + '.dat')
        self.vely = np.fromfile(fname + 'vy' + str(number) + '.dat')
        self.velz = np.fromfile(fname + 'vz' + str(number) + '.dat')

        self.dens = self.dens.reshape(self.nz, self.ny, self.nx)
        self.velx = self.velx.reshape(self.nz, self.ny, self.nx)
        self.vely = self.vely.reshape(self.nz, self.ny, self.nx)
        self.velz = self.velz.reshape(self.nz, self.ny, self.nx)

direc = '/Users/sjp/Codes/fargo3d/public/outputs/psi_linearA/'

vx = Scalar(direc + 'monitor/gas/momx.dat')
vy = Scalar(direc + 'monitor/gas/momy.dat')
vz = Scalar(direc + 'monitor/gas/momz.dat')

plt.plot(vz.t, vz.val)
#plt.plot(vy.t, vy.val)

plt.plot(vz.t, 1.0e-11*np.exp(0.84*vz.t))
plt.yscale('log')
plt.show()

exit()

coord = Coordinates(direc)

gas = Fluid(direc)
gas.read(0)
dust = Fluid(direc, number=1)
dust.read(1)

#nx = 1
#ny = len(coord.y)
#nz = len(coord.z)

#rho = np.fromfile(direc + 'gasvy0.dat').reshape(nz,ny,nx)

plt.contourf(coord.y, coord.z, dust.dens[:,:,0], levels=100)
plt.colorbar()

plt.show()
