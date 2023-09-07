import numpy as np
from scipy import linalg

try:
    from psitools.psi_mode import PSIMode
except ImportError as error:
    # as of python 3.6 now throws ModuleNotFoundError
    print('Will run despite error importing psitools:', error)


def PSI_eigen(dust_to_gas_ratio, stokes_range, wave_number_x, wave_number_z,
              viscous_alpha=0.0):
    '''Calculate PSI eigenfunctions

    Args:
        dust_to_gas_ratio: Background dust to gas ratio
        stokes_range: minimum and maximum Stokes number
        wave_number_x: Kx
        wave_number_z: Kz
        viscous_alpha (optional): viscosity parameter, defaults to zero.
    '''
    # Use PSIMode to calculate eigenvalue
    np.random.seed(0)
    pm = PSIMode(dust_to_gas_ratio=dust_to_gas_ratio,
                 stokes_range=stokes_range,
                 real_range=[-2, 2],
                 imag_range=[1.0e-8, 1])

    roots = pm.calculate(wave_number_x=wave_number_x,
                         wave_number_z=wave_number_z,
                         viscous_alpha=viscous_alpha,
                         constant_schmidt=True)

    # Get corresponding eigenfunctions
    rhog, vg, sigma, u = \
      pm.eigenvector(roots[0], wave_number_x=wave_number_x,
                     wave_number_z=wave_number_z,
                     viscous_alpha=viscous_alpha,
                     constant_schmidt=True)

    # Set local so we can access eigenvalue if needed
    PSI_eigen.eigenvalue = roots[0]
    print('Calculated PSI eigenvalue: ', roots[0])

    return rhog, vg, sigma, u

class Perturbation():
    '''Base class for PSI perturbations

    Args:
        n_dust: number of dust collocation points
        amplitude: amplitude of perturbation
        Kx: wave number x
        Kz: wave number z
    '''
    def __init__(self, n_dust, amplitude, Kx, Kz):
        self.Kx = Kx
        self.Kz = Kz

        self.n_dust = n_dust
        self.amp = amplitude

        self.eig = np.zeros((4*self.n_dust+4), dtype=np.complex128)
        self.set_eigenvector()

    def set_eigenvector(self):
        '''Set eigenvector dummy, should be overridden'''
        self.eig = 0*self.eig

    def string_single(self, eig, stagger=None):
        '''Return string amp*(r*cos(Kx*x+Kz*z) - i*sin(Kx*x+Kz*z))'''
        k = '{}*Ymed(j) + {}*Zmed(k)'.format(self.Kx, self.Kz)
        # Stagger in y/z if required
        if stagger == 'y':
            k = '{}*Ymin(j) + {}*Zmed(k)'.format(self.Kx, self.Kz)
        if stagger == 'z':
            k = '{}*Ymed(j) + {}*Zmin(k)'.format(self.Kx, self.Kz)

        ret = '+ {}*({}*cos('.format(self.amp, np.real(eig)) + k
        ret += ') - ({})*sin('.format(np.imag(eig)) + k + '))'
        return ret

    def to_string(self):
        '''Return list of strings to be used in condinit.c file'''
        ret = []

        for i in range(0, self.n_dust+1):
            ret.append(self.string_single(self.eig[4*i]))
            ret.append(self.string_single(self.eig[4*i+1]))
            ret.append(self.string_single(self.eig[4*i+2],stagger='y'))
            ret.append(self.string_single(self.eig[4*i+3],stagger='z'))

        return ret

class GasEpicycle(Perturbation):
    '''Gas-only epicycle.

    Args:
        Kx: wave number x
        Kz: wave number z
        sound_speed: sound speed
        amplitude: epicycle amplitude
    '''

    def __init__(self, Kx, Kz, sound_speed, amplitude):
        self.c = sound_speed
        Perturbation.__init__(self, 0, amplitude, Kx, Kz)

    def set_eigenvector(self):
        A = [[0, self.Kx, 0, self.Kz],
             [self.c*self.c*self.Kx, 0, 2j, 0],
             [0, -0.5j, 0, 0],
             [self.Kz*self.c*self.c, 0, 0, 0]]

        # Select eigenvalue with real part closest to unity
        ev = linalg.eigvals(np.asarray(A))
        n = np.abs(ev.real - 1.0).argmin()

        val, vec = linalg.eig(np.asarray(A))
        print('Eigenvector for eigenvalue {}: {}'.format(val[n], vec[:,n]))

        v = vec[:, n]

        # Convert to FARGO standard where y=x and x=y....
        v[2], v[1] = v[1], v[2]

        self.eig = v

class LinearA(Perturbation):
    '''Class for LinearA SI test problem.

    Args:
        amplitude: perturbation amplitude
    '''
    def __init__(self, amplitude):
        # LinearA: 1 dust fluid, wave number = 30
        Perturbation.__init__(self, 1, amplitude, 30, 30)

    def set_eigenvector(self):
        drhog = 7.46086390e-06+7.06921219e-06j
        dvgx = -5.63886397e-02+1.20336579e-02j
        dvgy =  4.45522906e-02+1.97146371e-02j
        dvgz =  5.63883479e-02-1.20337366e-02j
        drhod = 1.0
        dvdx =  -4.66288875e-02+1.24126320e-02j
        dvdy =  4.35180720e-02+2.13451298e-02j
        dvdz =  5.46594210e-02-7.75824909e-03j

        # Convert to FARGO standard where y=x and x=y....
        dvgx, dvgy = dvgy, dvgx
        dvdx, dvdy = dvdy, dvdx

        self.eig = [drhog, dvgx, dvgy, dvgz, drhod, dvdx, dvdy, dvdz]

class LinearB(Perturbation):
    '''Class for LinearB SI test problem.

    Args:
        amplitude: perturbation amplitude
    '''
    def __init__(self, amplitude):
        # LinearB: 1 dust fluid, wave number = 6
        Perturbation.__init__(self, 1, amplitude, 6, 6)

    def set_eigenvector(self):
        # Eigenvector LinB
        drhog = -0.0000337227 - 0.0003456248j
        dvgx = -0.0870451125 - 1.3851731095j
        dvgy = +1.3839936168 - 0.0937424679j
        dvgz = +0.0870497444 + 1.3852113520j
        drhod = 1.0
        dvdx =  +0.2314730923 - 1.3715260043j
        dvdy =  +1.3696536978 + 0.0196879160j
        dvdz =  +0.0416164539 + 1.3844311928j

        # Convert to FARGO standard where y=x and x=y....
        dvgx, dvgy = dvgy, dvgx
        dvdx, dvdy = dvdy, dvdx

        self.eig = [drhog, dvgx, dvgy, dvgz, drhod, dvdx, dvdy, dvdz]

class Linear3(Perturbation):
    '''Class for Linear3 MSI test problem.

    Args:
        amplitude: perturbation amplitude
    '''
    def __init__(self, amplitude):
        # LinearB: 2 dust fluids, wave number = 50
        Perturbation.__init__(self, 2, amplitude, 50, 50)

    def set_eigenvector(self):
        self.eig = [+0.0000061052 + 0.0000080743j,
                    +0.1327989476 + 0.0674232641j,
                    -0.1587288108 + 0.0213251096j,
                    +0.1587286212 - 0.0213252588j,
                    1.0,
                    +0.1325843682 + 0.0691301709j,
                    -0.1461274403 + 0.0234873672j,
                    +0.1571142133 - 0.0174328415j,
                    +0.1522281314 + 0.1836379253j,
                    +0.1092222067 + 0.0952973332j,
                    -0.1335593453 + 0.0025396632j,
                    +0.1485545469 + 0.0200753935j]

class RandomFixedK(Perturbation):
    '''Class for random eigenvector at specific wave number.

    Args:
        n_dust: number of dust collocation points
        amplitude: amplitude of perturbation
        Kx: wave number x
        Kz: wave number z
    '''
    def set_eigenvector(self):
        rng = np.random.default_rng(12345)

        for i in range(0, len(self.eig)):
            real_part = rng.random()
            imag_part = rng.random()

            self.eig[i] = real_part + 1j*imag_part

        # Make sure gas density perturbation is tiny (almost incompressible)
        self.eig[0] = 1.0e-4*self.eig[0]

        # Normalized to first dust density perturbation
        self.eig[4] = 1.0

class PSI_pert(Perturbation):
    '''Class for exact PSI perturbation

    Args:
        poly_dust: PolyDust object
        amplitude: perturbation amplitude
        Kx: wave number x
        Kz: wave number z
    '''
    def __init__(self, poly_dust, amplitude, Kx, Kz, viscous_alpha=0.0):
        self.pd = poly_dust
        self.viscous_alpha = viscous_alpha

        Perturbation.__init__(self, self.pd.N, amplitude, Kx, Kz)

    def set_eigenvector(self):
        # Dust to gas ratio
        mu = self.pd.dust_density/self.pd.gas_density

        # PSI eigenvector
        rhog, vg, sigma, u = \
          PSI_eigen(mu, self.pd.stokes_range, self.Kx, self.Kz,
                    viscous_alpha=self.viscous_alpha)

        tau = self.pd.dust_nodes()

        weights = np.empty(len(tau))
        weights.fill(1)

        sigma_norm = lambda x: sigma(x)
        if self.pd.gauss_legendre is True:
            fac = 0.5*np.log(self.pd.stokes_range[1]/self.pd.stokes_range[0])
            sigma_norm = lambda x: fac*x*sigma(x)
            xi, weights = self.pd.nodes_and_weights()
        else:
            xi = np.log(tau)
            dxi = xi[1] - xi[0]
            sigma_norm = lambda x: x*sigma(x)*dxi

        # NOTE: convert to FARGO standard where x=y...
        # NOTE 2: WHY THE FACTOR PI????
        self.eig = [rhog, vg[1], vg[0], vg[2]]

        for x,w in zip(tau,weights):
            self.eig.extend([np.pi*w*sigma_norm(x), u[1](x), u[0](x), u[2](x)])
