import numpy as np

class Perturbation():
    def __init__(self, n_dust, amplitude, Kx, Kz):
        self.Kx = Kx
        self.Kz = Kz

        self.n_dust = n_dust
        self.amp = amplitude

        self.eig = np.zeros((4*self.n_dust+4), dtype=np.complex128)
        self.set_eigenvector()

    def set_eigenvector(self):
        self.eig = 0*self.eig

    def string_single(self, eig, stagger=None):
        '''Return string amp*(r*cos(Kx*x+Kz*z) - i*sin(Kx*x+Kz*z))'''
        k = '{}*Ymed(j) + {}*Zmed(k)'.format(self.Kx, self.Kz)
        # Stagger in y
        if stagger == 'y':
            k = '{}*Ymin(j) + {}*Zmed(k)'.format(self.Kx, self.Kz)
        if stagger == 'z':
            k = '{}*Ymed(j) + {}*Zmin(k)'.format(self.Kx, self.Kz)

        ret = '+ {}*({}*cos('.format(self.amp, np.real(eig)) + k
        ret += ') - ({})*sin('.format(np.imag(eig)) + k + '))'
        return ret

    def to_string(self):
        ret = []

        for i in range(0, self.n_dust+1):
            ret.append(self.string_single(self.eig[4*i]))
            ret.append(self.string_single(self.eig[4*i+1]))
            ret.append(self.string_single(self.eig[4*i+2],stagger='y'))
            ret.append(self.string_single(self.eig[4*i+3],stagger='z'))

        return ret

class LinearA(Perturbation):
    def __init__(self, amplitude):
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
    def __init__(self, amplitude):
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
    def __init__(self, amplitude):
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
    def set_eigenvector(self):
        rng = np.random.default_rng(12345)

        for i in range(0, len(self.eig)):
            real_part = rng.random()
            imag_part = rng.random()

            self.eig[i] = real_part + 1j*imag_part

        # Make sure gas density perturbation is tiny (almost incompressible)
        self.eig[0] = 1.0e-6*self.eig[0]

        # Normalized to first dust density perturbation
        self.eig[4] = 1.0
