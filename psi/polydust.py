import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad

class SizeDistribution():
    def __init__(self, stokes_range):
        self.taumin = stokes_range[0]
        self.taumax = stokes_range[1]

    def sigma(self, x):
        '''MRN size distribution, normalized to unity'''
        return 0.5/(np.sqrt(x)*(np.sqrt(self.taumax) - np.sqrt(self.taumin)))

class Polydust():
    def __init__(self, n_dust, stokes_range,
                 dust_density, gas_density,
                 size_distribution=None,
                 gauss_legendre=True):
        self.N = n_dust
        self.stokes_range = np.asarray(stokes_range)
        self.dust_density = dust_density
        self.gas_density = gas_density
        self.gauss_legendre = gauss_legendre

        self.sigma = None
        # Continuous size distribution
        if isinstance(size_distribution, SizeDistribution):
            print('Continuous size distribution with {} nodes'.format(self.N))
            if self.gauss_legendre is True:
                print('Using Gauss-Legendre quadrature')
            else:
                print('Using equidistant nodes')
            self.sigma = size_distribution.sigma
        else:
            if size_distribution is None:
                # MRN, normalized to unity
                print('Continuous MRN size distribution with {} nodes'.format(self.N))
                if self.gauss_legendre is True:
                    print('Using Gauss-Legendre quadrature')
                else:
                    print('Using equidistant nodes')
                size_distribution = SizeDistribution(stokes_range)
            else:
                print('Discrete dust sizes with {} dust species'.format(self.N))
                self.N = len(stokes_range)
                # Discrete sizes passed, size should be the same as stokes_range
                self.sigma = np.asarray(size_distribution)
                # Normalize to total dust density unity
                self.sigma = self.sigma/np.sum(self.sigma)

    def dust_nodes(self):
        if callable(self.sigma):
            if self.gauss_legendre is True:
                xi, weights = roots_legendre(self.N)

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                q = self.stokes_range[1]/self.stokes_range[0]

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))
            else:
                # Constant spacing in log tau space
                #xi1 = np.log10(self.stokes_range[0])
                #xi2 = np.log10(self.stokes_range[1])

                #tau = np.logspace(xi1, xi2, self.N)

                xi1 = np.log(self.stokes_range[0])
                xi2 = np.log(self.stokes_range[1])

                tau = np.exp(np.linspace(xi1, xi2, self.N))

                xi_edge = np.linspace(xi1, xi2, self.N + 1)
                tau = xi_edge + 0.5*(xi_edge[1] - xi_edge[0])
                tau = np.exp(tau)[0:-1]

        else:
            # Discrete dust sizes
            tau = self.stokes_range

        return tau

    def dust_densities(self):
        if callable(self.sigma):
            if self.gauss_legendre is True:
                xi, weights = roots_legendre(self.N)

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                q = self.stokes_range[1]/self.stokes_range[0]

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))

                # 'density' at nodes
                dens=0.5*np.log(q)*weights*tau*self.dust_density*self.sigma(tau)
            else:
                # Constant spacing in log tau space
                xi1 = np.log(self.stokes_range[0])
                xi2 = np.log(self.stokes_range[1])

                xi_edge = np.linspace(xi1, xi2, self.N + 1)
                dxi = xi_edge[1] - xi_edge[0]
                tau = np.exp(xi_edge + 0.5*dxi)[0:-1]

                dens = self.dust_density*self.sigma(tau)*dxi*tau
       else:
            # Discrete dust sizes
            tau = self.stokes_range
            # Normalize so total dust_density is correct
            dens = self.sigma*self.dust_density

        print('Total dust density: ', np.sum(dens))

        return tau, dens

    def equilibrium_velocities(self, tau):
        if callable(self.sigma):
            f = lambda x: self.sigma(x)/(1.0 + x*x)
            J0 = quad(f, self.stokes_range[0], self.stokes_range[1])[0]
            f = lambda x: self.sigma(x)*x/(1.0 + x*x)
            J1 = quad(f, self.stokes_range[0], self.stokes_range[1])[0]

            J0 = J0*self.dust_density/self.gas_density
            J1 = J1*self.dust_density/self.gas_density

            denom = (1 + J0)*(1 + J0) + J1*J1
            v = []
            v.append(2*J1/denom)      # Gas vx
            v.append(-(1 + J0)/denom) # Gas vy
            v.append(0.0)             # Gas vz

            for t in tau:
                v.append(2*(J1 - t*(1 + J0))/((1 + t*t)*denom))
                v.append(-(1 + J0 + t*J1)/((1 + t*t)*denom))
                v.append(0.0)
        else:
            tau = np.asarray(self.stokes_range)

            mu = self.dust_density/self.gas_density
            AN = mu*np.sum(self.sigma*tau/(1 + tau*tau))
            BN = 1.0 + mu*np.sum(self.sigma/(1 + tau*tau))

            v = []
            v.append(2*AN/(AN*AN + BN*BN))   # Gas vx
            v.append(-BN/(AN*AN + BN*BN))    # Gas vy
            v.append(0.0)                    # Gas vz

            for t in tau:
                v.append((v[0] + 2*t*v[1])/(1 + t*t))
                v.append((v[1] - 0.5*t*v[0])/(1 + t*t))
                v.append(0.0)

        return np.asarray(v)

    def initial_conditions(self):
        tau, dens = self.dust_densities()

        v = self.equilibrium_velocities(tau)

        return tau, dens, v
