import numpy as np
import os
from scipy.special import roots_legendre
from scipy.integrate import quad

#from multistream import MSI

class SizeDistribution():
    '''Class holding a size (or stopping time) distribution

    Args:
        stokes_range: minimum and maximum Stokes number
    '''
    def __init__(self, stokes_range):
        self.taumin = stokes_range[0]
        if len(stokes_range) > 1:
            self.taumax = stokes_range[1]
        else:
            self.taumax = self.taumin

    def sigma(self, x):
        '''MRN size distribution, normalized to unity'''
        return 0.5/(np.sqrt(x)*(np.sqrt(self.taumax) - np.sqrt(self.taumin)))

class Polydust():
    def __init__(self, n_dust, stokes_range,
                 dust_density, gas_density,
                 size_distribution=None,
                 gauss_legendre=True,
                 discrete_equilibrium=False):
        '''Class for polydisperse dust component in FARGO3D simulations.

        Args:
            n_dust: Number of collocation points ('dust fluids')
            stokes_range: Range of Stokes numbers to consider. In case of discrete sizes, it can be a list of Stokes numbers.
            dust_density: background dust density
            gas_density: background gas density
            size_distribution (optional): Either a SizeDistribution object, or a list of dust densities, in which case the numer of elements should be equal to the number of elements of stokes_range. Defaults to None, which is a continuous MRN size distribution.
            gauss_legendre (optional): Flag whether to use Gauss-Legendre collocation points. If False, use equidistant nodes in log(tau) space. Defaults to True.
            discrete_equilibrium (optional): In case of equidistant nodes, force a discrete equilibrium a la Krapp et al. and Chao-Chin & Zhu. Defaults to False.
        '''

        self.N = n_dust
        self.stokes_range = np.asarray(stokes_range)
        self.dust_density = dust_density
        self.gas_density = gas_density
        self.gauss_legendre = gauss_legendre
        self.discrete_equilibrium = discrete_equilibrium

        # Set the size density function
        self.sigma = None
        if isinstance(size_distribution, SizeDistribution):
            # Continuous size distribution
            print('Continuous size distribution with {} nodes'.format(self.N))
            if self.gauss_legendre is True:
                print('Using Gauss-Legendre quadrature')
            else:
                print('Using equidistant nodes')
                if self.discrete_equilibrium is True:
                    print('Using discrete equilibrium')
            self.sigma = size_distribution.sigma
        else:
            if size_distribution is None:
                # MRN, normalized to unity
                print('Continuous MRN size distribution with {} nodes'.format(self.N))
                if self.gauss_legendre is True:
                    print('Using Gauss-Legendre quadrature')
                else:
                    print('Using equidistant nodes')
                    if self.discrete_equilibrium is True:
                        print('Using discrete equilibrium')
                size_distribution = SizeDistribution(stokes_range)
                self.sigma = size_distribution.sigma
            else:
                print('Discrete dust sizes with {} dust species'.format(self.N))
                self.N = len(stokes_range)
                # Discrete sizes passed, size should be the same as stokes_range
                self.sigma = np.asarray(size_distribution)
                # Normalize to total dust density unity
                self.sigma = self.sigma/np.sum(self.sigma)

    @classmethod
    def from_output_direc(cls, direc):
        n_dust = len([name for name in os.listdir(direc) if name.endswith('dens0_0.dat')]) - 1

        print('Number of dust species:', n_dust)

        # NOTE: Stokes range and gas and dust densities not read!
        return cls(n_dust, [0.01, 0.1], 3.0, 1.0)

    def nodes_and_weights(self):
        '''Return Gauss-Legendre nodes and weights'''
        return roots_legendre(self.N)

    def dust_nodes(self):
        '''Return collocation points'''
        if callable(self.sigma):
            # Continuous size distribution
            if self.gauss_legendre is True:
                # Gauss-Legendre collocation points
                xi, weights = self.nodes_and_weights()

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                if len(self.stokes_range) > 1:
                    q = self.stokes_range[1]/self.stokes_range[0]
                else:
                    q = 1.0

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))
            else:
                # Constant spacing in log tau space
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
        '''Return densities at collocation points'''
        if callable(self.sigma):
            # Continuous size distribution
            if self.gauss_legendre is True:
                # Gauss-Legende collocation points
                xi, weights = roots_legendre(self.N)

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                if len(self.stokes_range) > 1:
                    q = self.stokes_range[1]/self.stokes_range[0]
                else:
                    q = 1.0

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

                # Set dust density so mass in size bin is exact MRN
                if (self.discrete_equilibrium is True):
                    t_l = np.exp(0.5*xi_edge)[0:-1]
                    t_r = np.exp(0.5*np.roll(xi_edge, -1))[0:-1]
                    dens = self.dust_density*(t_r - t_l)/(t_r[-1] - t_l[0])

                    # Make sure self.sigma is no longer callable: discrete!
                    self.sigma = dens/self.dust_density
                    self.stokes_range = tau

                    #print(self.sigma, tau)

        else:
            # Discrete dust sizes
            tau = self.stokes_range
            # Normalize so total dust_density is correct
            dens = self.sigma*self.dust_density

        print('Total dust density: ', np.sum(dens))

        return tau, dens

    def equilibrium_velocities(self, tau):
        '''Return equilibrium dust and gas velocities

        Args:
           tau: list of stopping times
        '''
        if callable(self.sigma):
            # Continuous size distribution
            f = lambda x: np.exp(x)*self.sigma(np.exp(x))/(1.0 + np.exp(2*x))
            J0 = quad(f, np.log(self.stokes_range[0]),
                         np.log(self.stokes_range[1]))[0]
            f = lambda x: np.exp(2*x)*self.sigma(np.exp(x))/(1.0 + np.exp(2*x))
            J1 = quad(f, np.log(self.stokes_range[0]),
                         np.log(self.stokes_range[1]))[0]

            J0 = J0*self.dust_density/self.gas_density
            J1 = J1*self.dust_density/self.gas_density

            denom = (1 + J0)*(1 + J0) + J1*J1
            v = []
            v.append(2*J1/denom)      # Gas vx
            v.append(-(1 + J0)/denom) # Gas vy
            v.append(0.0)             # Gas vz

            # Dust velocities
            for t in tau:
                v.append(2*(J1 - t*(1 + J0))/((1 + t*t)*denom))
                v.append(-(1 + J0 + t*J1)/((1 + t*t)*denom))
                v.append(0.0)
        else:
            # Discrete dust sizes, override argument
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
        '''Return FARGO initial conditions (equilibrium)'''

        # Collocation points and densities
        tau, dens = self.dust_densities()

        # Equilibrium velocities
        v = self.equilibrium_velocities(tau)

        return tau, dens, v

    #def growth_rate(self, Kx, Kz):
    #    tau, dens = self.dust_densities()
    #
    #    return MSI(dens, tau).max_growth(Kx, Kz)
