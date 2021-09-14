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
    def __init__(self, n_dust, stokes_range, size_distribution=None):
        self.N = n_dust
        self.stokes_range = np.asarray(stokes_range)

        self.sigma = None
        # Continuous size distribution
        if isinstance(size_distribution, SizeDistribution):
            self.sigma = size_distribution.sigma
        else:
            if size_distribution is None:
                # MRN, normalized to unity
                size_distribution = SizeDistribution(stokes_range)
            else:
                # Discrete sizes passed, size should be the same as stokes_range
                self.sigma = np.asarray(size_distribution)

    def dust_nodes(self):
        if callable(self.sigma):
            xi, weights = roots_legendre(self.N)

            xi = np.asarray(xi)
            weights = np.asarray(weights)

            q = self.stokes_range[1]/self.stokes_range[0]

            # Stopping time nodes
            tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))
        else:
            # Discrete dust sizes
            tau = self.stokes_range

        return tau

    def dust_densities(self, dust_total_density):
        if callable(self.sigma):
            xi, weights = roots_legendre(self.N)

            xi = np.asarray(xi)
            weights = np.asarray(weights)

            q = self.stokes_range[1]/self.stokes_range[0]

            # Stopping time nodes
            tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))

            # 'density' at nodes
            dens = 0.5*np.log(q)*weights*tau*dust_total_density*self.sigma(tau)
        else:
            # Discrete dust sizes
            tau = self.stokes_range
            dens = self.sigma

        return tau, dens

    def equilibrium_velocities(self, tau, dust_to_gas_ratio):
        if callable(self.sigma):
            f = lambda x: self.sigma(x)/(1.0 + x*x)
            J0 = quad(f, self.stokes_range[0], self.stokes_range[1])[0]
            f = lambda x: self.sigma(x)*x/(1.0 + x*x)
            J1 = quad(f, self.stokes_range[0], self.stokes_range[1])[0]

            J0 = J0*dust_to_gas_ratio
            J0 = J0*dust_to_gas_ratio

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

            AN = np.sum(self.sigma*tau/(1 + tau*tau))
            BN = 1.0 + np.sum(self.sigma/(1 + tau*tau))

            v = []
            v.append(2*AN/(AN*AN + BN*BN))   # Gas vx
            v.append(-BN/(AN*AN + BN*BN))    # Gas vy
            v.append(0.0)                    # Gas vz

            for t in tau:
                v.append((v[0] + 2*t*v[1])/(1 + t*t))
                v.append((v[1] - 0.5*t*v[0])/(1 + t*t))
                v.append(0.0)

        return np.asarray(v)

    def initial_conditions(self, dust_total_density, gas_density):
        tau, dens = self.dust_densities(dust_total_density)

        mu = dust_total_density/gas_density
        v = self.equilibrium_velocities(tau, mu)

        return tau, dens, v
