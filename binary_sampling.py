#author: Hsiang-Chih Hwang
#September 2021

import numpy as np
import astropy.units as u
import copy
from numpy import cos, sin
from astropy import constants as const

def random_sample_from_a_power_law(x0, x1, gamma, N_random):
    
    return (np.random.uniform(size=N_random) * (x1**(gamma+1) - x0**(gamma+1)) + x0**(gamma+1))**(1./(gamma+1.))

def random_sample_from_a_straight_line(k, N):
    
    U = np.random.uniform(size=N)
    
    if np.abs(k) < 1e-5:
        return U
    
    elif k < -2: 
                
        return (-np.sqrt(-2.*k) + np.sqrt(-2.*k + 2*k*U)) / k 
    
    elif k > 2.:
                
        return 1. - (-np.sqrt(2.*k) + np.sqrt(2.*k - 2*k*U)) / (-k)
    
    else:
        return (-(1 - 0.5 * k) + np.sqrt((1 - 0.5 * k)**2 + 2 * k * U)) / k
    

def random_true_anomaly(e, N, u_precision = 0.5 * np.pi / 180., equal_spaced_M=False):

    #M: mean anomaly
    if equal_spaced_M:
        M = 2 * np.pi * np.linspace(0., 1., N, endpoint=False)
    else:
        M = 2 * np.pi * np.random.uniform(size=N)
    
    #u: eccentric anomaly, initialized by mean anomaly
    u = copy.copy(M)
    
    #solve for eccentric anomaly iteratively
    #u_precision = 0.5 * np.pi / 180.
    
    for i in range(1, N):
        
        if type(e) is int or type(e) is float or type(e) is np.float64:
            ee = e
        else:
            ee = e[i]
        
        correction = (M[i] - u[i] + ee * sin(u[i])) / (1. - ee * cos(u[i]))
        
        while np.abs(correction) > u_precision:
            u[i] = u[i] + correction
            correction = (M[i] - u[i] + ee * sin(u[i])) / (1. - ee * cos(u[i]))
    
    
    #f: true anomaly
    cosf = (cos(u) - e) / (1 - e * cos(u))
    sinf = (np.sqrt(1-e**2) * sin(u)) / (1 - e * cos(u))
    f = np.mod(np.arctan2(sinf, cosf), 2*np.pi)
    
    return f

def period_cal(a=1*u.au, Mtot=2*u.Msun):
    return (2. * np.pi * np.sqrt(a**3 / const.G / Mtot)).to(u.year)

def period_to_a(P=1.*u.year, Mtot=2*u.Msun):
    return ((P/2./np.pi)**(2./3) * (const.G * Mtot)**(1./3) ).to(u.au)

class two_body_system():
    def compute_components(self):

        self.P = (2. * np.pi * np.sqrt(self.a**3 / (const.G * self.m))).to(u.yr)

        self.p = self.a * (1. - self.e**2)

        self.r0 = self.p / (1. + self.e * cos(self.f))

        self.v0 = (np.sqrt(const.G * self.m / self.p)).to(u.km/u.s)
        self.rx = self.r0 * (cos(self.Omega) * cos(self.omega + self.f) - cos(self.iota) * sin(self.Omega) * sin(self.omega + self.f))
        self.ry = self.r0 * (sin(self.Omega) * cos(self.omega + self.f) + cos(self.iota) * cos(self.Omega) * sin(self.omega + self.f))
        self.rz = self.r0 * sin(self.iota) * sin(self.omega + self.f)

        self.vx = -self.v0 * (cos(self.Omega) * (sin(self.omega+self.f) + self.e * sin(self.omega)) +
               cos(self.iota) * sin(self.Omega) * (cos(self.omega + self.f) + self.e * cos(self.omega))
              )

        self.vy = -self.v0 * (sin(self.Omega) * (sin(self.omega+self.f) + self.e * sin(self.omega)) -
               cos(self.iota) * cos(self.Omega) * (cos(self.omega + self.f) + self.e * cos(self.omega))
              )

        self.vz = self.v0 * sin(self.iota) * (cos(self.omega+self.f) + self.e*cos(self.omega))

        self.cos_phi = (self.rx*self.vx + self.ry*self.vy) / np.sqrt(self.rx**2 + self.ry**2) / np.sqrt(self.vx**2 + self.vy**2)

        self.vr_angle = (np.arccos(self.cos_phi) * 180/np.pi) * u.deg

        self.projected_separation = np.sqrt(self.rx**2 + self.ry**2)
        self.angular_separation = (self.projected_separation / self.distance).to(u.arcsec, equivalencies=u.dimensionless_angles())

        self.projected_velocity = np.sqrt(self.vx**2 + self.vy**2)
        self.pm_diff = (self.projected_velocity / self.distance).to(u.mas/u.yr, equivalencies=u.dimensionless_angles())

 

class binary(two_body_system):
    def __init__(self,
        m1=1.*u.Msun, m2=1.*u.Msun,
        a=1.*u.AU, e=0., distance=1.*u.kpc, Nphase=100,
        u_precision = 0.5 * np.pi / 180.,
        faceon=True
        ):
        self.m1 = m1
        self.m2 = m2
        self.m = self.m1 + self.m2
        self.a = a
        self.e = e
        self.distance = distance
        self.Nphase = Nphase

        self.f = random_true_anomaly(self.e, Nphase, u_precision, equal_spaced_M=True)

        

        if faceon == True:
            self.Omega = np.zeros(Nphase)
            self.omega = np.zeros(Nphase)
            self.iota = np.zeros(Nphase)
        else:
            self.Omega = 2. * np.pi * np.random.uniform(size=Nphase)
            self.omega = 2. * np.pi * np.random.uniform(size=Nphase)
            V = np.random.uniform(size=Nphase)
            self.iota = np.arccos(2. * V - 1)

        self.compute_components()




class binaries(two_body_system):
    def __init__(self, Nbinary, 
        a=None, a_init=None, 
        e=None, e_alpha=None, e_k=None,
        m=None,
        distance=None,
        u_precision = 0.5 * np.pi / 180.,
        faceon=False,
        ):

        self.Nbinary = Nbinary

        if a_init is not None:
            if len(a_init) != 3:
                raise ValueError('len(a_init) != 3')
            a0, a1, gamma = a_init
            self.a = random_sample_from_a_power_law(a0, a1, gamma, Nbinary)
        elif a is None:
            self.a = np.ones(Nbinary) * u.au
        elif isinstance(a, (u.quantity.Quantity)):
            self.a = np.ones(Nbinary) * a
        elif isinstance(a, (list, np.ndarray)):
            if len(a) == Nbinary:
                self.a = a
            else:
                raise ValueError('len(a) not equal to Nbinary')

        if e_alpha is not None:
            if isinstance(e_alpha, (list, np.ndarray)):
                if len(e_alpha) == Nbinary:
                    self.e = random_sample_from_a_power_law(0., 1., e_alpha, Nbinary)
                else:
                    raise ValueError('len(e_alpha) != Nbinary')
            elif isinstance(e_alpha, (int, float)):
                if -1 < e_alpha:
                    self.e = random_sample_from_a_power_law(0., 1., e_alpha, Nbinary)
                else:
                    raise ValueError('e_alpha=%f not >-1' %(e_alpha))

        elif e_k is not None:
            if isinstance(e_k, (list, np.ndarray)):
                if len(e_k) == Nbinary:
                    self.e = random_sample_from_a_straight_line(e_k, Nbinary)
                else:
                    raise ValueError('len(e_k) != Nbinary')
            elif isinstance(e_k, (int, float)):
                self.e = random_sample_from_a_straight_line(e_k, Nbinary)

        elif e is None:
            self.e = np.zeros(Nbinary) 
        else:
            if isinstance(e, (int, float)):
                if 0. <= e <= 1.:
                    self.e = np.ones(Nbinary) * e
                else: 
                    raise ValueError('e=%f not in [0, 1]' %(e))
            elif isinstance(e, (list, np.ndarray)):
                if len(e) == Nbinary:
                    self.e = e
                else:
                    raise ValueError('len(a) not equal to Nbinary')

        if m is None:
            self.m = 2. * u.Msun
        elif isinstance(m, (u.quantity.Quantity)):
            self.m = m
        elif isinstance(m, (list, np.ndarray)):
            if len(m) == Nbinary:
                self.m = m
            else:
                raise ValueError('len(m) not equal to Nbinary')
        


        if distance is None:
            self.distance = np.ones(Nbinary) * u.kpc
        elif isinstance(distance, (u.quantity.Quantity)):
            self.distance = np.ones(Nbinary) * distance
        elif isinstance(distance, (list, np.ndarray)):
            if len(distance) == Nbinary:
                self.distance = distance
            else:
                raise ValueError('len(distance) != Nbinary')
        


        self.f = random_true_anomaly(self.e, Nbinary, u_precision)
        
        if faceon is True:
            self.Omega = np.zeros(Nbinary)
            self.omega = np.zeros(Nbinary)
            self.iota = np.zeros(Nbinary)
        else:
            self.Omega = 2. * np.pi * np.random.uniform(size=Nbinary)
            self.omega = 2. * np.pi * np.random.uniform(size=Nbinary)

            V = np.random.uniform(size=Nbinary)
            self.iota = np.arccos(2. * V - 1)

        self.compute_components()
