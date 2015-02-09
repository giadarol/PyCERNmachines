from __future__ import division
import numpy as np
from scipy.constants import c, e, m_p

from PyHEADTAIL.trackers.transverse_tracking import TransverseMap
from PyHEADTAIL.trackers.detuners import Chromaticity, AmplitudeDetuning
from PyHEADTAIL.trackers.simple_long_tracking import LinearMap, RFSystems
from PyHEADTAIL.particles.generators import CutRFBucket6D


class synchrotron(object):
	
	def __init__(self, *args, **kwargs):
		
			
		for attr in kwargs.keys():
			if kwargs[attr] is not None:
				print 'Synchrotron init. From kwargs: %s = %s'%(attr, repr(kwargs[attr]))
				setattr(self, attr, kwargs[attr])
			
		self.create_transverse_map()
		self.create_longitudinal_map()

		self.one_turn_map = []
		for m in self.transverse_map:
			self.one_turn_map.append(m)
		self.one_turn_map.append(self.longitudinal_map)

		
	
	def install_after_each_transverse_segment(self, element_to_add):
		one_turn_map_new = []
		for element in self.one_turn_map:
			one_turn_map_new.append(element)
			if element in self.transverse_map:
				one_turn_map_new.append(element_to_add)
		self.one_turn_map = one_turn_map_new
		
	@property
	def gamma(self):
		return self._gamma
	@gamma.setter
	def gamma(self, value):
		self._gamma = value
		self._beta = np.sqrt(1 - self.gamma**-2)
		self._betagamma = np.sqrt(self.gamma**2 - 1)
		self._p0 = self.betagamma * m_p * c

	@property
	def beta(self):
		return self._beta
	@beta.setter
	def beta(self, value):
		self.gamma = 1. / np.sqrt(1 - value ** 2)

	@property
	def betagamma(self):
		return self._betagamma
	@betagamma.setter
	def betagamma(self, value):
		self.gamma = np.sqrt(value**2 + 1)

	@property
	def p0(self):
		return self._p0
	@p0.setter
	def p0(self, value):
		self.gamma = value / (m_p * self.beta * c)
		
	@property
	def eta(self):
		return self.alpha - self.gamma**-2
	
	@property
	def Q_s(self):
		if self.V2!=0. or self.p_increment!=0 or self.dphi1!=0:
			raise ValueError('Formula not valid in this case!!!!') 
		return np.sqrt( e*np.abs(self.eta)*(self.h1*self.V1 + self.h2*self.V2)
						/(2*np.pi*self.p0*self.beta*c) )

	def track(self, bunch):
		for m in self.one_turn_map:
			m.track(bunch)
			
	def create_transverse_map(self):
		self.transverse_map = TransverseMap(
			self.circumference, self.s, self.alpha_x, self.beta_x, self.D_x, self.alpha_y, self.beta_y, 
			self.D_y, self.Q_x, self.Q_y,
			Chromaticity(self.Qp_x, self.Qp_y),
			AmplitudeDetuning(self.app_x, self.app_y, self.app_xy))

		

	def create_longitudinal_map(self):

		if self.longitudinal_focusing == 'linear':
			raise ValueError('''longitudinal_focusing='linear' not tested!!!''')
			#~ eta = self.alpha - self.gamma**-2
			#~ p0 = np.sqrt(self.gamma**2 - 1) * m_p * c
			#~ Q_s = np.sqrt( e*np.abs(eta)*(self.h1*self.V1 + self.h2*self.V2)
						 #~ /(2*np.pi*p0*beta*c) )
			#~ self.longitudinal_map = LinearMap([self.alpha], self.circumference, self.Q_s())
		elif self.longitudinal_focusing == 'non-linear':
			self.longitudinal_map = RFSystems(self.circumference, [self.h1, self.h2], [self.V1, self.V2], [self.dphi1, self.dphi2],
										[self.alpha], self.gamma, self.p_increment)
		else:
			raise ValueError('ERROR: unknown focusing', self.longitudinal_focusing)

	def generate_6D_Gaussian_bunch(self, n_macroparticles, intensity, epsn_x, epsn_y, sigma_z):
		if self.longitudinal_focusing == 'linear':
			raise ValueError('This mode is not compatible with the current implementations')
		beta_z    = self.eta*self.circumference/2./np.pi/self.Q_s
		sigma_dp  = sigma_z/beta_z
		return CutRFBucket6D(macroparticlenumber=n_macroparticles, intensity=intensity, charge=self.charge, mass=self.mass,
                 circumference = self.circumference, gamma_reference=self.gamma,
                 transverse_map=self.transverse_map, epsn_x=epsn_x, epsn_y=epsn_y,
                 sigma_z=sigma_z, sigma_dp=sigma_dp,
                 is_accepted=self.longitudinal_map.rfbucket.make_is_accepted(margin = 0.05)).generate()
      
