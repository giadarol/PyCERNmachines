from __future__ import division

import numpy as np
from scipy.constants import c, e, m_p

from machines import synchrotron


class PSB(synchrotron):

    def __init__(self, *args, **kwargs):

        if 'n_segments' not in kwargs.keys():
            raise ValueError('Number of segments must be specified')

        if 'machine_configuration' not in kwargs.keys():
            raise ValueError('machine_configuration must be specified')

        self.n_segments = kwargs['n_segments']
        self.machine_configuration = kwargs['machine_configuration']

        self.circumference  = 50*np.pi
        self.s = (np.arange(0, self.n_segments + 1)
                  * self.circumference / self.n_segments)

        if self.machine_configuration == '160MeV':
            self.charge = e
            self.mass = m_p

            self.gamma = 160e6*e/(self.mass*c**2) + 1

            self.Q_x     = 4.23
            self.Q_y     = 4.37

            self.Qp_x    = -1*self.Q_x
            self.Qp_y    = -2*self.Q_y

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha_x = 0 * np.ones(self.n_segments)
            self.beta_x  = self.circumference/(2*np.pi*self.Q_x) * np.ones(self.n_segments)
            self.D_x     = 0 * np.ones(self.n_segments)
            self.alpha_y = 0 * np.ones(self.n_segments)
            self.beta_y  = self.circumference/(2*np.pi*self.Q_y) * np.ones(self.n_segments)
            self.D_y     = 0 * np.ones(self.n_segments)

            self.alpha       = 0.06
            self.h1          = 1
            self.h2          = 2
            self.V1          = 8e3
            self.V2          = 0
            self.dphi1       = 0
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'

        elif self.machine_configuration == '1GeV':
            self.charge = e
            self.mass = m_p

            self.gamma = 1e9*e/(self.mass*c**2) + 1

            self.Q_x     = 4.23
            self.Q_y     = 4.37

            self.Qp_x    = -1*self.Q_x
            self.Qp_y    = -2*self.Q_y

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha_x = 0 * np.ones(self.n_segments)
            self.beta_x  = self.circumference/(2*np.pi*self.Q_x) * np.ones(self.n_segments)
            self.D_x     = 0 * np.ones(self.n_segments)
            self.alpha_y = 0 * np.ones(self.n_segments)
            self.beta_y  = self.circumference/(2*np.pi*self.Q_y) * np.ones(self.n_segments)
            self.D_y     = 0 * np.ones(self.n_segments)

            self.alpha       = 0.06
            self.h1          = 1
            self.h2          = 2
            self.V1          = 8e3
            self.V2          = 0
            self.dphi1       = 0
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'

        elif self.machine_configuration == '1.4GeV':
            self.charge = e
            self.mass = m_p

            self.gamma = 1.4e9*e/(self.mass*c**2) + 1

            self.Q_x     = 4.23
            self.Q_y     = 4.37

            self.Qp_x    = -1*self.Q_x
            self.Qp_y    = -2*self.Q_y

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha_x = 0 * np.ones(self.n_segments)
            self.beta_x  = self.circumference/(2*np.pi*self.Q_x) * np.ones(self.n_segments)
            self.D_x     = 0 * np.ones(self.n_segments)
            self.alpha_y = 0 * np.ones(self.n_segments)
            self.beta_y  = self.circumference/(2*np.pi*self.Q_y) * np.ones(self.n_segments)
            self.D_y     = 0 * np.ones(self.n_segments)

            self.alpha       = 0.06
            self.h1          = 1
            self.h2          = 2
            self.V1          = 8e3
            self.V2          = 0
            self.dphi1       = 0
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'
        else:
            raise ValueError('Unknown configuration '+self.machine_configuration)

        super(PSB, self).__init__(*args, **kwargs)


class SPS(synchrotron):

    def __init__(self, *args, **kwargs):

        if 'n_segments' not in kwargs.keys():
            raise ValueError('Number of segments must be specified')

        if 'machine_configuration' not in kwargs.keys():
            raise ValueError('machine_configuration must be specified')

        self.n_segments = kwargs['n_segments']
        self.machine_configuration = kwargs['machine_configuration']

        self.circumference = 1100*2*np.pi

        self.s = np.arange(0, self.n_segments + 1) * self.circumference / self.n_segments

        if self.machine_configuration == 'Q20-injection':
            self.charge = e
            self.mass = m_p

            self.alpha_x        = 0 * np.ones(self.n_segments)
            self.beta_x         = 54.6 * np.ones(self.n_segments)
            self.D_x            = 0 * np.ones(self.n_segments)
            self.alpha_y        = 0 * np.ones(self.n_segments)
            self.beta_y         = 54.6 * np.ones(self.n_segments)
            self.D_y            = 0 * np.ones(self.n_segments)

            self.Q_x            = 20.13
            self.Q_y            = 20.18

            self.Qp_x           = 0
            self.Qp_y           = 0

            self.app_x          = 0.0000e-9
            self.app_y          = 0.0000e-9
            self.app_xy         = 0

            self.alpha     = 0.00308

            self.h1, self.h2       = 4620, 4620*4
            self.V1, self.V2       = 5.75e6, 0
            self.dphi1, self.dphi2 = 0, np.pi
            self.p_increment       = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'
            self.gamma = 27.7

            self.p_increment       = 0 * e/c * self.circumference/(self.beta*c)

        elif self.machine_configuration =='Q26-injection':
            self.charge = e
            self.mass = m_p

            self.alpha_x        = 0 * np.ones(self.n_segments)
            self.beta_x         = 42. * np.ones(self.n_segments)
            self.D_x            = 0 * np.ones(self.n_segments)
            self.alpha_y        = 0 * np.ones(self.n_segments)
            self.beta_y         = 42. * np.ones(self.n_segments)
            self.D_y            = 0 * np.ones(self.n_segments)

            self.Q_x            = 26.13
            self.Q_y            = 26.18

            self.Qp_x           = 0
            self.Qp_y           = 0

            self.app_x          = 0.0000e-9
            self.app_y          = 0.0000e-9
            self.app_xy         = 0

            self.alpha     = 0.00192

            self.h1, self.h2       = 4620, 4620*4
            self.V1, self.V2       = 2.e6, 0
            self.dphi1, self.dphi2 = 0, np.pi


            self.longitudinal_focusing = 'non-linear'

            self.gamma = 27.7

            self.p_increment       = 0 * e/c * self.circumference/(self.beta*c)

        else:
            raise ValueError('ERROR: unknown machine configuration', machine_configuration)

        super(SPS, self).__init__(*args, **kwargs)


class HLLHC(synchrotron):

    def __init__(self, *args, **kwargs):

        if 'n_segments' not in kwargs.keys():
            raise ValueError('Number of segments must be specified')

        if 'machine_configuration' not in kwargs.keys():
            raise ValueError('machine_configuration must be specified')

        self.n_segments = kwargs['n_segments']
        self.machine_configuration = kwargs['machine_configuration']

        self.circumference  = 26658.883
        self.s = (np.arange(0, self.n_segments + 1)
                  * self.circumference / self.n_segments)

        if self.machine_configuration == '7TeV':
            self.charge = e
            self.mass = m_p

            self.gamma = np.sqrt( (7000e9*e/(self.mass*c**2))**2 + 1 )

            self.alpha_x = 0 * np.ones(self.n_segments)
            self.beta_x  = 65.9756337546 * np.ones(self.n_segments)
            self.D_x     = 0 * np.ones(self.n_segments)
            self.alpha_y = 0 * np.ones(self.n_segments)
            self.beta_y  = 71.5255058456 * np.ones(self.n_segments)
            self.D_y     = 0 * np.ones(self.n_segments)

            self.Q_x     = 60.31
            self.Q_y     = 60.32

            self.Qp_x    = 0
            self.Qp_y    = 0

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha       = 0.0003225
            self.h1          = 35640
            self.h2          = 71280
            self.V1          = 16e6
            self.V2          = 0
            self.dphi1       = 0
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'

        super(HLLHC, self).__init__(*args, **kwargs)
