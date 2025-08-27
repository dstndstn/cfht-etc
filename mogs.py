import numpy as nm
from scipy.special import gammaincinv, gamma

class mogs:

  def __init__(self,half_light_radius=1.0,sersic_n=4,moffat_alpha=1.0,pickledir='./mog-pickles'):

    # initialize moffat mixture amplitudes and variances

    # Read mog params
    self.moffat_K         = 6
    self.moffat_pickle    = pickledir+'/moffat_K%02d_MR08.pickle'%self.moffat_K
    moffat_pars           = nm.load(self.moffat_pickle)
    self.moffat_amplitudes = moffat_pars[0:self.moffat_K]
    self.moffat_variances  = moffat_pars[self.moffat_K:self.moffat_K*2]

    # Scaling
    self.moffat_alpha      = moffat_alpha
    self.moffat_variances *= moffat_alpha**2
    # amplitudes are automatically taken care of since each gaussian is normalized 

    # initialize Sersic amplitudes and variances
    self.sersic_K          = 12
    self.sersic_n          = sersic_n
    if (sersic_n==1): # Exponential profile
      self.sersic_pickle='exp_K%02d_MR08.pickle'%self.sersic_K
    elif (sersic_n==2):
      self.sersic_pickle='ser2_K%02d_MR08.pickle'%self.sersic_K
    elif (sersic_n==3):
      self.sersic_pickle='ser3_K%02d_MR08.pickle'%self.sersic_K
    elif (sersic_n==4):
      self.sersic_pickle='dev_K%02d_MR08.pickle'%self.sersic_K
    elif (sersic_n==5):
      self.sersic_pickle='ser5_K%02d_MR08.pickle'%self.sersic_K
    else:
      print('Unsupported Sersic profile')
      return

    self.sersic_alpha = gammaincinv(2.*self.sersic_n, 0.5)
    self.sersic_norm  = nm.exp(self.sersic_alpha)*2.*nm.pi*self.sersic_n*gamma(2*self.sersic_n)/self.sersic_alpha**(2*self.sersic_n)
    self.sersic_pickle = pickledir+'/'+self.sersic_pickle
    self.sersic_pars = nm.load(self.sersic_pickle)
    self.sersic_amplitudes = self.sersic_pars[0:self.sersic_K]
    # Normalise integral of sersic mog to 1
    # The next line assumes that the fit is good enough for a theoretical normalization
    # unfortunately for high enough sersic_n, this becomes untrue, with normalization errors as big as 5%
    # for a deVaucouleurs profile. Hence we will just enforce integral of mog to be 1
    # self.sersic_amplitudes /= self.sersic_norm
    self.sersic_amplitudes /= nm.sum(self.sersic_amplitudes)
    # Will do the same for the Moffat, so that convolution is precisely normalized to 1
    self.moffat_amplitudes /= nm.sum(self.moffat_amplitudes)
    #
    self.sersic_variances = self.sersic_pars[self.sersic_K:self.sersic_K*2]

    # Scale them by half_light_radius
    self.half_light_radius = half_light_radius
    self.sersic_variances *= half_light_radius**2

    # Convolution
    self.convolution_K = self.moffat_K * self.sersic_K
    # amplitudes get multiplied
    self.convolution_amplitudes=nm.outer(self.moffat_amplitudes,self.sersic_amplitudes).flatten()
    # variances get added
    self.convolution_variances=nm.add.outer(self.moffat_variances,self.sersic_variances).flatten()

    return

  def moffat_mog(self,x):

    """ MoG approximation to (normalized) Moffat """
    return self.mog(x,self.moffat_amplitudes,self.moffat_variances,self.moffat_K)

  def moffat_mog_cumul(self,x):

    """ MoG approximation to (normalized) cumulative (disc) integral of Moffat """
    return self.mog_cumul(x,self.moffat_amplitudes,self.moffat_variances,self.moffat_K)

  def moffat(self,x):
    
    """
    Moffat PSF profile, normalized to 1 in 2-d integration
    """
    alpha=self.moffat_alpha
    beta=3.0
    betafac = 2.0**(1./beta)-1.0
    prefac = (beta-1.)*betafac/(nm.pi*alpha**2)
    return prefac * (1.0+betafac*(x/alpha)**2)**(-beta)
  
  def moffat_cumul(self,x):
    """
    Moffat cumulative (disc) integral
    """
    alpha=self.moffat_alpha
    beta=3.0
    betafac = 2.0**(1./beta)-1.0
    return 1. - (1. + betafac*(x/alpha)**2)**(1.-beta)

  def sersic_mog(self,x):

    """ MoG approximation of Sersic profile """
    return self.mog(x,self.sersic_amplitudes,self.sersic_variances,self.sersic_K)

  def sersic_mog_cumul(self,x):

    """ MoG approximation of cumulative (disc) integral of Sersic profile """
    return self.mog_cumul(x,self.sersic_amplitudes,self.sersic_variances,self.sersic_K)

  def sersic(self,x):

    return nm.exp( -self.sersic_alpha*(x**(1./self.sersic_n) -1.0) ) / self.sersic_norm

  def __call__(self,x):

    """ Convolution of Moffat and Sersic MoG approximations """
    return self.convolution_mog(x)

  def convolution_mog(self,x):

    """ Convolution of Moffat and Sersic MoG approximations """
    return self.mog(x,self.convolution_amplitudes,self.convolution_variances,self.convolution_K)
   
  def convolution_mog_cumul(self,x):

    """ MoG approximation of cumulative (disc) integral of convolution of Moffat and Sersic profiles """
    return self.mog_cumul(x,self.convolution_amplitudes,self.convolution_variances,self.convolution_K)

  def mog(self,x,amps,vars,K):

    ''' Computes corresponding mogs at (vector) values x'''
    y=0.
    for i in range(K):
      y += amps[i]*self.not_normal(x,vars[i])

    return y

  def mog_cumul(self,x,amps,vars,K):
    
    ''' Cumulative (disc) integral of mog '''
    y=0.
    for i in range(K):
      y += amps[i]*self.not_normal_cumul(x,vars[i])
    return y

  def not_normal(self,x, V):
    
    """
    Make a one-dimensional profile of a two-dimensional Gaussian.

    Note strange normalization because this is for 2-d Gaussians
    (but only ever called in 1-d).
    """
    exparg = -0.5 * x * x / V
    if (nm.isscalar(x)):
      result = 1. / (2.*nm.pi*V)*nm.exp(exparg)
    else:
      result = nm.zeros_like(x)
      I = ((exparg > -1000) * (exparg < 1000))
      result[I] = 1. / (2. * nm.pi * V) * nm.exp(exparg[I])
    return result

  def not_normal_cumul(self,x,V):
    
    """ 2D Cumulative (disc) integral of 2D gaussian """
    return 1.0-nm.exp(-0.5*x*x/V)


