# Joint effort of S. PRUNET, P. FOUQUE and B. MAHONEY

import numpy as nm
import mogs
import numpy as np

import pylab as plt
plotnum = 1

# Light speed in m/s
lightspeed = 3e8;
# n takes into account non-Poissonian noise sources, such as fringing
# the original value was 1.5 but 1.0 is the natural value
n=1.0
# Sky level definitions in e-/s/pixel
#dark_sky = {'uS':0.58,'gS':3.34,'rS':4.17,'iS':8.35,'zS':8.35}
dark_sky = {'uS':0.58,'gS':3.34,'rS':4.17,'iS':8.35,'zS':8.35,
            'u' :0.58,'g' :3.34,'r' :4.17,'i' :8.35,'z' :8.35,
            'CaHK':0.22,
# Dustin's number
#'CaHK':0.109,
            'Ha':0.28,'HaOFF':0.28,'OIII':0.22,'OIIIOFF':0.22,'gri':16.0,
            'Y':120,'J':400,'H':2500,'K':1500,'lowoh1':3.2,'lowoh2':4.6,'ch4on':1100,'ch4off':930,'H2':110,'KCont':160,'BrG':110,'W':360,'CO':110}
#grey_sky = {'uS':2.34,'gS':6.68,'rS':6.26,'iS':8.35,'zS':8.35,'Y':150,'J':500,'H':2500,'K':1165,'lowoh1':2.8,'lowoh2':5.0,'ch4on':710.0,'ch4off':973.0,'h2':120.0,'KCont':120.0,'BrG':120.0}
grey_sky = {'uS':2.34,'gS':6.68,'rS':6.26,'iS':8.35,'zS':8.35,
            'u' :2.34,'g' :6.68,'r' :6.26,'i' :8.35,'z' :8.35,
            'CaHK':0.45,'Ha':0.42,'HaOFF':0.42,'OIII':0.45,'OIIIOFF':0.45,'gri':21.0,
            'Y':140,'J':480,'H':2900,'K':1700,'lowoh1':3.8,'lowoh2':5.1,'ch4on':1250,'ch4off':1200,'H2':130,'KCont':170,'BrG':125,'W':410,'CO':130}
#bright_sky = {'uS':9.35,'gS':13.36,'rS':9.39,'iS':8.35,'zS':8.35,'Y':190,'J':750,'H':3250,'K':1500,'lowoh1':5.0,'lowoh2':5.0,'ch4on':1250.0,'ch4off':1250.0,'h2':150.0,'KCont':260.0,'BrG':260.0}
bright_sky = {'uS':9.35,'gS':13.36,'rS':9.39,'iS':8.35,'zS':8.35,
              'u' :9.35,'g' :13.36,'r' :9.39,'i' :8.35,'z' :8.35,
              'CaHK':0.89,'Ha':0.63,'HaOFF':0.63,'OIII':0.89,'OIIIOFF':0.89,'gri':30.,
              'Y':160,'J':560,'H':3400,'K':1900,'lowoh1':4.6,'lowoh2':5.7,'ch4on':1400,'ch4off':1300,'H2':150,'KCont':180,'BrG':140,'W':490,'CO':150}
skies = {'dark' : dark_sky, 'grey' : grey_sky, 'bright' : bright_sky}
# Zero points (in e-/s) for AB magnitudes
#zpts = {'uS':25.68,'gS':26.84,'rS':26.35,'iS':26.14,'zS':25.15,'Y':25.81,'J':26.05,'H':26.65,'K':26.36,'lowoh1':22.40,'lowoh2':22.48,'ch4on':25.31,'ch4off':25.30,'h2':23.75,'KCont':23.71,'BrG':23.71}
#zpts = {'uS':25.68,'gS':26.84,'rS':26.35,'iS':26.14,'zS':25.15,
#        'u' :25.94,'g' :27.19,'r' :26.81,'i' :26.41,'z' :25.05,
#        'CaHK':24.07,'Ha':23.86,'HaOFF':23.90,'OIII':24.19,'OIIIOFF':24.20,
#        'Y':25.87,'J':26.06,'H':26.66,'K':26.38,'lowoh1':23.08,'lowoh2':23.38,'ch4on':25.56,'ch4off':25.47,'h2':23.80,'KCont':23.63,'BrG':23.72}
### Values in adequation with Megacam website, measured in deep fields for all wide bands (2015/06/15) 
zpts = {'uS':25.74,'gS':27.00,'rS':26.50,'iS':26.38,'zS':25.34,
        'u' :25.78,'g' :27.11,'r' :26.74,'i' :26.22,'z' :25.02,
        'CaHK':24.07,'Ha':23.86,'HaOFF':23.90,'OIII':24.19,'OIIIOFF':24.20,'gri':27.95,
        'Y':25.87,'J':26.06,'H':26.66,'K':26.38,'lowoh1':23.08,'lowoh2':23.38,'ch4on':25.56,'ch4off':25.47,'H2':23.80,'KCont':23.63,'BrG':23.72,'W':25.10,'CO':23.76}
# Extinction per airmass variation
d_ext_d_am = {'uS':0.35,'gS':0.15,'rS':0.1,'iS':0.04,'zS':0.03,
              'u' :0.35,'g' :0.15,'r' :0.1,'i' :0.04,'z' :0.03,
              'CaHK':0.15,'Ha':0.1,'HaOFF':0.1,'OIII':0.15,'OIIIOFF':0.15,'gri':0.1,
              'Y':0.02,'J':0.05,'H':0.03,'K':0.05,'lowoh1':0.05,'lowoh2':0.05,'ch4on':0.03,'ch4off':0.03,'H2':0.05,'KCont':0.05,'BrG':0.05,'W':0.04,'CO':0.05}
# Sky background per airmass variation
d_Se_d_am = {'uS':0.34,'gS':0.34,'rS':2.46,'iS':11.6,'zS':11.9,
             'u' :0.34,'g' :0.34,'r' :2.46,'i' :11.6,'z' :11.9,
             'CaHK':0.02,'Ha':0.16,'HaOFF':0.16,'OIII':0.02,'OIIIOFF':0.02,'gri':14.0,  
             'Y':20,'J':40,'H':180,'K':170,'lowoh1':20,'lowoh2':20,'ch4on':42,'ch4off':42,'H2':13,'KCont':23,'BrG':23,'W':30,'CO':23}
# Central filter wavelength and bandwidth in nm (EB 02/08/2023)
# Origin: https://www.cfht.hawaii.edu/Instruments/Imaging/MegaPrime/specsinformation.html
# and https://www.cfht.hawaii.edu/Instruments/Filters/wircam.html
# These are for the filters alone, so for MegaPrime the bandwidths are off
# in the UV and in the NIR because of the triangular shape of the combined response.
# Hence for u and z filters I use instead the FWHMs from the CADC filter curves at
# https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
filter_lambda = {'uS': 375.0, 'gS': 487.0, 'rS': 630.0, 'iS': 770.0, 'zS': 900.0,
                 'u': 355.0, 'g': 475.0, 'r': 640.0, 'i': 776.0, 'z': 925.0,
        'CaHK': 395.1, 'Ha': 659.0, 'HaOFF': 671.8, 'OIII': 500.6, 'OIIIOFF': 510.5, 'gri': 600.0,
        'Y': 1035.0, 'J': 1253.0, 'H': 1631.0, 'K': 2146.0, 'lowoh1': 1061.0, 'lowoh2': 1187.0, 'ch4on': 1690.0, 'ch4off': 1580.0, 'H2': 2122.0, 'KCont': 2218.0, 'BrG': 2166.0, 'W': 1543.0, 'CO': 2320.0}
filter_dlambda = {'uS': 66.0, 'gS': 143.0, 'rS': 124.0, 'iS': 159.0, 'zS': 87.0,
        'u': 52.0, 'g': 154.0, 'r': 148.0, 'i': 155.0, 'z': 65.0,
        'CaHK': 9.6, 'Ha': 10.4, 'HaOFF': 10.7, 'OIII': 9.9, 'OIIIOFF': 9.5, 'gri': 400.0,
        'Y': 100.0, 'J': 88.0, 'H': 99.0, 'K': 98.0, 'lowoh1': 10.0, 'lowoh2': 10.0, 'ch4on': 100.0, 'ch4off': 100.0, 'H2': 32.0, 'KCont': 33.0, 'BrG': 30.0, 'W': 88.4, 'CO': 40.0}

mp_config = {'rpix':0.187,'nccd':5,'gain':1.6,'saturation':65535.0}
#wc_config = {'rpix':0.307,'nccd':30,'gain':3.7,'saturation':32767.0}
wc_config = {'rpix':0.307,'nccd':30,'gain':3.7,'saturation':28000.0}

# Convert flux in erg/s/cm2 to mAB magnitude, assuming the source has a flat
# fnu spectrum within the filter. The wider the filter, the worse this assumption
# EB 02/08/2023
def flux2mag(flux, filter):
    # Frequency range in Hz ~ c * dlambda/(lambda**2) * 1e9 (with lambda in nm)
    dnu = lightspeed * 1e9 * filter_dlambda[filter] / (filter_lambda[filter]**2)
    # Sum over dnu and convert to Jy
    flux /= dnu * 1e-23
    # Use the AB magnitude definition to convert to magnitude
    mag = -2.5 * nm.log10(flux / 3631.0)
    return mag


class snr:

  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,beta=3.0,background='dark',nccd=5.,texp=3600.,rpix=0.187,zpt=None,sky=None):
    
    self.filter=filter
    self.background=background
    self.texp=texp
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.beta=beta #Point source psf profile

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5)
    print(self.Fe)
    print("<br>")
    self.Ftot = self.Fe * self.texp # in e-
    print(self.Ftot)
    print("<br>")
    self.Se = self.sky 
    print(self.Se)
    print("<br>")
    self.Stot = self.Se * self.texp
    print(self.Stot)
    print("<br>")
    # Normalization factor that ensures 2*alpha=FWHM
    self.fac = 2.0**(1./beta)-1.
    return

  def __call__(self):
    return self.SNR(self.Ropt())

  def F(self,R):
    return self.Ftot*(1-(1+self.fac*(R/self.alpha)**2)**(-self.beta+1.))

  def dF(self,R):
    return 2*(self.beta-1.)*self.fac/self.alpha**2 * self.Ftot *R*(1.+self.fac*(R/self.alpha)**2)**(-self.beta)

  def dSNR(self,R):
    df = self.dF(R)
    f = self.F(R)
    return ( df*f + 2*(df*n*nm.pi*R**2 - f*n*nm.pi*R)/self.rpix**2*(self.Stot+self.nccd**2) ) / (2 * (f + n*nm.pi*(R/self.rpix)**2*(self.Stot+self.nccd**2)))**(3./2.)

  def SNR(self,R):
    f = self.F(R)
    return f/nm.sqrt(f + n*nm.pi*(R/self.rpix)**2*(self.Stot+self.nccd**2))

  def frac(self,R):
    return self.F(R)/self.Ftot*100. #in %

  def Ropt(self):
    import scipy.optimize as so
    ##return so.diagbroyden(self.dSNR,self.alpha*1.4,f_tol=1e-6) # Root finder, starting at R=alpha
    return so.brentq(self.dSNR,self.alpha/10.,self.alpha*10.) # Root finder, starting in [alpha/10,10*alpha]. More stable than above

  def R_of_frac(self,fraction):
    # inverse of frac 
    import scipy.optimize as so
    radius=float(so.brentq(lambda x: self.frac(x)-fraction,0.0,1000.*self.alpha))
    return radius

class extsnr:

  ## Computes SNR for extended sources, per arcsec^2 
  def __init__(self,mAB=26.0,fluxormag='mag',filter='r',am=1.2,trans=1.0,background='dark',nccd=5.,texp=3600.,rpix=0.187,zpt=None,sky=None):
    
    self.filter=filter
    self.background=background
    self.texp=texp

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5)
    self.Ftot = self.Fe * self.texp # in e- per arcsec^2
    self.Se = self.sky
    self.Stot = self.Se * self.texp # background in e-  per pixel
    return

  def __call__(self):
    return self.SNR()

  def SNR(self):
    # SNR in 1 arcsec^2. Beware that Stot and nccd are per pixel of area rpix^2
    # while Ftot is a flux per arcsec^2, hence the (1/rpix)^2 factor.
    return self.Ftot/nm.sqrt(self.Ftot + (1.0/self.rpix)**2*(self.Stot+self.nccd**2))

class extexptime:

  ## Computes exposure time for a given SNR for extended sources, per arcsec^2
  def __init__(self,mAB=26.0,fluxormag='mag',filter='r',am=1.2,trans=1.0,background='dark',nccd=5.0,snr=7.0,rpix=0.187,zpt=None,sky=None):

    self.filter=filter
    self.background=background
    self.snr=snr

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)
    
    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5) # in e-/s per arcsec^2
    self.Se = self.sky # in e-/s per pixel

  def __call__(self):
    return self.time()

  def time(self):
    # Solves the quadratic equation root
    delta = (self.snr**4)*(self.Fe +self.Se/self.rpix**2)**2 + 4.0*(self.Fe*self.snr*self.nccd/self.rpix)**2
    tt = (self.snr**2*(self.Fe+self.Se/self.rpix**2) + nm.sqrt(delta)) / (2.0*self.Fe**2)
    return tt


class sattime:

  def __init__(self,mAB=10.,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,beta=3.0,background='dark',rpix=0.187,gain=1.6,fsatADU=65535.0,zpt=None,sky=None):

    self.filter=filter

    self.seeing=seeing
    self.am=am
    self.trans=trans
    self.beta=beta
    self.seeing=seeing
    self.alpha=seeing/2.0
    self.background=background
    self.rpix=rpix
    self.gain=gain #e-/ADU
    self.fsatADU=fsatADU

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt-self.mAB-self.k*(self.am-1))/2.5)
    self.Se = self.sky 
    self.fac = 2.0**(1./beta)-1.
    return

  def ifrac(self,R):
    return (1-(1+self.fac*(R/self.alpha)**2)**(-self.beta+1.))

  def __call__(self):
    FepixADU = self.Fe * self.ifrac(self.rpix/nm.sqrt(nm.pi)) / self.gain
    SeADU = self.Se / self.gain
    tsatobj = self.fsatADU / FepixADU
    tsatsky = self.fsatADU / SeADU
    return tsatobj,tsatsky

class exptime:

  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,beta=3.0,background='dark',nccd=5.,rpix=0.187,snr=10.0,zpt=None,sky=None):
    
    self.snr=snr
    self.filter=filter
    self.background=background
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.beta=beta #Point source psf profile

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0) 

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5) # in e-/s
    self.Se = self.sky # in e-/s
    # Normalization factor that ensures 2*alpha=FWHM
    self.fac = 2.0**(1./beta)-1.

    return

  def __call__(self):
    return self.trmin()

  def ifrac(self,R):
    return (1-(1+self.fac*(R/self.alpha)**2)**(-self.beta+1.))

  def F(self,R,texp): 
    return self.Fe*texp*self.ifrac(R)

  # Exposure time for a given snr (self.snr) as a function of radius
  # Given by one of the roots of a second degree equation in exposure time
  def t(self,R):
    delta = self.snr**4*( self.Fe*self.ifrac(R) + n*nm.pi*(R/self.rpix)**2*self.Se)**2
    delta += 4*self.Fe**2*self.ifrac(R)**2*self.snr**2*n*nm.pi*(R/self.rpix)**2*self.nccd**2
    return ( self.snr**2*(self.Fe*self.ifrac(R) + n*nm.pi*(R/self.rpix)**2*self.Se) + nm.sqrt(delta) )/(2*self.Fe**2*self.ifrac(R)**2)

  def trmin(self):
    ''' returns optimal radius and corresponding minimal exposure time for a given snr '''
    import scipy.optimize as so
    ##ropt = float(so.diagbroyden(lambda x: self.dSNRdR(x,self.t(x)),self.alpha*1.4,f_tol=1e-6)) ## Uses implicit funtion theorem
    ropt = float(so.brentq(lambda x: self.dSNRdR(x,self.t(x)),self.alpha/10.,self.alpha*10.)) ## Uses implicit funtion theorem, more stable than above
    topt = self.t(ropt)
    return (topt,ropt)

  def dFdR(self,R,texp):
    return 2*(self.beta-1.)*self.fac/self.alpha**2 * self.Fe*texp *R*(1.+self.fac*(R/self.alpha)**2)**(-self.beta)

  def dSNRdR(self,R,texp):
    df = self.dFdR(R,texp)
    f = self.F(R,texp)
    return ( df*f + 2*(df*n*nm.pi*R**2 - f*n*nm.pi*R)/self.rpix**2*(self.Se*texp+self.nccd**2) ) / (2 * (f + n*nm.pi*(R/self.rpix)**2*(self.Se*texp+self.nccd**2)))**(3./2.)

  def SNR(self,R,texp):
    f = self.F(R,texp)
    return f/nm.sqrt(f + n*nm.pi*(R/self.rpix)**2*(self.Se*texp+self.nccd**2))

  def Ropt(self):
    import scipy.optimize as so
    ropt = float(so.brentq(lambda x: self.dSNRdR(x,self.t(x)),self.alpha/10.,self.alpha*10.)) ## Uses implicit funtion theorem
    return ropt
    
  def frac(self,R):
    return self.ifrac(R)*100. #in %

  def R_of_frac(self,fraction):
    # inverse of frac 
    import scipy.optimize as so
    radius=float(so.brentq(lambda x: self.frac(x)-fraction,0.0,1000.*self.alpha))
    return radius

class psfsnr:

  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,beta=3.0,background='dark',nccd=5.,texp=3600.,rpix=0.187,gain=1.6,zpt=None,sky=None,xs=0,ys=0):
    
    self.filter=filter
    self.background=background
    self.texp=texp
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.psfaperture=2.5 # radius, in arcsec
    self.beta=beta #Point source psf profile

    self.am=am
    self.trans=trans
    self.gain = gain # in e-/ADU
    self.nccd_ADU=nccd/gain # in ADU/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec
    self.xs=xs # Source coordinates, in arcsec
    self.ys=ys # Source coordinates, in arcsec


    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe_ADU = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5) / self.gain # in ADU/s
    self.Ftot_ADU = self.Fe_ADU * self.texp # in ADU
    self.Se_ADU = self.sky / self.gain # in ADU/s 
    self.Stot_ADU = self.Se_ADU * self.texp

    print('zeropoint: %.3f' % self.zpt)
    print('zeropoint with airmass: %.3f' % (self.zpt - self.k*(self.am-1)))
    print('gain: %.3f' % self.gain)
    print('signal / second [ADU]: %.3f' % self.Fe_ADU)
    print('signal / second [e-]:  %.3f' % (self.Fe_ADU * self.gain))
    print('sky / second [ADU]: %.3f' % self.Se_ADU)
    print('sky / second [e-] : %.3f' % (self.Se_ADU * self.gain))
    
    # Normalization factor that ensures 2*alpha=FWHM
    self.fac = 2.0**(1./beta)-1.
    return

  def modify_texp(self,newtexp):
    self.texp=newtexp
    self.Ftot_ADU = self.Fe_ADU * self.texp
    self.Stot_ADU = self.Se_ADU * self.texp

  def __call__(self):
    return self.SNR()

  def M(self,xi,yi): # Moffat PSF
    norm = (self.beta-1.)*self.fac/(nm.pi*self.alpha**2) * (1. + self.fac*( (xi-self.xs)**2 + (yi-self.ys)**2 )/self.alpha**2)**(-self.beta)
    return norm

  def dMdx(self,xi,yi): #first derivative w.r.t. x
    return 2.*self.beta*(self.beta-1.)*self.fac**2/(nm.pi*self.alpha**4) * (xi-self.xs) * (1. + self.fac*( (xi-self.xs)**2 + (yi-self.ys)**2 )/self.alpha**2)**(-self.beta-1.)

  def dMdy(self,xi,yi): #first derivative w.r.t. x
    return 2.*self.beta*(self.beta-1.)*self.fac**2/(nm.pi*self.alpha**4) * (yi-self.ys) * (1. + self.fac*( (xi-self.xs)**2 + (yi-self.ys)**2 )/self.alpha**2)**(-self.beta-1.)

  def sigma2pix(self,xi,yi):
    sigma2 = self.Stot_ADU/self.gain + self.nccd_ADU**2 + self.Ftot_ADU/self.gain*self.rpix**2*self.M(xi,yi)
    return sigma2
  
  def moment_derivatives(self,xi,yi):
    npix = nm.size(xi)
    moff = self.M(xi,yi).squeeze()

    global plotnum
    plt.clf()
    plt.scatter(xi, yi, c=moff)
    cb = plt.colorbar()
    cb.set_label('Moffat value')
    plt.savefig('plot-%02i.png' % plotnum)
    plotnum += 1

    dmu = nm.ndarray((npix,3),dtype=nm.double)
    ds2 = nm.ndarray((npix,3),dtype=nm.double)
    dmu[:,0] = self.rpix**2 * moff # dmudfs
    dmu[:,1] = self.Ftot_ADU * self.rpix**2 * self.dMdx(xi,yi).squeeze() # dmudxs
    dmu[:,2] = self.Ftot_ADU * self.rpix**2 * self.dMdy(xi,yi).squeeze() # dmudys

    ds2[:,0] = self.rpix**2/self.gain * moff #dsigma2dfs
    ds2[:,1] = dmu[:,1] / self.gain # dsigma2dxs
    ds2[:,2] = dmu[:,2] / self.gain # dsigma2dys

    return (dmu,ds2)

  def Fisher_matrix(self,xi,yi):
    (dmu,ds2) = self.moment_derivatives(xi,yi)
    s2 = self.sigma2pix(xi,yi)
    fish = nm.dot(dmu.T,dmu/s2) + 0.5*nm.dot(ds2.T,ds2/s2**2)

    print('n pix:', len(xi))
    print('Fisher: moment derivs dmu shape:', dmu.shape, 'ds2 shape', ds2.shape)
    print('s2 shape:', s2.shape)
    print('dmu.T dot dmu/s2:', nm.dot(dmu.T, dmu/s2))
    print('0.5 * ds2.T dot ds2/s2**2:', 0.5 * nm.dot(ds2.T,ds2/s2**2))

    print('dmu:', dmu[:10, :])
    
    print('fish:', fish)

    return fish

  def SNR(self):
    global plotnum

    xi,yi = self.pixel_coordinates()
    fish = self.Fisher_matrix(xi,yi)
    print('Pixel coordinates:', len(xi), ', fisher', fish)

    plt.clf()
    plt.plot(xi, yi, '.')
    plt.savefig('plot-%02i.png' % plotnum)
    plotnum += 1

    cov = nm.linalg.inv(fish)
    print('cov', cov)
    sig2flux = nm.sqrt(cov[0,0])


    print('sig2flux:', sig2flux)
    print('sqrt(1/fish[0,0]:', 1./nm.sqrt(fish[0,0]))

    (dmu,ds2) = self.moment_derivatives(xi,yi)
    s2 = self.sigma2pix(xi,yi)
    fish1 = nm.dot(dmu.T,dmu/s2)
    fish2 = 0.5*nm.dot(ds2.T,ds2/s2**2)
    print('fish1', fish1[0,0])
    print('fish2', fish2[0,0])
    f1 = fish1[0,0]
    f2 = fish2[0,0]
    #S = 1./nm.sqrt(f1 + f2)
    f1 = nm.sum(dmu[:,0] * dmu[:,0] / s2[:,0])
    print('Approx f1:', f1)

    moff = self.M(xi,yi).squeeze()
    dmu0 = self.rpix**2 * moff # dmudfs
    print('dmu0 shape', dmu0.shape)
    f1 = nm.sum(dmu0 * dmu0 / s2[:,0])
    print('Approx f1:', f1)

    print('s2: Stot_ADU/gain:', (self.Stot_ADU/self.gain))
    print('s2: nccd_ADU**2:', (self.nccd_ADU**2))
    fterm = (self.Ftot_ADU/self.gain*self.rpix**2*self.M(xi,yi))
    print('s2: Ftot_ADU/gain*rpix**2*M: range', fterm.min(), 'to', fterm.max())
    # approx s2:
    s2x = self.Stot_ADU/self.gain
    f1 = nm.sum(dmu0 * dmu0 / s2x)
    S = 1./nm.sqrt(f1)
    print('Approx f1(w/approx s2):', f1)
    print('Approx: S=', S)

    print('sum of moff:', nm.sum(moff))
    # weird, their 'moff' has a sum ~= 1./pixscale**2 !!

    #print('dmu0 range:', dmu0.min(), dmu0.max())
    #print('Moffat NEA:', 1./(self.rpix**2 * nm.sqrt(nm.sum(moff ** 2))))

    print('sum of dmu0:', nm.sum(dmu0))
    
    umoff = moff * self.rpix**2
    print('Moffat NEA:', 1./nm.sum(umoff**2))

    S = nm.sqrt(s2x) / (self.rpix**2 * nm.sqrt(nm.sum(moff ** 2)))
    print('Approx: S=', S)

    print('Approx SNR:', self.Ftot_ADU/S)

    # Ftot_ADU = Fe_ADU * texp
    # Stot_ADU = Se_ADU * texp
    #    Se_ADU = self.sky / self.gain
    # Stot_ADU = self.sky * texp / gain

    # s2x = self.sky * texp / gain**2
    print('s2x:', s2x)
    s2x = self.sky * self.texp / self.gain**2
    print('s2x:', s2x)

    print('real SNR:', self.Ftot_ADU/S)

    print('PSF norm:', nm.sqrt(nm.sum(dmu0 ** 2)))

    npixlin = int(nm.ceil(self.psfaperture/self.rpix))
    xpix,ypix = np.meshgrid(np.arange(-npixlin, npixlin+1),
                            np.arange(-npixlin, npixlin+1))
    umoff = self.M(xpix * self.rpix, ypix * self.rpix)
    umoff /= np.sum(umoff)
    print('PSF (umoff) norm:', nm.sqrt(nm.sum(umoff ** 2)))

    # fwhm in arcsec
    seeing = self.alpha * 2.0
    # fwhm in pixels
    seeing_pix = seeing / self.rpix
    # psf sigma in pixels
    sig_pix = seeing_pix / 2.35
    #xpix = xi / self.rpix
    #ypix = yi / self.rpix
    G = nm.exp(-0.5 * (xpix**2 + ypix**2) / sig_pix**2)
    G /= nm.sum(G)
    print('PSF norm for Gaussian:', nm.sqrt(nm.sum(G ** 2)))
    print('Gaussian NEA:', 1./nm.sum(G**2))

    plt.clf()
    mx = max(umoff.max(), G.max())
    ima = dict(interpolation='nearest', origin='lower', vmin=0, vmax=mx)
    plt.subplot(2,2,1)
    plt.imshow(umoff, **ima)
    plt.title('Moffat')
    plt.subplot(2,2,2)
    plt.imshow(G, **ima)
    plt.title('Gaussian')
    ima = dict(interpolation='nearest', origin='lower', vmin=np.log10(mx) - 3, vmax=np.log10(mx))
    plt.subplot(2,2,3)
    plt.imshow(np.log10(np.maximum(umoff, 1e-10)), **ima)
    plt.title('Moffat (log)')
    plt.subplot(2,2,4)
    plt.imshow(np.log10(np.maximum(G, 1e-10)), **ima)
    plt.title('Gaussian (log)')
    plt.savefig('plot-%02i.png' % plotnum)
    plotnum += 1

    #snr = self.Fe_ADU * self.texp * (self.rpix**2 * nm.sqrt(nm.sum(moff ** 2))) / nm.sqrt(s2x)
    #snr = self.Fe_ADU * self.texp * nm.sqrt(nm.sum(dmu0 ** 2)) / nm.sqrt(s2x)
    #snr = self.Fe_ADU * self.texp * nm.sqrt(nm.sum(dmu0 ** 2)) / nm.sqrt(self.sky * self.texp) * self.gain
    snr = self.Fe_ADU * self.gain * self.texp * nm.sqrt(nm.sum(dmu0 ** 2)) / nm.sqrt(self.sky * self.texp)
    snr_g = self.Fe_ADU * self.gain * self.texp * nm.sqrt(nm.sum(G ** 2)) / nm.sqrt(self.sky * self.texp)
    print('snr:', snr)
    print('snr(Gaussian):', snr_g)

    # approx s2y: including "nccd" term (~dark current??)
    # the units on this make *no* sense
    s2y = self.Stot_ADU/self.gain + self.nccd_ADU**2
    # Stot_ADU = self.sky * texp / gain
    print('s2y:', s2y)
    #snr = self.Fe_ADU * self.texp * nm.sqrt(nm.sum(dmu0 ** 2)) / nm.sqrt(s2y)
    snr_2 = self.Fe_ADU * self.texp * nm.sqrt(nm.sum(dmu0 ** 2)) / nm.sqrt(self.sky * self.texp / self.gain**2 + self.nccd_ADU**2)
    print('snr_2:', snr_2)
    snr_2g = self.Fe_ADU * self.texp * nm.sqrt(nm.sum(G ** 2)) / nm.sqrt(self.sky * self.texp / self.gain**2 + self.nccd_ADU**2)

    f1 = nm.sum(dmu0 * dmu0) / s2x
    S = 1./nm.sqrt(f1)
    
    print('t_exp', self.texp, 'Ftot_ADU:', self.Ftot_ADU, 'sig2flux:', sig2flux, '-> SNR', self.Ftot_ADU/sig2flux)
    #return self.Ftot_ADU/sig2flux
    #return self.Ftot_ADU/S
    #return snr_g
    return snr_2

  def pixel_coordinates(self):
    # First compute the number of pixels around the source center, within a large aperture
    npixlin = int(nm.ceil(self.psfaperture/self.rpix))
    tmp = nm.arange(-npixlin,npixlin+1)*self.rpix
    tmp.shape = ((2*npixlin+1,1))
    ones = nm.ones(2*npixlin+1)
    ones.shape = ((2*npixlin+1,1))
    xs_int = nm.floor(self.xs/self.rpix)*self.rpix # Coordinates of pixel where source center falls
    ys_int = nm.floor(self.ys/self.rpix)*self.rpix
    xi = nm.kron(tmp+xs_int,ones.T) # 2D maps of pixel center coordinates, around xs_int
    yi = (nm.kron(tmp+ys_int,ones.T)).T
    print('pixel_coordinates: rpix=%f, npixlin=%i, psfaperture=%f' % (self.rpix, npixlin, self.psfaperture))
    dist = ( (xi-self.xs)**2+(yi-self.ys)**2 ) < self.psfaperture**2
    ou = nm.argwhere(dist.flatten())
    xi = xi.flatten()[ou]
    yi = yi.flatten()[ou]
    return (xi,yi)

class psfexptime:

  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,beta=3.0,background='dark',nccd=5.,rpix=0.187,gain=1.6,snr=10.0,zpt=None,sky=None):
    
    self.snr=snr
    self.filter=filter
    self.background=background
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.beta=beta #Point source psf profile

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec
    self.gain=gain
    self.snr=snr
#    self.zpt=zpt
#    self.sky=sky

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0) 
    print('sky (including airmass term):', self.sky)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
      print('Flux2mag:', self.mAB)
      self.fluxormag = 'mag'
    else:
      self.mAB=mAB
        
    # Compute bracket search, based on aperture photometry exposure time
    tt=exptime(self.mAB,self.fluxormag,self.filter,self.am,self.trans,self.seeing,self.beta,self.background,self.nccd,self.rpix,self.snr,zpt=self.zpt,sky=self.sky)
    texp_aper = tt()[0]
    self.tmin = texp_aper/3.
    self.tmax = texp_aper*3.
    print('Aperture photom exposure time:', texp_aper)
    # Initiate psfsnr instance
    self.ps=psfsnr(self.mAB,self.fluxormag,self.filter,self.am,self.trans,self.seeing,self.beta,self.background,self.nccd,texp_aper,self.rpix,self.gain,zpt=self.zpt,sky=self.sky)
    return

  def __call__(self):
    return self.exptime_compute()

  def snrdiff(self,texp):
    self.ps.modify_texp(texp)
    return self.ps.SNR()-self.snr

  def exptime_compute(self):
    import scipy.optimize as so
    return so.brentq(self.snrdiff,self.tmin,self.tmax)

class galsnr:

  ## Galactic profile
  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,half_light_radius=1.0,sersic_index=4,background='dark',nccd=5.,texp=3600.,rpix=0.187,zpt=None,sky=None): 
    

    self.filter=filter
    self.background=background
    self.texp=texp
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.beta=3.0 # Fixed for now
    self.half_light_radius=half_light_radius # Galaxy R(1/2)
    self.sersic_index=sersic_index

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0)

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5)
    self.Ftot = self.Fe * self.texp # in e-
    self.Se = self.sky 
    self.Stot = self.Se * self.texp
    # Instantiate Mixture of Gaussian class
    self.mg=mogs.mogs(half_light_radius=self.half_light_radius,sersic_n=self.sersic_index,moffat_alpha=self.alpha)
    ## Define sum of R(1/2) and alpha
    self.Reff = (self.alpha+self.half_light_radius)

    return

  def __call__(self):
    return self.SNR(self.Ropt())

  def F(self,R):
    return self.Ftot*self.mg.convolution_mog_cumul(R)

  def dF(self,R):
    return self.Ftot * 2.*nm.pi*R * self.mg.convolution_mog(R)

  def dSNR(self,R):
    df = self.dF(R)
    f = self.F(R)
    return ( df*f + 2*(df*n*nm.pi*R**2 - f*n*nm.pi*R)/self.rpix**2*(self.Stot+self.nccd**2) ) / (2 * (f + n*nm.pi*(R/self.rpix)**2*(self.Stot+self.nccd**2)))**(3./2.)

  def SNR(self,R):
    f = self.F(R)
    return f/nm.sqrt(f + n*nm.pi*(R/self.rpix)**2*(self.Stot+self.nccd**2))

  def frac(self,R):
    return self.F(R)/self.Ftot*100. #in %

  def R_of_frac(self,fraction):
    # inverse of frac 
    import scipy.optimize as so
    radius=float(so.brentq(lambda x: self.frac(x)-fraction,0.0,1000.*self.Reff))
    return radius

  def Ropt(self):
    import scipy.optimize as so
    return so.brentq(self.dSNR,self.Reff/10.,self.Reff*10.) # Root finder, starting in [Reff/10,10*Reff]

class galexptime:

  ## Galactic profile
  def __init__(self,mAB=24.9,fluxormag='mag',filter='r',am=1.2,trans=1.0,seeing=0.69,half_light_radius=1.0,sersic_index=4,background='dark',nccd=5.,rpix=0.187,snr=10.0,zpt=None,sky=None):
    
    self.snr=snr
    self.filter=filter
    self.background=background
    self.seeing=seeing
    self.alpha=seeing/2.0 #HWHM, in arcsec
    self.beta=3.0 #Point source psf profile, fixed for now
    self.half_light_radius=half_light_radius # galaxy R(1/2)
    self.sersic_index=sersic_index # galaxy Sersic profile index

    self.am=am
    self.trans=trans
    self.nccd=nccd # in e-/pixel/sqrt(second)
    self.rpix=rpix # Linear size of pixel in arcsec

    if (zpt is not None):
      self.zpt=zpt
    else:
      self.zpt=zpts[filter]
    self.k=d_ext_d_am[filter]
    if (sky is not None):
      self.sky=sky
    else:
      self.sky=skies[background][filter] + d_Se_d_am[filter] * (self.am-1.0) 

    self.fluxormag=fluxormag
    if (fluxormag=='flux'):
      self.mAB=flux2mag(mAB, filter)
    else:
      self.mAB=mAB
        
    self.Fe = self.trans * 10.0**((self.zpt - self.mAB - self.k*(self.am-1))/2.5) # in e-/s
    self.Se = self.sky # in e-/s
    # Instantiate Mixture of Gaussian class
    self.mg=mogs.mogs(half_light_radius=self.half_light_radius,sersic_n=self.sersic_index,moffat_alpha=self.alpha)
    ## Define sum of R(1/2) and alpha
    self.Reff = (self.alpha+self.half_light_radius)

    return

  def __call__(self):
    return self.trmin()

  def ifrac(self,R):
    return self.mg.convolution_mog_cumul(R)

  def F(self,R,texp): 
    return self.Fe*texp*self.ifrac(R)

  # Exposure time for a given snr (self.snr) as a function of radius
  # Given by one of the roots of a second degree equation in exposure time
  def t(self,R):
    delta = self.snr**4*( self.Fe*self.ifrac(R) + n*nm.pi*(R/self.rpix)**2*self.Se)**2
    delta += 4*self.Fe**2*self.ifrac(R)**2*self.snr**2*n*nm.pi*(R/self.rpix)**2*self.nccd**2
    return ( self.snr**2*(self.Fe*self.ifrac(R) + n*nm.pi*(R/self.rpix)**2*self.Se) + nm.sqrt(delta) )/(2*self.Fe**2*self.ifrac(R)**2)

  def trmin(self):
    ''' returns optimal radius and corresponding minimal exposure time for a given snr '''
    import scipy.optimize as so
    ##ropt = float(so.diagbroyden(lambda x: self.dSNRdR(x,self.t(x)),self.alpha*1.4,f_tol=1e-6)) ## Uses implicit funtion theorem
    ropt = float(so.brentq(lambda x: self.dSNRdR(x,self.t(x)),self.Reff/10.,self.Reff*10.)) ## Uses implicit funtion theorem, more stable than above
    topt = self.t(ropt)
    return (topt,ropt)

  def Ropt(self):
    import scipy.optimize as so
    ropt = float(so.brentq(lambda x: self.dSNRdR(x,self.t(x)),self.Reff/10.,self.Reff*10.)) ## Uses implicit funtion theorem
    return ropt
    
  def dFdR(self,R,texp):
    return self.Fe*texp * 2.*nm.pi*R * self.mg.convolution_mog(R)

  def dSNRdR(self,R,texp):
    df = self.dFdR(R,texp)
    f = self.F(R,texp)
    return ( df*f + 2*(df*n*nm.pi*R**2 - f*n*nm.pi*R)/self.rpix**2*(self.Se*texp+self.nccd**2) ) / (2 * (f + n*nm.pi*(R/self.rpix)**2*(self.Se*texp+self.nccd**2)))**(3./2.)

  def SNR(self,R,texp):
    f = self.F(R,texp)
    return f/nm.sqrt(f + n*nm.pi*(R/self.rpix)**2*(self.Se*texp+self.nccd**2))

  def frac(self,R):
    return self.ifrac(R)*100. #in %

  def R_of_frac(self,fraction):
    # inverse of frac 
    import scipy.optimize as so
    radius=float(so.brentq(lambda x: self.frac(x)-fraction,0.0,1000.*self.Reff))
    return radius

