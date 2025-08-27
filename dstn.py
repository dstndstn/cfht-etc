# Compute etime,snr:5,sourcetype:pointsource,filter:g,mAB:26,photometry:psf,aperture:optimalaperture,sky:dark 0,airmass:1.2,transparency:1,seeing:0.7,etime: 638.09,radius:,flux frac:,obj sat:,sky sat:undefined

# sourcetype:pointsource
# photometry:psf
#	keywords['gain']=diet.mp_config['gain']
#		tt=diet.psfexptime(**keywords)
#		print ('PSF Etime:   %.2f<br>'%tt())
#		ss=diet.psfsnr(**keywords)
#		print ('PSF SNR:   %.2f<br>'%ss())


import diet

#mAB=27.15,
#fluxormag='mag',

                     #filter='CaHK',


if False:
    tt = diet.psfexptime(gain=diet.mp_config['gain'],
                         snr=5.0,
                         filter='CaHK',
                         mAB=1.0e-16,
                         fluxormag='flux',
                         background='dark',
                         am=1.0,
                         trans=1.0,
                         seeing=0.8)

'''
3147172p.fits.fz
QSOGRADE=                    1 / QSO grade (1=good 5=unusable)
IQMEASUR=                 0.79 / Measured IQ in arcsec
BGMEASUR=                 3.64 / Measured skybg in ADU/pix/sec on chip

AIRMASS =                1.351 / Airmass at start of observation

PHOT_C  =              26.4770 / Elixir zero point - measured for camera run
PHOT_CS =               0.0087 / Elixir zero point - scatter
PHOT_NS =                15535 / Elixir zero point - N stars
PHOT_NM =                   11 / Elixir zero point - N images
PHOT_C0 =              26.6600 / Elixir zero point - nominal
PHOT_X  =              -0.0870 / Elixir zero point - color term (slope)
PHOT_DX =               0.0000 / Elixir zero point - mean color term
PHOT_K  =              -0.1500 / Elixir zero point - airmass term
PHOT_C1 = 'g_SNLS            ' / Elixir zero point - color 1
PHOT_C2 = 'r_SNLS            ' / Elixir zero point - color 2
COMMENT   Formula for Photometry, based on keywords given in this header:
COMMENT   m = -2.5*log(DN) + 2.5*log(EXPTIME)
COMMENT   M = m + PHOT_C + PHOT_K*(AIRMASS - 1)
COMMENT         + PHOT_X*(PHOT_C1 - PHOT_C2 - PHOT_DX)

'''
#mag = 25.0
#airmass = 1.351
#seeing = 0.79

mag = 25.0
airmass = 1.0
seeing = 0.8
#gain=diet.mp_config['gain']
gain = 1.5

print('mag:', mag)
tt = diet.psfexptime(snr=5.0,
                     filter='g',
                     mAB=mag,
                     fluxormag='mag',
                     background='dark',
                     am=airmass,
                     trans=1.0,
                     gain=gain,
                     seeing=seeing)

#exptime = 120.
exptime = 144.86
print('Exptime:', exptime)
tt.ps.modify_texp(exptime)
snr = tt.ps.SNR()
print('SNR:', snr)

# -> creates a "psfsnr" object.
#   - calls modify_texp(), SNR() to find exptime.

tt = tt()
print('Exptime', tt)
