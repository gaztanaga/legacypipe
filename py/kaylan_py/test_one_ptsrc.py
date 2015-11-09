from __future__ import print_function
import os
import numpy as np

import pylab as plt

import astropy.time
import fitsio

from legacypipe.common import *
from legacypipe.image import LegacySurveyImage
from legacypipe.runbrick import run_brick
from astrometry.util.util import Tan

from tractor.psf import GaussianMixturePSF
from tractor.tractortime import TAITime

class FakeBrick(object):
    def __init__(self):
        self.brickid = -1
        self.brickname = 'fakey'
        self.ra = 0.
        self.dec = 0.
        self.ra1  = -0.1
        self.ra2  =  0.1
        self.dec1 = -0.1
        self.dec2 =  0.1
        
    def about(self):
        print('I am a FakeBrick')
    pass

#class TestingDecals(Decals):
class TestingDecals(object):
    def __init__(self, **kwargs):
        self.decals_dir = None
        
    def get_brick_by_name(self, brickname):
        return FakeBrick()

    def ccds_touching_wcs(self, wcs, **kwargs):
        ccds = fits_table()
        ccds.xxx = [0.]
        ccds.filter = ['r']
        ccds.propid = ['FundMe']
        ccds.to_np_arrays()
        return ccds

    def photometric_ccds(self, ccds):
        return np.arange(len(ccds))

    def get_image_object(self, ccd):
        return TestingImage(ccd)

    def drop_cache(self):
        pass
    
class TestingImage(LegacySurveyImage):
    def __init__(self, ccd):
        self.band = ccd.filter
        self.exptime = 1000000
        self.fwhm = 3.0
        self.name = 'TestingImage'
        self.ccd = ccd
        
    def run_calibs(self, **kwargs):
        pass

    def get_tractor_image(self, **kwargs):
        imagew, imageh = 100,100
        psf_sigma = self.fwhm / 2.35
        pixscale= 0.262/3600.
        crval1, crval2, crpix1, crpix2= 0., 0., imagew/2., imageh/2.
        cd11, cd12, cd21, cd22= -pixscale,0.,0.,pixscale
        tanwcs = Tan(crval1, crval2, crpix1, crpix2, cd11, cd12, cd21, cd22, imagew, imageh)
        sig1 = 0.01

        image = np.zeros((imageh, imagew), np.float32)
        xx,yy = np.meshgrid(np.arange(imagew), np.arange(imageh))
        cx,cy = imagew/2., imageh/2.
        flux = 1.
        image += flux * 1./(2.*np.pi * psf_sigma**2) * np.exp(-0.5 * ((xx - cx)**2 + (yy - cy)**2) / psf_sigma**2)
        image += np.random.normal(size=image.shape) * sig1

        plt.clf()
        plt.imshow(image)
        plt.colorbar()
        plt.savefig('image.png')
        
        tim = Image(data=image,
                    inverr=np.ones((imageh, imagew), np.float32) / sig1,
                    name=self.name,
                    psf=GaussianMixturePSF(1., 0.,0.,
                                           psf_sigma**2, psf_sigma**2, 0.),
                    photocal=LinearPhotoCal(1., band=self.band),
                    wcs=ConstantFitsWcs(tanwcs),
                    )
        tim.skyver = ('1', '1')
        tim.psfver = ('1', '1')
        tim.wcsver = ('1', '1')
        tim.band = self.band
        tim.psf_fwhm = self.fwhm
        tim.psf_sigma = psf_sigma
        tim.psfnorm = self.psf_norm(tim)
        tim.galnorm = self.galaxy_norm(tim)
        tim.plver = '1'
        tim.subwcs = tanwcs
        tim.sig1 = sig1
        tim.dq = np.zeros((imageh,imagew), np.uint8)
        tim.dq_bits = dict(satur=1)
        tim.primhdr = fitsio.FITSHDR()
        tim.x0 = 0
        tim.y0 = 0
        tim.propid = self.ccd.propid
        tim.imobj = self
        self.psfnorm = tim.psfnorm
        self.galnorm = tim.galnorm
        mjd_tai = astropy.time.Time('2000-01-01T00:00:00.000000').tai.mjd
        tim.time = TAITime(None, mjd=mjd_tai)
        return tim

# render testing image in memory

fakedecals = TestingDecals()
run_brick(None, decals=fakedecals, blacklist=False,
          writePickles=False, forceAll=True,
          sdssInit=False, wise=False,
          width=200, height=200, outdir='testing',
          coadd_bw=True)

