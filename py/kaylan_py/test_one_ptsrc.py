from __future__ import print_function
import os
import numpy as np

from legacypipe.common import *
from legacypipe.image import LegacySurveyImage
from legacypipe.runbrick import run_brick

class FakeBrick(object):
    def __init__(self):
        self.brickid = -1
        self.brickname = 'fakey'
        self.ra = 0.
        self.dec = 0.
        self.ra1 = 0.
        self.ra2 = 0.
        self.dec1 = 0.
        self.dec2 = 0.
        
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
    
class TestingImage(object): #LegacySurveyImage):
    def __init__(self, ccd):
        self.band = ccd.filter
        self.exptime = 1000000

    def run_calibs(self, **kwargs):
        pass

    def get_tractor_image(self, **kwargs):
        tim = Image(np.zeros((100,100)), inverr=np.ones((100,100)))
        tim.skyver = ('1', '1')
        tim.psfver = ('1', '1')
        tim.wcsver = ('1', '1')
        tim.band = self.band
        tim.plver = '1'
		tim.subwcs = #wcs_for_brick(self.get_brick_by_name(), W=3600, H=3600, pixscale=0.262)
        return tim

# render testing image in memory

fakedecals = TestingDecals()
run_brick(None, decals=fakedecals, blacklist=False)

