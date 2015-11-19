from __future__ import print_function
import os
import numpy as np

from legacypipe.image import LegacySurveyImage
from legacypipe.common import *

import fitsio
from astrometry.util.file import trymakedirs
from astrometry.util.fits import fits_table
from astrometry.util.util import Tan, Sip, anwcs_t

'''
Code specific to images from the (intermediate) Palomar Transient Factory (iPTF/PTF), bands = g,R.
11 CCDs and 1.2m telescope at Palomar Observatory.
'''


class PtfImage(LegacySurveyImage):
    '''
   
    A LegacySurveyImage subclass to handle images from the Dark Energy
    Camera, DECam, on the Blanco telescope.

    '''
    def __init__(self, decals, t):
        super(PtfImage, self).__init__(decals, t)
        if self.band == 'g':
            ind=self.imgfn.rfind('images/')
            beg=self.imgfn[:ind+len('images/')]
            end=self.imgfn[ind+len('images/'):]
            mid='ptf/ptf_22694/gband/science/'
            self.imgfn=beg+mid+end 
        else: 
            print('band != g, band= ',self.band)
            raise ValueError
        self.name= self.imgfn
        print('KJB: dir(self)= ',dir(self))
        print('KJB: dir(decals)= ',dir(decals))
        #for i in dir(self):
        #    if i.startswith('__'): continue
        #    else: print('self.%s= ' % i,getattr(self, i))

        #self.dqfn = self.imgfn.replace('_ooi_', '_ood_')
        #self.wtfn = self.imgfn.replace('_ooi_', '_oow_')

        #for attr in ['imgfn', 'dqfn', 'wtfn']:
        #    fn = getattr(self, attr)
        #    if os.path.exists(fn):
        #        continue
        #    if fn.endswith('.fz'):
        #        fun = fn[:-3]
        #        if os.path.exists(fun):
        #            print('Using      ', fun)
        #            print('rather than', fn)
        #            setattr(self, attr, fun)
        print('KJB: calib_dir= ',self.decals.get_calib_dir(),'camera=',self.camera)
        #calibdir = os.path.join(self.decals.get_calib_dir(), self.camera)
        #self.pvwcsfn = os.path.join(calibdir, 'astrom-pv', self.calname + '.wcs.fits')
        #self.sefn = os.path.join(calibdir, 'sextractor', self.calname + '.fits')
        #self.psffn = os.path.join(calibdir, 'psfex', self.calname + '.fits')
        #self.skyfn = os.path.join(calibdir, 'sky', self.calname + '.fits')

    def get_image_shape(self):
        return self.height, self.width

    def shape(self):
        return self.get_image_shape()

    def get_tractor_image(self, **kwargs):
        tim = super(PtfImage, self).get_tractor_image(**kwargs)
        return tim

    def __str__(self):
        return 'PTF ' + self.name
    
    #override funcs get_tractor_image calls
    def get_wcs(self):
        return self.read_pv_wcs()

    def read_pv_wcs(self):
        '''extract wcs from fits header directly'''
        hdr = fitsio.read_header(self.imgfn, self.hdu)
        H,W = self.get_image_shape()
        wcs= Tan(hdr['CRVAL1'], hdr['CRVAL2'],hdr['CRPIX1'],hdr['CRPIX2'],\
                     hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2'],\
                     float(W),float(H))
        return wcs
    #    wcs.version = '0' #done in bok.py 
    #    wcs.plver = '0'
    #    return wcs
        #from astrometry.util.util import Sip
        #print('Reading WCS from', self.pvwcsfn)
        #wcs = Sip(self.pvwcsfn)
        #dra,ddec = self.decals.get_astrometric_zeropoint_for(self)
        #r,d = wcs.get_crval()
        #print('Applying astrometric zeropoint:', (dra,ddec))
        #wcs.set_crval((r + dra, d + ddec))
        #hdr = fitsio.read_header(self.pvwcsfn)
        #wcs.version = hdr.get('LEGPIPEV', '')
        #if len(wcs.version) == 0:
        #    wcs.version = hdr.get('TRACTORV', '').strip()
        #    if len(wcs.version) == 0:
        #        wcs.version = str(os.stat(self.pvwcsfn).st_mtime)
        #wcs.plver = hdr.get('PLVER', '').strip()
        #return wcs

    def get_good_image_subregion(self):
        pass
       
    def read_image(self,fname):
        '''return numpy array of pixels given filename'''
        F=fitsio.FITS(fname)
        return F[0].read() 

    def read_dq(self, fname, header=False, **kwargs):
        '''return "data quality" (flags) image'''
        from distutils.version import StrictVersion
        #print('Reading data quality from', self.dqfn, 'hdu', self.hdu)
        #dq,hdr = self._read_fits(self.dqfn, self.hdu, header=True, **kwargs)
        F=fitsio.FITS(fname)
        dq=F[0].read() 
        '''
        0 = good, 2 = detected object, all others = bad
        in mask image i see values 0,2,512 not 1,3-15!?
        BIT00   =                    0 / AIRCRAFT/SATELLITE TRACK
        BIT01   =                    1 / OBJECT (detected by SExtractor)
        BIT02   =                    2 / HIGH DARK-CURRENT
        BIT03   =                    3 / RESERVED FOR FUTURE USE
        BIT04   =                    4 / NOISY
        BIT05   =                    5 / GHOST
        BIT06   =                    6 / CCD BLEED
        BIT07   =                    7 / RAD HIT
        BIT08   =                    8 / SATURATED
        BIT09   =                    9 / DEAD/BAD
        BIT10   =                   10 / NAN (not a number)
        BIT11   =                   11 / DIRTY (10-sigma below coarse local median)
        BIT12   =                   12 / HALO
        BIT13   =                   13 / RESERVED FOR FUTURE USE
        BIT14   =                   14 / RESERVED FOR FUTURE USE
        BIT15   =                   15 / RESERVED FOR FUTURE USE
        #1 = bad
        #2 = no value (for remapped and stacked data)
        #3 = saturated
        #4 = bleed mask
        #5 = cosmic ray
        #6 = low weight
        #7 = diff detect (multi-exposure difference detection from median)
        #8 = long streak (e.g. satellite trail)
        '''
        #dqbits = np.zeros(dq.shape, np.int16)
        #dqbits[dq == 1] |= CP_DQ_BITS['badpix']
        #dqbits[dq == 2] |= CP_DQ_BITS['badpix']
        #dqbits[dq == 3] |= CP_DQ_BITS['satur']
        #dqbits[dq == 4] |= CP_DQ_BITS['bleed']
        #dqbits[dq == 5] |= CP_DQ_BITS['cr']
        #dqbits[dq == 6] |= CP_DQ_BITS['badpix']
        #dqbits[dq == 7] |= CP_DQ_BITS['trans']
        #dqbits[dq == 8] |= CP_DQ_BITS['trans']

        #dq = dqbits
        return dq

    def read_invvar(self, f_mask,f_img,clip=True, **kwargs):
        #print('Reading inverse-variance from', self.wtfn, 'hdu', self.hdu)
        #invvar = self._read_fits(self.wtfn, self.hdu, **kwargs)
        mask=self.read_dq(f_mask) 
        img=self.read_image(f_img)
        invvar=np.zeros(img.shape)
        igood= np.where(np.logical_or(mask == 2,mask == 0) == True)
        #ibad= np.where(np.logical_and(mask != 2,mask != 0) == True)
        invvar[igood]= np.power(img[igood],-0.5) 
        #if clip:
        #    # Clamp near-zero (incl negative!) invvars to zero.
        #    # These arise due to fpack.
        #    med = np.median(invvar[invvar > 0])
        #    thresh = 0.2 * med
        #    invvar[invvar < thresh] = 0
        return invvar

    def read_sky_model(self, **kwargs):
        #from bok.py
        img = self.read_image()
        sky = np.median(img)
        print('Median "sky" model:', sky)
        sky = ConstantSky(sky)
        sky.version = '0'
        sky.plver = '0'
        return sky

    def run_calibs(self):
        '''
        Run calibration pre-processing steps.
        '''
        pass

class PtfDecals(Decals):
    def __init__(self, **kwargs):
        super(PtfDecals, self).__init__(**kwargs)
        self.image_typemap.update({'ptf' : PtfImage})
 
def main():
    from runbrick import run_brick, get_parser, get_runbrick_kwargs
    
    parser = get_parser()
    opt = parser.parse_args()
    if opt.brick is None and opt.radec is None:
        parser.print_help()
        return -1

    print('Forcing --no-blacklist')
    opt.blacklist = False

    kwargs = get_runbrick_kwargs(opt)
    if kwargs in [-1,0]:
        return kwargs

    decals = PtfDecals(decals_dir=opt.decals_dir)
    kwargs['decals'] = decals
    
    # runbrick...
    run_brick(opt.brick, **kwargs)
    return 0
    
if __name__ == '__main__':
    import sys
    sys.exit(main())
