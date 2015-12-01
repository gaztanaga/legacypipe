from __future__ import print_function
import os
import numpy as np

from legacypipe.image import LegacySurveyImage
from legacypipe.common import *

import fitsio
from astrometry.util.file import trymakedirs
from astrometry.util.fits import fits_table
from astrometry.util.util import Tan, Sip, anwcs_t
from tractor.tractortime import TAITime

'''
Code specific to images from the (intermediate) Palomar Transient Factory (iPTF/PTF), bands = g,R.
11 CCDs and 1.2m telescope at Palomar Observatory.
'''

#### test code
class TestCode(object):
    def __init__(self):
        pass

    def read_image(self,imgfn,**kwargs):
        '''return numpy array of pixels given filename'''
        return fitsio.read(imgfn, ext=0, header=True) 

    def read_dq(self,dqfn,**kwargs):
        '''return bit mask which Tractor calls "data quality" image
        PTF DMASK BIT DEFINITIONS
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
        INFOBITS=                    0 / Database infobits (2^2 and 2^3 excluded)
        '''
        dq= fitsio.read(dqfn, ext=0, header=False)
        return dq.astype(np.int16)
        
    def read_invvar(self,imgfn,dqfn, clip=False, clipThresh=0.2, **kwargs):
        print('*** No Weight Map *** computing invvar with image and data quality mapd')
        dq=read_dq(dqfn) 
        img,hdr=read_image(imgfn,header=True)
        assert(dq.shape == img.shape)
        invvar=np.zeros(img.shape)
        invvar[dq == 0]= np.power(img[dq == 0],-0.5)
        invvar[dq == 2]= np.power(img[dq == 2],-0.5) #SExtractor ojbect if binary(00010) = 2
        if clip:
            # Clamp near-zero (incl negative!) invvars to zero.
            # These arise due to fpack.
            if clipThresh > 0.:
                med = np.median(invvar[invvar > 0])
                thresh = clipThresh * med
            else:
                thresh = 0.
            invvar[invvar < thresh] = 0
        return invvar
####

class PtfImage(LegacySurveyImage):
    '''
   
    A LegacySurveyImage subclass to handle images from the Dark Energy
    Camera, DECam, on the Blanco telescope.

    '''
    def __init__(self, decals, t):
        super(PtfImage, self).__init__(decals, t)
        self.dqfn= os.path.join(os.path.dirname(self.imgfn),'../mask',os.path.basename(self.imgfn))
        self.dqfn = self.dqfn.replace('_scie_', '_mask_')
        #self.wtfn = self.imgfn.replace('_ooi_', '_oow_')

        self.name= self.imgfn
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
        #calibdir = os.path.join(self.decals.get_calib_dir(), self.camera)
        #self.pvwcsfn = os.path.join(calibdir, 'astrom-pv', self.calname + '.wcs.fits')
        #self.sefn = os.path.join(calibdir, 'sextractor', self.calname + '.fits')
        #self.psffn = os.path.join(calibdir, 'psfex', self.calname + '.fits')
        #self.skyfn = os.path.join(calibdir, 'sky', self.calname + '.fits')

    #def get_image_shape(self):
    #    return self.height, self.width

    #def shape(self):
    #    return self.get_image_shape()

    #def get_tractor_image(self, **kwargs):
    #    tim = super(PtfImage, self).get_tractor_image(**kwargs)
    #    return tim

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
       
    def read_image(self,**kwargs):
        '''return numpy array of pixels given filename'''
        print('Reading image from', self.imgfn, 'hdu', self.hdu)
        return fitsio.read(self.imgfn, ext=self.hdu, header=True) 

    def read_dq(self,**kwargs):
        '''return bit mask which Tractor calls "data quality" image
        PTF DMASK BIT DEFINITIONS
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
        INFOBITS=                    0 / Database infobits (2^2 and 2^3 excluded)
        '''
        print('Reading data quality image from', self.dqfn, 'hdu', self.hdu)
        dq= fitsio.read(self.dqfn, ext=self.hdu, header=False)
        return dq.astype(np.int16)
   
    def read_invvar(self, clip=False, clipThresh=0.2, **kwargs):
        print('*** No Weight Map *** computing invvar with image and data quality mapd')
        dq=self.read_dq() 
        img,hdr=self.read_image(header=True)
        assert(dq.shape == img.shape)
        invvar=np.zeros(img.shape)
        invvar[dq == 0]= np.power(img[dq == 0],-0.5)
        invvar[dq == 2]= np.power(img[dq == 2],-0.5) #SExtractor ojbect if binary(00010) = 2
        if clip:
            # Clamp near-zero (incl negative!) invvars to zero.
            # These arise due to fpack.
            if clipThresh > 0.:
                med = np.median(invvar[invvar > 0])
                thresh = clipThresh * med
            else:
                thresh = 0.
            invvar[invvar < thresh] = 0
        return invvar

    def read_sky_model(self, **kwargs):
        print('Constant sky model, median of ', self.imgfn)
        img,hdr = self.read_image(header=True)
        sky = np.median(img)
        print('Median "sky" =', sky)
        sky = ConstantSky(sky)
        sky.version = '0'
        sky.plver = '0'
        return sky

    def run_calibs(self, **kwargs):
        # def run_calibs(self, pvastrom=True, psfex=True, sky=True, se=False,
        #           funpack=False, fcopy=False, use_mask=True,
        #           force=False, just_check=False, git_version=None,
        #           splinesky=False):
        '''
        Run calibration pre-processing steps.
        '''
        print('doing Nothing for Calibrations')
        pass

    def get_tractor_image(self, slc=None, radecpoly=None,
                          gaussPsf=False, const2psf=False, pixPsf=False,
                          splinesky=False,
                          nanomaggies=True, subsky=True, tiny=5,
                          dq=True, invvar=True, pixels=True):
        '''
        Returns a tractor.Image ("tim") object for this image.
        
        Options describing a subimage to return:

        - *slc*: y,x slice objects
        - *radecpoly*: numpy array, shape (N,2), RA,Dec polygon describing bounding box to select.

        Options determining the PSF model to use:

        - *gaussPsf*: single circular Gaussian PSF based on header FWHM value.
        - *const2Psf*: 2-component general Gaussian fit to PsfEx model at image center.
        - *pixPsf*: pixelized PsfEx model at image center.

        Options determining the sky model to use:
        
        - *splinesky*: median filter chunks of the image, then spline those.

        Options determining the units of the image:

        - *nanomaggies*: convert the image to be in units of NanoMaggies;
          *tim.zpscale* contains the scale value the image was divided by.

        - *subsky*: instantiate and subtract the initial sky model,
          leaving a constant zero sky model?

        '''
        from astrometry.util.miscutils import clip_polygon
        print('#### IN KBJs get tractor image code! ###')
        get_dq = dq
        get_invvar = invvar
        
        band = self.band
        imh,imw = self.get_image_shape()
        wcs = self.get_wcs()
        x0,y0 = 0,0
        x1 = x0 + imw
        y1 = y0 + imh
        #if don't comment out tim = NoneType b/c clips all pixels out
        #if slc is None and radecpoly is not None:
        #    imgpoly = [(1,1),(1,imh),(imw,imh),(imw,1)]
        #    ok,tx,ty = wcs.radec2pixelxy(radecpoly[:-1,0], radecpoly[:-1,1])
        #    tpoly = zip(tx,ty)
        #    clip = clip_polygon(imgpoly, tpoly)
        #    clip = np.array(clip)
        #    if len(clip) == 0:
        #        return None
        #    x0,y0 = np.floor(clip.min(axis=0)).astype(int)
        #    x1,y1 = np.ceil (clip.max(axis=0)).astype(int)
        #    slc = slice(y0,y1+1), slice(x0,x1+1)
        #    if y1 - y0 < tiny or x1 - x0 < tiny:
        #        print('Skipping tiny subimage')
        #        return None
        #if slc is not None:
        #    sy,sx = slc
        #    y0,y1 = sy.start, sy.stop
        #    x0,x1 = sx.start, sx.stop

        #old_extent = (x0,x1,y0,y1)
        #new_extent = self.get_good_image_slice((x0,x1,y0,y1), get_extent=True)
        #if new_extent != old_extent:
        #    x0,x1,y0,y1 = new_extent
        #    print('Applying good subregion of CCD: slice is', x0,x1,y0,y1)
        #    if x0 >= x1 or y0 >= y1:
        #        return None
        #    slc = slice(y0,y1), slice(x0,x1)
        if pixels:
            print('Reading image slice:', slc)
            img,imghdr = self.read_image(header=True, slice=slc)
            #print('SATURATE is', imghdr.get('SATURATE', None))
            #print('Max value in image is', img.max())
            # check consistency... something of a DR1 hangover
            #e = imghdr['EXTNAME']
            #assert(e.strip() == self.ccdname.strip())
        else:
            img = np.zeros((imh, imw))
            imghdr = dict()
            if slc is not None:
                img = img[slc]
            
        if get_invvar:
            invvar = self.read_invvar(slice=slc, clipThresh=0.)
        else:
            invvar = np.ones_like(img)
            
        if get_dq:
            dq = self.read_dq(slice=slc)
            invvar[dq != 0] = 0.
        if np.all(invvar == 0.):
            print('Skipping zero-invvar image')
            return None
        assert(np.all(np.isfinite(img)))
        assert(np.all(np.isfinite(invvar)))
        assert(not(np.all(invvar == 0.)))

        # header 'FWHM' is in pixels
        # imghdr['FWHM']
        psf_fwhm = self.fwhm 
        psf_sigma = psf_fwhm / 2.35
        primhdr = self.read_image_primary_header()

        sky = self.read_sky_model(splinesky=splinesky, slc=slc)
        midsky = 0.
        if subsky:
            print('Instantiating and subtracting sky model...')
            from tractor.sky import ConstantSky
            skymod = np.zeros_like(img)
            sky.addTo(skymod)
            img -= skymod
            midsky = np.median(skymod)
            zsky = ConstantSky(0.)
            zsky.version = sky.version
            zsky.plver = sky.plver
            del skymod
            del sky
            sky = zsky
            del zsky

        magzp = self.decals.get_zeropoint_for(self)
        if isinstance(magzp,str): 
            print('WARNING: no ZeroPoint in header for image: ',self.imgfn)
            magzp= 23.
        orig_zpscale = zpscale = NanoMaggies.zeropointToScale(magzp)
        if nanomaggies:
            # Scale images to Nanomaggies
            img /= zpscale
            invvar *= zpscale**2
            if not subsky:
                sky.scale(1./zpscale)
            zpscale = 1.

        assert(np.sum(invvar > 0) > 0)
        if get_invvar:
            sig1 = 1./np.sqrt(np.median(invvar[invvar > 0]))
        else:
            # Estimate from the image?
            # # Estimate per-pixel noise via Blanton's 5-pixel MAD
            slice1 = (slice(0,-5,10),slice(0,-5,10))
            slice2 = (slice(5,None,10),slice(5,None,10))
            mad = np.median(np.abs(img[slice1] - img[slice2]).ravel())
            sig1 = 1.4826 * mad / np.sqrt(2.)
            print('sig1 estimate:', sig1)
            invvar *= (1. / sig1**2)
            
        assert(np.all(np.isfinite(img)))
        assert(np.all(np.isfinite(invvar)))
        assert(np.isfinite(sig1))

        if subsky:
            ##
            imgmed = np.median(img[invvar>0])
            if np.abs(imgmed) > sig1:
                print('WARNING: image median', imgmed, 'is more than 1 sigma away from zero!')
                # Boom!
                assert(False)

        twcs = ConstantFitsWcs(wcs)
        if x0 or y0:
            twcs.setX0Y0(x0,y0)

        #print('gaussPsf:', gaussPsf, 'pixPsf:', pixPsf, 'const2psf:', const2psf)
        psf = self.read_psf_model(x0, y0, gaussPsf=gaussPsf, pixPsf=pixPsf,
                                  const2psf=const2psf, psf_sigma=psf_sigma)

        tim = Image(img, invvar=invvar, wcs=twcs, psf=psf,
                    photocal=LinearPhotoCal(zpscale, band=band),
                    sky=sky, name=self.name + ' ' + band)
        assert(np.all(np.isfinite(tim.getInvError())))

        # PSF norm
        psfnorm = self.psf_norm(tim)
        print('PSF norm', psfnorm, 'vs Gaussian',
              1./(2. * np.sqrt(np.pi) * psf_sigma))

        # Galaxy-detection norm
        tim.band = band
        galnorm = self.galaxy_norm(tim)
        print('Galaxy norm:', galnorm)
        
        # CP (DECam) images include DATE-OBS and MJD-OBS, in UTC.
        import astropy.time
        #mjd_utc = mjd=primhdr.get('MJD-OBS', 0)
        mjd_tai = astropy.time.Time(primhdr['DATE-OBS']).tai.mjd
        tim.slice = slc
        tim.time = TAITime(None, mjd=mjd_tai)
        tim.zr = [-3. * sig1, 10. * sig1]
        tim.zpscale = orig_zpscale
        tim.midsky = midsky
        tim.sig1 = sig1
        tim.psf_fwhm = psf_fwhm
        tim.psf_sigma = psf_sigma
        tim.propid = self.propid
        tim.psfnorm = psfnorm
        tim.galnorm = galnorm
        tim.sip_wcs = wcs
        tim.x0,tim.y0 = int(x0),int(y0)
        tim.imobj = self
        tim.primhdr = primhdr
        tim.hdr = imghdr
        tim.plver = str(primhdr['PTFVERSN']).strip()
        tim.skyver = (sky.version, sky.plver)
        tim.wcsver = ('-1','-1') #wcs.version, wcs.plver)
        tim.psfver = (psf.version, psf.plver)
        if get_dq:
            tim.dq = dq
        tim.dq_bits = CP_DQ_BITS
        tim.saturation = imghdr.get('SATURATE', None)
        tim.satval = tim.saturation or 0.
        if subsky:
            tim.satval -= midsky
        if nanomaggies:
            tim.satval /= orig_zpscale
        subh,subw = tim.shape
        tim.subwcs = tim.sip_wcs.get_subimage(tim.x0, tim.y0, subw, subh)
        mn,mx = tim.zr
        tim.ima = dict(interpolation='nearest', origin='lower', cmap='gray',
                       vmin=mn, vmax=mx)
        return tim


class PtfDecals(Decals):
    def __init__(self, **kwargs):
        super(PtfDecals, self).__init__(**kwargs)
        self.image_typemap.update({'ptf' : PtfImage})

    def get_zeropoint_for(self,tractor_image):
        print('WARNING: zeropoints from header of ',tractor_image.imgfn)
        hdr=fitsio.read_header(tractor_image.imgfn)
        return hdr['IMAGEZPT']
    
    #def ccds_touching_wcs(self, wcs, **kwargs):
    #    '''PTF testing, continue even if no overlap with DECaLS bricks
    #    '''
    #    print('WARNING: ccds do not have to be overlapping with DECaLS bricks')
    #    T = self.get_ccds_readonly()
    #    I = ccds_touching_wcs(wcs, T, **kwargs)
    #    if len(I) == 0:
    #        return None
    #    T = T[I]
    #    return T

    def photometric_ccds(self, CCD):
        '''PTF testing process non-photometric ccds too'''
        print('WARNING: non-photometric ccds allowed')
        good = np.ones(len(CCD), bool)
        return np.flatnonzero(good)
    
    def get_ccds(self):
        '''
        Return SMALL CCD for testing: 2 rows
        '''
        fn = os.path.join(self.decals_dir, 'decals-ccds.fits')
        if not os.path.exists(fn):
            fn += '.gz'
        print('Reading CCDs from', fn)
        T = fits_table(fn)
        print('Got', len(T), 'CCDs')
        if 'ccdname' in T.columns():
            # "N4 " -> "N4"
            T.ccdname = np.array([s.strip() for s in T.ccdname])
        T= T[ [np.where(T.filter == 'R')[0][0],np.where(T.filter == 'g')[0][0]] ] #1 R and 1 g band
        print('ccd ra= ',T.ra,'ccd dec= ',T.dec) 
        return T


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
    
    kwargs['bands'] = 'gR'

    # runbrick...
    run_brick(opt.brick, **kwargs)
    return 0
    
if __name__ == '__main__':
    import sys
    sys.exit(main())
