'''PTF analysis flow:
1. legacypipe/cp_files.py -- copy all images want to analyze into path/to/images/
2. this script -- make 1 directory in path/to/ for each _scie_,_mask_ pair with that dir's ccd.fits file only containing that one file, script also puts decals_dir, soft links, psfex, etc in each dir
3. run tractor on all dir's pointing it to the decals_dir,output dir each time'''


import sys
import os
import numpy as np
import fitsio
import argparse
import glob
from astrometry.util.fits import fits_table
from astrometry.util.starutil_numpy import degrees_between, hmsstring2ra, dmsstring2dec
from astrometry.util.util import Tan, Sip, anwcs_t

from legacypipe.runptf import read_image,read_dq,read_invvar,ptf_zeropoint

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-images_dir",action="store",help='path/to/images/')
#parser.add_argument("-start_stop",nargs=2,type=int,action="store",help='start and stop index, no 0th index, so if 100 image files then [1,100] does all of them, [2,5] does the the 2nd through the 5th, ccd.fits table made for all fits files BUT SExtractor/PSFex only run on index [start,stop] of that fits file list')
parser.add_argument("-configdir",action="store",help='directory with SExtractor config files')
parser.add_argument("-bricks_table",action="store",help='all bricks table, not dr2 verions')
args = parser.parse_args()

def make_dir(name):
    if not os.path.exists(name): 
        os.makedirs(name)
    else: print 'WARNING path exists: ',name

def ptf_exposure_metadata(filenames, hdus=None, trim=None):
    nan = np.nan
    primkeys = [('FILTER',''),
                ('RA', nan),
                ('DEC', nan),
                ('AIRMASS', nan),
                ('DATE-OBS', ''),
                ('EXPTIME', nan),
                ('EXPNUM', 0),
                ('MJD-OBS', 0),
                ('PROPID', ''),
                ]
    hdrkeys = [('AVSKY', nan),
               ('ARAWGAIN', nan),
               ('FWHM', nan),
               ('CRPIX1',nan),
               ('CRPIX2',nan),
               ('CRVAL1',nan),
               ('CRVAL2',nan),
               ('CD1_1',nan),
               ('CD1_2',nan),
               ('CD2_1',nan),
               ('CD2_2',nan),
               ('EXTNAME',''),
               ('CCDNAME',''),
               ('CCDNUM',''),
               ('CAMERA',''),
               ]

    otherkeys = [('IMAGE_FILENAME',''), ('IMAGE_HDU',0),
                 ('HEIGHT',0),('WIDTH',0),
                 ]

    allkeys = primkeys + hdrkeys + otherkeys
    #for each hdu and file, append ptf header info to vals
    vals = dict([(k,[]) for k,d in allkeys])
    for i,fn in enumerate(filenames):
        print('Reading', (i+1), 'of', len(filenames), ':', fn)
        F = fitsio.FITS(fn)
        cpfn = fn
        if trim is not None:
            cpfn = cpfn.replace(trim, '')

        if hdus is not None:
            hdulist = hdus
        else:
            hdulist = range(0, len(F))
        #loop through hdus
        for hdu in hdulist:
            for k,d in allkeys: #d is not used below
                if k in F[hdu].read_header().keys():
                    vals[k].append(F[hdu].read_header()[k])
                else: 
                    continue #will be set below 
            #special handling
            H,W = F[hdu].get_info()['dims'] 
            vals['HEIGHT'].append(H)
            vals['WIDTH'].append(W)
            vals['IMAGE_HDU'].append(hdu)
            vals['AVSKY'].append(0.) #estimate of sky level, ok to set to 0 in legacypipe
            vals['EXTNAME'].append(os.path.basename(fn)[-8:-5]) #CCD id, ex) N16 for DECam
            vals['CCDNAME'][-1]= os.path.basename(fn)[-8:-5]
            vals['FILTER'][-1]= vals['FILTER'][-1].strip().lower()  #unique ptf band name
            vals['EXPNUM'].append(0) #default value givein in function "exposure_metadata"
            vals['CAMERA'].append('ptf') #default value givein in function "exposure_metadata"
            vals['IMAGE_FILENAME'].append(os.path.basename(fn))
            #diff naming convenctions between (DECaLS,PTF)
            map=[
                    ('RA', 'OBJRA'),
                    ('DEC', 'OBJDEC'),
                    ('MJD-OBS', 'OBSMJD'),
                    ('PROPID', 'PTFPID'),
                   ('ARAWGAIN', 'GAIN'),
                   ('FWHM', 'MEDFWHM'),
                   ('CCDNUM','CCDID'),
                    ]
            for decal,ptf in map: 
                try: vals[decal].append(F[hdu].read_header()[ptf])
                except AttributeError: print "WARNING, could not find ",decal,"or",ptf
            #for k,d in allkeys: print "FINAL vals: k= ",k,"vals[k]= ",vals[k]
    #header info now stord in val dict
    T = fits_table()
    for k,d in allkeys:
        T.set(k.lower().replace('-','_'), np.array(vals[k]))


    #finish naming conventions
    T.filter = np.array([s.split()[0] for s in T.filter])
#KJB LEFT OFF HERE
    T.ra_bore  = np.array([hmsstring2ra (s) for s in T.ra ])
    T.dec_bore = np.array([dmsstring2dec(s) for s in T.dec])

    T.ra  = np.zeros(len(T))
    T.dec = np.zeros(len(T))
    for i in range(len(T)):
        W,H = T.width[i], T.height[i]

        wcs = Tan(T.crval1[i], T.crval2[i], T.crpix1[i], T.crpix2[i],
                  T.cd1_1[i], T.cd1_2[i], T.cd2_1[i], T.cd2_2[i], float(W), float(H))
        
        xc,yc = W/2.+0.5, H/2.+0.5
        rc,dc = wcs.pixelxy2radec(xc,yc)
        T.ra [i] = rc
        T.dec[i] = dc
    #anything with type int64 crashes fits.write()
    T.expnum= T.expnum.astype(np.int16)
    T.image_hdu= T.image_hdu.astype(np.int16)
    T.height= T.height.astype(np.int16)
    T.width= T.width.astype(np.int16)
    #for c in T.columns(): 
        #if type(T.get(c)[0]) is np.int64: T.c= T.c.astype(np.int16) #answer='yes'
        #else:answer='no'
        #print('column= ',c,"type= ",type(T.get(c)),"type is int64?",answer) 
    #sanity check  
    for c in T.columns(): print (c,T.get(c))
    return T


scie_files= glob.glob(os.path.join(args.images_dir,'PTF*_scie_*.fits'))
if len(scie_files) == 0: raise ValueError
for full_fn in scie_files:
    print 'making dir structure for:\n'
    print os.path.basename(full_fn)
    #make dir structure
    fn= os.path.basename(full_fn)
    ustr,field= fn[34:44],fn[49:56]
    rootdir=os.path.join(args.images_dir,'../',field+'_'+ustr)
    make_dir(rootdir)
    #decal_dir
    calib_dir= 'decals_dir/calib/ptf/'
    for dn in [calib_dir+'se_config',calib_dir+'se',calib_dir+'psfex']: make_dir(os.path.join(rootdir,dn))
    cmd= ' '.join(['cp',args.bricks_table,os.path.join(rootdir,'decals_dir/')])
    if os.system(cmd): raise ValueError
    #soft link to images/
    home=os.getcwd()
    os.chdir(os.path.join(rootdir,'decals_dir/'))
    if not os.path.exists('images'):
        cmd= ' '.join(['ln','-s','../../images/','images'])
        if os.system(cmd): raise ValueError
    os.chdir(home)
    print('made soft links')
    #make ccd table
    ccd_fname= os.path.join(rootdir,'decals_dir/','decals-ccds.fits')
    if not os.path.exists(ccd_fname): 
        T=ptf_exposure_metadata([full_fn])
        T.writeto(ccd_fname)
        print('Wrote ccd fits table, now SExtractor + PSFex')
    #SExtractor, PSFex 
    make_dir('junk') #need temp dir for mask-2 and invvar map
    #[start,stop] files only
    #subset= scie_files[args.start_stop[0]-1:args.start_stop[1]]
    for cnt,imgfn in enumerate([full_fn]):
        print('SExtractor,PSFex on %d of %d images' % (cnt,len([fn])))
        #SExtractor
        hdu=0
        maskfn= imgfn.replace('_scie_','_mask_')
        print('imgfn= ',imgfn)
        invvar= read_invvar(imgfn,maskfn,hdu) #note, all post processing on image,mask done in read_invvar
        mask= read_dq(maskfn,hdu)
        maskfn= os.path.join('junk',os.path.basename(maskfn))
        invvarfn= maskfn.replace('_mask_','_invvar_')
        fitsio.write(maskfn, mask)
        fitsio.write(invvarfn, invvar)
        print('wrote mask-2 to %s, invvar to %s' % (maskfn,invvarfn))
        #run se
        magzp  = ptf_zeropoint(imgfn)
        hdr=fitsio.read_header(imgfn,ext=hdu)
        seeing = hdr['PIXSCALE'] * hdr['MEDFWHM']
        gain= hdr['GAIN']
        sefn= os.path.join(rootdir,calib_dir,'se/',os.path.basename(imgfn).replace('.fits','.se_cat'))
        cmd = ' '.join(['sex','-c', os.path.join(args.configdir, 'DECaLS.se'),
                        '-WEIGHT_IMAGE %s' % invvarfn, '-WEIGHT_TYPE MAP_WEIGHT',
                        '-GAIN %f' % gain,
                        '-FLAG_IMAGE %s' % maskfn,
                        '-FLAG_TYPE OR',
                        '-SEEING_FWHM %f' % seeing,
                        '-DETECT_MINAREA 3',
                        '-PARAMETERS_NAME', os.path.join(args.configdir, 'DECaLS.param'),
                        '-FILTER_NAME', os.path.join(args.configdir, 'gauss_3.0_5x5.conv'),
                        '-STARNNW_NAME', os.path.join(args.configdir, 'default.nnw'),
                        '-PIXEL_SCALE 0',
                        # SE has a *bizarre* notion of "sigma"
                        '-DETECT_THRESH 1.0',
                        '-ANALYSIS_THRESH 1.0',
                        '-MAG_ZEROPOINT %f' % magzp,
                        '-CATALOG_NAME', sefn,
                        imgfn])
        if os.system(cmd):
            raise RuntimeError('Command failed: ' + cmd)
        #PSFex
        cmd= ' '.join(['psfex',sefn,'-c', os.path.join(args.configdir,'DECaLS.psfex'),
                        '-PSF_DIR',os.path.join(rootdir,calib_dir,'psfex/')])
        if os.system(cmd):
            raise RuntimeError('Command failed: ' + cmd)                   
        #delete temporary mask-2 and invvar
        os.remove(maskfn)
        os.remove(invvarfn)
    print('finished dir: %s' % rootdir)
print('done')
