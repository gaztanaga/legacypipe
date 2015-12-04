import optparse
import fitsio
import os
import numpy as np
import glob
from astrometry.util.fits import fits_table
from astrometry.util.starutil_numpy import degrees_between, hmsstring2ra, dmsstring2dec
from astrometry.util.util import Tan, Sip, anwcs_t
import argparse

def exposure_metadata(filenames, hdus=None, trim=None):
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

#parser = optparse.OptionParser(epilog=ep)
#parser.add_option('--splinesky', action='store_true', default=False,
#                  help='Use flexible sky model?')
#
#parser.add_option('--coadd-bw', action='store_true', default=False,
#                  help='Create grayscale coadds if only one band is available?')
#
#parser.add_option('--on-bricks', action='store_true', default=False,
#                  help='Tractor-on-bricks?')
#
#print()
#print('runbrick.py starting at', datetime.datetime.now().isoformat())
#print('Command-line args:', sys.argv)
#print()
#
#opt,args = parser.parse_args()
#

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="test")
    parser.add_argument("-pimage_search",action="store",help='full path + search string for ptf images')
    args = parser.parse_args()
    
    files=glob.glob(args.pimage_search)
    T=exposure_metadata(files)
#save table to .fits
    outfn = 'decals-ccds.fits'
    T.writeto(outfn)
    print('Wrote', outfn)
