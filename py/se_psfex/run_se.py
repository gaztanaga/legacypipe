import os
import numpy as np
import fitsio
import argparse
import glob

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-img_search",action="store",help='search string for processed images')
parser.add_argument("-configdir",action="store",help='directory with SExtractor config files')
parser.add_argument("-savefn",action="store",help='output file name')
args = parser.parse_args()

# grab header values...
imgfns= glob.glob(os.path.join('/project/projectdirs/desi/imaging/data/ptf/cosmos/pimage',args.img_search))
print('found %d files' % len(imgfns))
if len(imgfns) == 0: raise ValueError
for cnt,imgfn in enumerate(imgfns):
    print('reading %d of %d images: %s' % (cnt,len(imgfns),imgfn))
    hdr=fitsio.read_header(imgfn,ext=0)
    magzp  = hdr['IMAGEZPT']
    seeing = hdr['PIXSCALE'] * hdr['MEDFWHM']
    gain= hdr['GAIN']
    print('magzp= ',magzp,'seeing= ',seeing,'gain= ',gain)
    #read and save bit mask with 2 subtracted off
    maskfn= os.path.join(os.path.dirname(imgfn).replace('pimage','mask'),\
                            os.path.basename(imgfn).replace('_scie_','_mask_'))
    mask,mhd= fitsio.read(maskfn,ext=0,header=True)
    mask[mask > 0]= mask[mask > 0]-2
    new_maskfn= os.path.join(os.path.dirname(args.savefn),os.path.basename(maskfn))
    fitsio.write(new_maskfn,mask,header=mhd)
    #
    cmd = ' '.join(['sex','-c', os.path.join(args.configdir, 'DECaLS.se'),
                    '-FLAG_IMAGE %s' % new_maskfn,
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
                    '-CATALOG_NAME', args.savefn,
                    '-GAIN %f' % gain,
                    imgfn])
    print(cmd)
    if os.system(cmd):
        raise RuntimeError('Command failed: ' + cmd)

