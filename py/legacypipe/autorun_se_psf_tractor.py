import os
import numpy as np
import fitsio
import argparse
import glob

from legacypipe.runptf import read_image,read_dq,read_invvar,ptf_zeropoint

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-img_search",action="store",help='full path + search string for ptf images')
parser.add_argument("-configdir",action="store",help='directory with SExtractor config files')
parser.add_argument("-se_outdir",action="store",help='set to run SExtractor, output dir',required=False)
parser.add_argument("-psfx_outdir",action="store",help='set to run PSFex, output dir', required=False)
parser.add_argument("-radec",nargs=2,type=float,action="store",help='for legacypipe')
parser.add_argument("-decals_dir",action="store",help='for legacypipe')
parser.add_argument("-tractor_outdir",action="store",help='set to run legacypipe/tractor', required=False)
parser.add_argument("-wh",nargs=2,type=int,action="store",help='for legacypipe')
parser.add_argument("-coadd_bw",action="store",help='for legacypipe, any value calls it',required=False)
args = parser.parse_args()

#make temporary dir 'junk' to write mask-2 and invvar (weight map) images for SExtractor
if not os.path.exists('junk'): os.makedirs('junk') 
# grab header values...
imgfns= glob.glob(args.img_search)
print('found %d files' % len(imgfns))
if len(imgfns) == 0: raise ValueError
for cnt,imgfn in enumerate(imgfns):
    print('reading %d of %d images: %s' % (cnt,len(imgfns),imgfn))
    hdu=0
    #write mask-2 and wt map to junk/
    maskfn= os.path.join(os.path.dirname(imgfn).replace('pimage','mask'),\
                            os.path.basename(imgfn).replace('_scie_','_mask_'))
    invvar= read_invvar(imgfn,maskfn,hdu) #note, all post processing on image,mask done in read_invvar
    mask= read_dq(maskfn,hdu)
    maskfn= os.path.join('junk',os.path.basename(maskfn))
    invvarfn= os.path.join('junk',os.path.basename(imgfn).replace('_scie_','_invvar_'))
    fitsio.write(maskfn, mask)
    fitsio.write(invvarfn, invvar)
    print('wrote mask-2 to %s, invvar to %s' % (maskfn,invvarfn))
    #run se
    magzp  = ptf_zeropoint(imgfn)
    hdr=fitsio.read_header(imgfn,ext=hdu)
    seeing = hdr['PIXSCALE'] * hdr['MEDFWHM']
    gain= hdr['GAIN']
    print('magzp= ',magzp,'seeing= ',seeing,'gain= ',gain)
    if args.se_outdir:
        sefn= os.path.join(args.se_outdir,os.path.basename(imgfn).replace('.fits','.se_cat'))
        #
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
        print('###RUNNING SExtractor with: ',cmd)
        if os.system(cmd):
            raise RuntimeError('Command failed: ' + cmd)
        print('SExtractor finished, now PSFex')
    #do PSFex
    if args.psfx_outdir:
        cmd= ' '.join(['psfex',sefn,'-c', os.path.join(args.configdir,'DECaLS.psfex'),
                        '-PSF_DIR',args.psfx_outdir])
        print('###RUNNING PSFex with: ',cmd)
        if os.system(cmd):
            raise RuntimeError('Command failed: ' + cmd)                   
        print('PSFex finished, now Tractor')
    #do legacypipe/tractor
    if args.tractor_outdir:
        cmd= ' '.join(['python','legacypipe/runptf.py','--radec %.2f %.2f' % (args.radec[0],args.radec[1]),\
                        '--decals-dir %s' % args.decals_dir,'--outdir %s' % args.tractor_outdir,\
                        '--pixpsf --splinesky --no-sdss --no-wise --force-all', '--pixscale %.2f' % hdr['PIXSCALE'],\
                        '--width %d --height %d' % (args.wh[0],args.wh[1])])
        if args.coadd_bw: cmd+= ' --coadd-bw'
        print('###RUNNING TRACTOR with: ',cmd)
        if os.system(cmd):
            raise RuntimeError('Command failed: ' + cmd)
print('tractor finished, ...removing junk/ directory')
for fn in glob.glob('junk/*'): os.remove(fn)
os.rmdir('junk')
print('done')
