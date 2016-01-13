import os
import fitsio
import argparse
import glob
import numpy as np

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-mask_search",action="store",help='relative path + search string for mask fits files')
parser.add_argument("-mask_outdir",action="store",help='relative path to output dir for mask - 2 images')
parser.add_argument("-img_search",action="store",help='relative path + search string for image files')
parser.add_argument("-invvar_outdir",action="store",help='path to output dir for invvar images')
args = parser.parse_args()

# write mask-2
maskfns= glob.glob(args.mask_search)
print('found %d mask files' % len(maskfns))
if len(maskfns) == 0: raise ValueError
for cnt,maskfn in enumerate(maskfns):
    print('reading %d of %d masks: %s' % (cnt,len(maskfns),maskfn))
    #read and save bit mask with 2 subtracted off
    mask,mhd= fitsio.read(maskfn,ext=0,header=True)
    mask[mask > 0]= mask[mask > 0]-2
    newfn= os.path.join(args.mask_outdir,os.path.basename(maskfn))
    fitsio.write(newfn,mask,header=mhd)
    print('wrote %d of %d new masks: %s' % (cnt,len(maskfns),newfn))
# write invvar
imgfns= glob.glob(args.img_search)
print('found %d image files' % len(imgfns))
if len(imgfns) == 0: raise ValueError
for cnt,imgfn in enumerate(imgfns):
    print('reading %d of %d images: %s' % (cnt,len(imgfns),imgfn))
    #get image and header
    img,hdr= fitsio.read(imgfn,ext=0,header=True)
    #get mask-2
    mask2fn = os.path.basename(imgfn).replace('_scie_', '_mask_')
    mask2= fitsio.read(os.path.join(args.mask_outdir,mask2fn),ext=0,header=False)
    #compute invvar and save to file
    invvar=np.zeros(img.shape)
    invvar[mask2 == 0]= hdr['GAIN']/img[mask2 == 0]
    name= os.path.basename(imgfn).replace('_scie_', '_invvar_') 
    fitsio.write(os.path.join(args.invvar_outdir,name),invvar,header=hdr)
    print('wrote %d of %d invvar: %s' % (cnt,len(imgfns),name))
 
