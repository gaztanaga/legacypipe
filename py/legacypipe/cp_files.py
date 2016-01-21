import os
import argparse
import glob
import fitsio
import numpy as np

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-Ncopy",type=int,action="store",help='# of images to copy from /project to /scratch')
parser.add_argument("-band",choices=['g','R'],action="store",help='which band')
parser.add_argument("-pimage_dir",action="store",help='/project processed image dir')
parser.add_argument("-mask_dir",action="store",help='same but for bit mask images')
parser.add_argument("-outdir",action="store",help='directory to make pimage/ and mask/ in for storing the copied images')
args = parser.parse_args()

for mydir in ['pimage','mask']:
    if not os.path.exists(os.path.join(args.outdir,mydir)): os.makedirs(os.path.join(args.outdir,mydir))

search_str= '*_f0%d_*.fits' % 2
if args.band == 'g': search_str= '*_f0%d_*.fits' % 1
copied=0
for fn in glob.glob(os.path.join(args.pimage_dir,search_str)):
    hdr=fitsio.read_header(fn,ext=0)
    if type(hdr['IMAGEZPT']) == np.float:    
        copied +=1
        cmd= ' '.join(['cp',fn,os.path.join(args.outdir,'pimage/',os.path.basename(fn))])
        if os.system(cmd): raise ValueError
        maskfn= os.path.join(args.mask_dir, os.path.basename(fn).replace('_scie_','_mask_'))
        cmd= ' '.join(['cp',maskfn,os.path.join(args.outdir,'mask/',os.path.basename(maskfn))])
        if os.system(cmd): raise ValueError
    if copied >= args.Ncopy: break
print "done"
