import os
import fitsio
import argparse
import glob

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-indir",action="store",help='inpute directory, path to mask files')
parser.add_argument("-search",action="store",help='search string including wildcard for mask fits files')
parser.add_argument("-outdir",action="store",help='output dir, where new mask files will be written')
args = parser.parse_args()

# grab header values...
maskfns= glob.glob(os.path.join(args.indir,args.search))
print('found %d mask files' % len(maskfns))
if len(maskfns) == 0: raise ValueError
for cnt,maskfn in enumerate(maskfns):
    print('reading %d of %d masks: %s' % (cnt,len(maskfns),maskfn))
    #read and save bit mask with 2 subtracted off
    mask,mhd= fitsio.read(maskfn,ext=0,header=True)
    mask[mask > 0]= mask[mask > 0]-2
    newfn= os.path.join(args.outdir,os.path.basename(maskfn))
    fitsio.write(newfn,mask,header=mhd)
    print('wrote %d of %d new masks: %s' % (cnt,len(maskfns),newfn))
 
