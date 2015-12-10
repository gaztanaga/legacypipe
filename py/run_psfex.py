import os
import numpy as np
import fitsio
import argparse
import glob

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-secat_search",action="store",help='path plus wildcard for sextractor catalogs')
parser.add_argument("-config_dir",action="store",help='config pile')
parser.add_argument("-outdir",action="store",help='psfex output file name')
args = parser.parse_args()

sefns= glob.glob(args.secat_search)
print('found %d files' % len(sefns))
if len(sefns) == 0: raise ValueError
for cnt,sefn in enumerate(sefns):
    print('reading %d of %d sextractor catalogues: %s' % (cnt,len(sefns),sefn))
    cmd= ' '.join(['psfex',sefn,'-c', os.path.join(args.config_dir,'DECaLS.psfex'),
                    '-PSF_DIR',args.outdir])
    print(cmd)
    if os.system(cmd):
        raise RuntimeError('Command failed: ' + cmd)                   

