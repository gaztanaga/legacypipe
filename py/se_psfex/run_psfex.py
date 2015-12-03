import os
import numpy as np
import fitsio
import argparse
import glob

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-secat_search",action="store",help='sextractor output catalog')
parser.add_argument("-psfex_config",action="store",help='config pile')
parser.add_argument("-outdir",action="store",help='psfex output file name')
args = parser.parse_args()

sefns= glob.glob(os.path.join('./se_output',args.secat_search))
print('found %d files' % len(sefns))
if len(sefns) == 0: raise ValueError
for cnt,sefn in enumerate(sefns):
    print('reading %d of %d sextractor catalogues: %s' % (cnt,len(sefns),sefn))
    cmd= ' '.join(['psfex',sefn,'-c', args.psfex_config,
                    '-PSF_DIR',args.outdir])
                    #'-SAMPLE_FLAGMASK 254',
                    #'-SAMPLE_IMAFLAGMASK 254'])
    print(cmd)
    if os.system(cmd):
        raise RuntimeError('Command failed: ' + cmd)                   

