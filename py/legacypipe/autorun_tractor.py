import os
import argparse
import glob

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-camera",choices=['ptf','decam'],action="store",help='sets pixscale')
parser.add_argument("-radec",nargs=2,type=float,action="store",help='for legacypipe')
parser.add_argument("-decals_dir",action="store",help='for legacypipe')
parser.add_argument("-outdir",action="store",help='set to run legacypipe/tractor', required=False)
parser.add_argument("-wh",nargs=2,type=int,action="store",help='for legacypipe')
parser.add_argument("-force_all",action="store",help='for legacypipe, any value calls it',required=False)
parser.add_argument("-coadd_bw",action="store",help='for legacypipe, any value calls it',required=False)
args = parser.parse_args()

def make_dir(name):
    if not os.path.exists(name): os.makedirs(name)

make_dir(args.outdir)
if args.camera == 'ptf': pixscale=1.01
else: raise ValueError
cmd= ' '.join(['python','legacypipe/runptf.py','--radec %.2f %.2f' % (args.radec[0],args.radec[1]),\
                '--decals-dir %s' % args.decals_dir,'--outdir %s' % args.outdir,\
                '--pixpsf --splinesky --no-sdss --no-wise',\
                '--width %d --height %d' % (args.wh[0],args.wh[1])])
if args.force_all: cmd+= ' --force-all'
if args.coadd_bw: cmd+= ' --coadd-bw'
print('###RUNNING TRACTOR with: ',cmd)
if os.system(cmd):
    raise RuntimeError('Command failed: ' + cmd)
print('done')
