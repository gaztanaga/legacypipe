from __future__ import print_function
from mpi4py.MPI import COMM_WORLD as comm
import argparse
import os
import numpy as np
from legacypipe.runbrick import main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="test")
    parser.add_argument("--brick_list",action="store",required=True,help='text file listing all bricks to run')
    parser.add_argument("--outdir",action="store",required=True)
    parser.add_argument("--threads", type=int, action="store",default=1,required=False)
    parser.add_argument("--zoom", type=int, nargs=4,action="store",required=False)
    opt = parser.parse_args()

    brick_list=np.loadtxt(opt.brick_list,dtype=str)
    split_list=np.split(brick_list,comm.size)
    brick_list= split_list[comm.rank]
    print('rank=%d, working on bricks:' % comm.rank,brick_list)
    for brick in brick_list:
        bri=brick[:3]
        check_fn=os.path.join(opt.outdir,'checkpoints/%s/checkpoint-%s.pickle' % (bri,brick))
        pick_fn=os.path.join(opt.outdir,"pickles/%s/" % bri,"runbrick-%(brick)s-%%(stage)s.pickle") 
        args_list= ['--brick', brick,\
                   '--skip', '--threads', str(opt.threads),\
                   '--checkpoint',check_fn,\
                   '--no-wise',\
                   '--force-all',\
                   '--pickle',pick_fn,\
                   '--outdir', opt.outdir]
        if opt.zoom:
           args_list+= ['--zoom', '%d' % opt.zoom[0], '%d' % opt.zoom[1], '%d' % opt.zoom[2], '%d' % opt.zoom[3]]
        main(args=args_list)
 



