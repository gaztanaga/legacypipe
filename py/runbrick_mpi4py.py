from __future__ import print_function
from mpi4py.MPI import COMM_WORLD as comm
import argparse
import os
import numpy as np
from glob import glob
import subprocess
from legacypipe.runbrick import main

######## 
## Ted's
import sys
import time
from contextlib import contextmanager

@contextmanager
def stdouterr_redirected(to=os.devnull, comm=None):
    '''
    Based on http://stackoverflow.com/questions/5081657
    import os
    with stdouterr_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    sys.stdout.flush()
    sys.stderr.flush()
    fd = sys.stdout.fileno()
    fde = sys.stderr.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd
        sys.stderr.close() # + implicit flush()
        os.dup2(to.fileno(), fde) # fd writes to 'to' file
        sys.stderr = os.fdopen(fde, 'w') # Python writes to fd
        
    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        if (comm is None) or (comm.rank == 0):
            print("Begin log redirection to {} at {}".format(to, time.asctime()))
        sys.stdout.flush()
        sys.stderr.flush()
        pto = to
        if comm is None:
           # if not os.path.exists(os.path.dirname(pto)):
           #     os.makedirs(os.path.dirname(pto))
            with open(pto, 'w') as file:
                _redirect_stdout(to=file)
        else:
            pto = "{}_{}".format(to, comm.rank)
            with open(pto, 'w') as file:
                _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different
            if comm is not None:
                # concatenate per-process files
                comm.barrier()
                if comm.rank == 0:
                    with open(to, 'w') as outfile:
                        for p in range(comm.size):
                            outfile.write("================= Process {} =================\n".format(p))
                            fname = "{}_{}".format(to, p)
                            with open(fname) as infile:
                                outfile.write(infile.read())
                            os.remove(fname)
                comm.barrier()

            if (comm is None) or (comm.rank == 0):
                print("End log redirection to {} at {}".format(to, time.asctime()))
            sys.stdout.flush()
            sys.stderr.flush()
            
    return
##############

def file_contains(fn,text):
    '''return True if fn contains the string text, False otherwise'''
    hasit=subprocess.Popen(['grep',text,fn], stdout= subprocess.PIPE)
    hasit=hasit.stdout.read()
    if text in hasit: 
        return True
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="test")
    parser.add_argument("--brick_list",action="store",required=True,help='text file listing all bricks to run')
    parser.add_argument("--outdir",action="store",required=True)
    parser.add_argument("--jobid", action="store",required=True)
    parser.add_argument("--threads", type=int, action="store",default=1,required=False)
    parser.add_argument("--zoom", type=int, nargs=4,action="store",required=False)
    opt = parser.parse_args()

    # Each process does this:
    brick_list=np.loadtxt(opt.brick_list,dtype=str)
    split_list=np.array_split(brick_list,comm.size)
    brick_list= split_list[comm.rank]
    print('rank=%d, working on bricks:' % comm.rank,brick_list)
    for brick in brick_list:
        bri=brick[:3]
        # Check if Tractor Cat exists?
        tractor_fn=os.path.join(opt.outdir,'tractor/%s/tractor-%s.fits' % (bri,brick))
        if os.path.exists(tractor_fn):
            print('Skipping Brick %s: tractor cat exists!' % brick)
            continue
        # Check if no ccds touch it
        outfn_pat=os.path.join(opt.outdir,'logs/','b%s_jobid*.o' % (brick,))
        fils=glob(outfn_pat)
        if len(fils) > 0: 
            fils=fils[0]
            if np.any((file_contains(fils,'No CCDs touching brick'),\
                       file_contains(fils,'No photometric CCDs touching brick')),axis=0):
                print('Skipping Brick %s: No CCDs touching and/or photometric' %  brick)
                continue
        # Run Legacypipe 
        check_fn=os.path.join(opt.outdir,'checkpoints/%s/%s.pickle' % (bri,brick))
        pick_fn=os.path.join(opt.outdir,"pickles/%s/" % bri,"runbrick-%(brick)s-%%(stage)s.pickle") 
        args_list= ['--brick', brick,\
                   '--skip', \
                   '--threads', str(opt.threads),\
                   '--checkpoint',check_fn,\
                   '--pickle',pick_fn,\
                   '--outdir', opt.outdir, '--nsigma', '6',\
                   '--no-wise']
        if opt.zoom:
           args_list+= ['--zoom', '%d' % opt.zoom[0], '%d' % opt.zoom[1], '%d' % opt.zoom[2], '%d' % opt.zoom[3]]
        outfn=os.path.join(opt.outdir,'logs/','b%s_jobid%s.o' % (brick,opt.jobid))
        if not os.path.exists(os.path.dirname(outfn)):
            try:
                os.makedirs(os.path.dirname(outfn))
            except OSError:
                print('good, %s already exists' % os.path.dirname(outfn))
        with stdouterr_redirected(to=outfn, comm=None):
            main(args=args_list)
 



