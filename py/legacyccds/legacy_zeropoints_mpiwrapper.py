from __future__ import division, print_function

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import os
import numpy as np
from glob import glob
import datetime
import sys

from astrometry.util.ttime import Time
from legacyccds.legacy_zeropoints import get_parser,ptime,read_lines
from legacyccds.legacy_zeropoints import main as legacy_main

######## 
# stdouterr_redirected() is from Ted Kisner
# Every mpi task (zeropoint file) gets its own stdout file
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
            if not os.path.exists(os.path.dirname(pto)):
                os.makedirs(os.path.dirname(pto))
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
            

if __name__ == "__main__":
    t0 = Time()
    tbegin = t0
    parser = get_parser() 
    args = parser.parse_args()

    if args.nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
    t0=ptime('parse-args',t0)

    images= read_lines(args.image_list) 
    if args.nproc > 1:
        image_list= np.array_split(images, comm.size)[comm.rank]
    else:
        image_list= np.array_split(images, 1)[0]
    if args.nproc > 1:
        # Log to unique file
        outfn=os.path.join(args.outdir,"logs",\
                           "%s" % datetime.datetime.now().strftime("%Y-%m-%d"),\
                           "log.%d_%s" % datetime.datetime.now().strftime("%H-min%M")) 
        with stdouterr_redirected(to=outfn, comm=None):  
            t0=ptime('before-%s' % brick,t0)
            legacy_main(image_list=image_list, args=args) 
            t0=ptime('after-%s' % brick,t0)
    else:
        legacy_main(image_list=image_list, args=args) 
    if args.nproc > 1:
        # Wait for all mpi tasks to finish 
        confirm_files = comm.gather( image_list, root=0 )
        if comm.rank == 0:
            tnow= Time()
            print("Done, total time = %s" % (tnow-tbegin,))
    else:
        tnow= Time()
        print("Done, total time = %s" % (tnow-tbegin,))

