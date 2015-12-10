import os
import numpy as np
import fitsio
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-se_cat",action="store",help='sextractor catalogue',required=False)
parser.add_argument("-psf_cat",action="store",help='psfex catalogue',required=False)
args = parser.parse_args()

if args.psf_cat or args.se_cat:
    pass #one of these was specified
else: 
    print "WARNING: must specify sextractor or psfex catalog or both"
    raise ValueError

if args.se_cat:
    cat=fitsio.FITS(args.se_cat)
    #h2= cat[2].read_header()  #numpy array shape = 1, it is a mess of words
    results= cat[2].read()
    print '%d sextractor objects, %d properties for each object' % (len(results),len(cat[2].get_colnames()))
    print 'the properties are: '
    for p in cat[2].get_colnames(): print p

    plt.subplots_adjust(hspace=0, wspace=0,
                        left=0.05, right=0.95, bottom=0.05, top=0.95)
    plt.clf()
    rows,cols = 5,5
    for i in range(rows*cols):
        plt.subplot(rows, cols, i+1)
        plt.imshow(results['VIGNET'][i+1000],cmap=cm.gray)
    plt.savefig('vignets.png')
    plt.clf()
if args.psf_cat:
    cat=fitsio.FITS(args.psf_cat)
    #h2= cat[2].read_header()  #numpy array shape = 1, it is a mess of words
    results= cat[1].read()
    results= results[0][0]
    print '%d psfex objects which is <= number of sextractor objects' % results.shape[0]

    #plt.subplots_adjust(hspace=0, wspace=0,
    #                    left=0.05, right=0.95, bottom=0.05, top=0.95)
    #plt.clf()
    #rows,cols = 5,5
    #for i in range(rows*cols):
    #    plt.subplot(rows, cols, i+1)
    #    plt.imshow(results['VIGNET'][i],cmap=cm.gray)
    #plt.savefig('vignets_%d.png' % rows*cols)
    #plt.clf()
   
