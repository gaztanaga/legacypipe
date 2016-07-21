#!/usr/bin/env python

"""combine list of tractor cats and match them along the way

Example of Matched Comparison:
from legacyanalysis.combine_cats import Matched_DataSet
import common_plots as plots
d= Matched_DataSet(list_of_ref_cats,list_of_test_cats)
plots.confusion_matrix(d.ref_matched.t,d.test_matched.t)
"""

from __future__ import division, print_function

import matplotlib
matplotlib.use('Agg') #display backend
import os
import sys
import logging
import numpy as np

from astropy.io import fits
from astropy.table import vstack, Table, Column
#from scipy.spatial import KDTree
from astrometry.libkd.spherematch import match_radec

from legacyanalysis.validation.pathnames import get_outdir

def read_lines(fn):
    fin=open(fn,'r')
    lines=fin.readlines()
    fin.close()
    return list(np.char.strip(lines))

def deg2_lower_limit(ra,dec):
    '''deg2 spanned by objects in each data set, lower limit'''
    ra_wid= ra.max()-ra.min()
    assert(ra_wid > 0.)
    dec_wid= abs(dec.max()-dec.min())
    return ra_wid*dec_wid

def kdtree_match(ref_ra,ref_dec, ra,dec, k=1, dsmax=1./3600):
    '''finds approx. distance between reference and test points and
    returns array of indices and distance where distance < dsmax (degrees)
    '''
    assert(len(ref_ra) == len(ref_dec))
    assert(len(ra) == len(dec))
    # Index the tree
    tree = KDTree(np.transpose([dec.copy(),ra.copy()])) 
    # Find the k nearest neighbors to each reference ra,dec
    ds, i_tree = tree.query(np.transpose([ref_dec.copy(),ref_ra.copy()]), k=k) 
    # Get indices where match
    i_ref,i_other={},{}
    ref_match= np.arange(len(ref_ra))[ds<=dsmax]
    other_match= i_tree[ds<=dsmax]
    assert(len(ref_ra[ref_match]) == len(ra[other_match]))
    return np.array(ref_match),np.array(other_match),np.array(ds[ds<=dsmax])

def combine_cats(cat_list, debug=False):
    '''return dict containing astropy Table of concatenated tractor cats'''
    #get lists of tractor cats to compare
    fns= read_lines(ref_cats_file) 
    print('Combining tractor catalogues: ')
    for fn in fns: print("%s" % fn) 
    #object to store concatenated matched tractor cats
    bigtractor = []
    deg2= 0.
    # One catalogue for quick debugging
    if debug: fns= fns[:1]
    # Loop over cats
    for cnt,fn in zip(range(len(fns)),fns):
        print('Reading %s' % fn)
        tractor = Table(fits.getdata(fn, 1), masked=True)
        # Build combined catalogs
        if len(bigtractor) == 0:
            bigtractor = tractor
        else:
            bigtractor = vstack((bigtractor, tractor), metadata_conflicts='error')
        deg2+= deg2_lower_limit(tractor['ra'],tractor['dec'])
    
    return bigtractor, dict(deg2=deg2)


def match_two_cats(ref_cats_file,test_cats_file, debug=False):
    '''return dict containing astropy Table of concatenated tractor cats
    one Table for matched and missed reference and test objects, each'''
    # Set the debugging level
    lvl = logging.INFO
    logging.basicConfig(format='%(message)s', level=lvl, stream=sys.stdout)
    log = logging.getLogger('__name__')

    #get lists of tractor cats to compare
    fns_1= read_lines(ref_cats_file) 
    fns_2= read_lines(test_cats_file) 
    print('Comparing tractor catalogues: ')
    for one,two in zip(fns_1,fns_2): print("%s -- %s" % (one,two)) 
    #if fns_1.size == 1: fns_1,fns_2= [fns_1],[fns_2]
    #object to store concatenated matched tractor cats
    ref_matched = []
    ref_missed = []
    test_matched = []
    test_missed = []
    d_matched= 0.
    deg2= dict(ref=0.,test=0.,matched=0.)
    # One catalogue for quick debugging
    if debug: fns_1,fns_2= fns_1[:1],fns_2[:1]
    # Loop over cats
    for cnt,cat1,cat2 in zip(range(len(fns_1)),fns_1,fns_2):
        print('Reading %s -- %s' % (cat1,cat2))
        ref_tractor = Table(fits.getdata(cat1, 1), masked=True)
        test_tractor = Table(fits.getdata(cat2, 1), masked=True)
        m1, m2, d12 = match_radec(ref_tractor['ra'].data.copy(), ref_tractor['dec'].data.copy(),\
                                  test_tractor['ra'].data.copy(), test_tractor['dec'].data.copy(), \
                                  1.0/3600.0)
        #m1, m2, d12= kdtree_match(ref_tractor['ra'].copy(), ref_tractor['dec'].copy(),\
        #                          test_tractor['ra'].copy(), test_tractor['dec'].copy(),\
        #                          k=1, dsmax=1./3600)
        print("Matched: %d/%d objects" % (m1.size,len(ref_tractor['ra'])))
        miss1 = np.delete(np.arange(len(ref_tractor)), m1, axis=0)
        miss2 = np.delete(np.arange(len(test_tractor)), m2, axis=0)

        # Build combined catalogs
        if len(ref_matched) == 0:
            ref_matched = ref_tractor[m1]
            ref_missed = ref_tractor[miss1]
            test_matched = test_tractor[m2]
            test_missed = test_tractor[miss2]
            d_matched= d12
        else:
            ref_matched = vstack((ref_matched, ref_tractor[m1]), metadata_conflicts='error')
            ref_missed = vstack((ref_missed, ref_tractor[miss1]), metadata_conflicts='error')
            test_matched = vstack((test_matched, test_tractor[m2]), metadata_conflicts='error')
            test_missed = vstack((test_missed, test_tractor[miss2]), metadata_conflicts='error')
            d_matched= np.concatenate([d_matched, d12])
        deg2['ref']+= deg2_lower_limit(ref_tractor['ra'],ref_tractor['dec'])
        deg2['test']+= deg2_lower_limit(test_tractor['ra'],test_tractor['dec'])
        deg2['matched']+= deg2_lower_limit(ref_matched['ra'],ref_matched['dec'])
    
    return dict(ref_matched = ref_matched,
                ref_missed = ref_missed,
                test_matched = test_matched,
                test_missed = test_missed), \
           dict(d_matched= d_matched, deg2= deg2)

class Single_DataSet(object):
    '''Has five things:
    1) self.tractor -- astropy table
    2) self.more -- astropy table with additional info e.g. decam_mag, decam_mag_ivar etc
    3) Select -- dict of indices for subsets of data, e.g. PSF
    4) self.t, self.m -- subset astropy tables of self.tractor, self.more using Select()
    and supporting functions
    '''
    def __init__(self,cat_list, comparison='test',debug=False):
        # Store tractor catalogue as astropy Table
        self.tractor, self.meta= combine_cats(cat_list, debug=debug)
        self.clean_up()
        self.more= self.get_more_table()
        self.keep= self.get_keep_dict()
        # Initialize without cuts
        # Creates self.t,self.m tables
        self.select(['all'])
        # Outdir for plots
        self.outdir= get_outdir(comparison)

    def clean_up(self):
        # "PSF" not 'PSF '
        self.tractor['type']= np.char.strip(self.tractor['type']) 
        
    def get_more_table(self):
        # AB Mags
        mag= self.get_magAB()
        mag_ivar= self.get_magAB_ivar()
    
    def get_magAB(self):
        mag= self.tractor['decam_flux']/self.tractor['decam_mw_transmission']
        return 22.5 -2.5*np.log10(mag)

    def get_magAB_ivar(self):
        return np.power(np.log(10.)/2.5*self.tractor['decam_flux'], 2)* \
                                    self.tractor['decam_flux_ivar']
 
    def get_keep_dict(self):
        '''permanent -- a special key in the "keep" dict which is always used if set, default="all"
        current -- the currently set boolean array, self.select() sets self.keep["current"] to what you want
        all other keys only used when explicitly set with "self.select(key_name)"'''
        # True where want to keep
        return dict(permanent= np.ones(len(self.tractor)).astype(bool),\
                    current= np.ones(len(self.tractor)).astype(bool),\
                    default= self.basic_cuts(),\
                    all= np.ones(len(self.tractor)).astype(bool),\
                    psf= self.tractor['type'].data == 'PSF',\
                    simp= self.tractor['type'].data == 'SIMP',\
                    exp= self.tractor['type'].data == 'EXP',\
                    dev= self.tractor['type'].data == 'DEV',\
                    comp= self.tractor['type'].data == 'COMP',\
                    extended= self.tractor['type'].data != 'PSF')
        #self.targets= self.get_TargetSelectin(data)

    def basic_cuts(self): 
        return  np.any((self.t['decam_flux'].data[:,1] > 0,\
                        self.t['decam_flux'].data[:,2] > 0,\
                        self.t['decam_flux'].data[:,4] > 0, \
                        self.t['decam_anymask'].data[:,1] == 0,\
                        self.t['decam_anymask'].data[:,2] == 0,\
                        self.t['decam_anymask'].data[:,4] == 0,\
                        self.t['decam_fracflux'].data[:,1] <= 0.05,\
                        self.t['decam_fracflux'].data[:,2] <= 0.05,\
                        self.t['decam_fracflux'].data[:,4] <= 0.05,\
                        self.t['brick_primary'].data == True),axis=0)

    def select(self,keep_list):
        '''keep_list -- list of names from "keep" dict to use as a subset
        permanent -- a special key in the "keep" dict which is always used if set, default="all"'''
        # Add permanent to list
        keep_list+= ['permanent'] 
        for name in list_of_names: assert(name in self.keep.keys())
        # Union
        keep= self.keep[keep_list[0]]
        if len(keep_list) > 1:
            for i,name in enumerate(keep_list[1:]): 
                keep= np.all((keep, self.keep[name]), axis=0)
        # Apply
        self.apply_selection(keep)

    def select_my_own(self, my_b_arr):
        '''keep your own boolean array'''
        # Add permanent to list
        keep= np.all((my_b_arr, self.keep['permanent']), axis=0)
        # Apply
        self.apply_selection(keep)
    
    def apply_selection(self,b_arr):
        self.keep['current']= b_arr
        self.t= self.tractor[b_arr] 
        self.o= self.tractor[b_arr] 

    def number_kept(self):
        return len(self.t)
    #def get_TargetSelectin(self, data):
    ##### store as self.t.masks['elg'], self.t.masks['lrg'], etc
    #    d={}
    #e    d['desi_target'], d['bgs_target'], d['mws_target']= \
    #                    cuts.apply_cuts( data )
    #    return d

 


class Matched_DataSet(object):
    '''a Matched_DataSet contains a dict of 4 concatenated, matched, tractor catalogue astropy Tables
    each table has additional columns for mags, target selection, masks, etc.'''
    def __init__(self,ref_cats_file,test_cats_file, comparison='test',debug=False):
        # Combine catalogues and get all info could possibly need
        self.ref= Single_DataSet(ref_cats_file, comparison=comparison,debug=debug)
        self.test= Single_DataSet(test_cats_file, comparison=comparison,debug=debug)
        # Match
        m1, m2, d12 = match_radec(self.ref.tractor['ra'].data.copy(), self.ref.tractor['dec'].data.copy(),\
                                  self.test.tractor['ra'].data.copy(), self.test.tractor['dec'].data.copy(),\
                                  1.0/3600.0)
        print("Matched: %d/%d objects" % (m1.size,len(self.ref.tractor)))
        miss1 = np.delete(np.arange(len(self.ref.tractor)), m1, axis=0)
        miss2 = np.delete(np.arange(len(self.test.tractor)), m2, axis=0)
        # Indices to bool array
        self.keep= dict(ref_match= self.indices2bool(m1,size=len(self.ref.tractor)),\
                        ref_miss= self.indices2bool(miss1,size=len(self.ref.tractor),\
                        test_match= self.indices2bool(m2,size=len(self.test.tractor),\
                        test_miss= self.indices2bool(miss2,size=len(self.test.tractor))
        # Use Matched by default 
        self.select('matched')
        
    def select(name):
        assert(name == 'matched' or name == 'missed')
        if name == 'matched':
            self.ref.keep['permanent']= self.keep['ref_match']
            self.test.keep['permanent']= self.keep['test_match']
        elif name == 'matched':
            self.ref.keep['permanent']= self.keep['ref_miss']
            self.test.keep['permanent']= self.keep['test_miss']

    def indices2bool(self,indices,size):
        '''return boolean array of len size corresponding to indices of an array of len size'''
        b_arr= np.zeros(size).astype(bool)
        b_arr[indices]= True
        return b_arr
