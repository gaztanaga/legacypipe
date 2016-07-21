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
        self.more= self.get_more()
        self.keep= self.get_keep_dict()
        # Initialize without cuts
        # Creates self.t,self.m tables
        self.select(['all'])
        # Outdir for plots
        self.outdir= get_outdir(comparison)

    def clean_up(self):
        # "PSF" not 'PSF '
        self.tractor['type']= np.char.strip(self.tractor['type']) 
        
    def get_more(self):
        # AB Mags
        mag= self.get_magAB()
        mag_ivar= self.get_magAB_ivar()
        return Table([mag, mag_ivar], names=('mag', 'mag_ivar'))
    
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
        # Matching indices, as dict of boolean arrays
        self.m_dict= self.do_matching()

    def do_matching(self):
        '''matches test ra,dec to ref ra,dec
        Returns dict of boolean arrays for the matched and missed samples'''    
        m1, m2, d12 = match_radec(self.ref.tractor['ra'].data.copy(), self.ref.tractor['dec'].data.copy(),\
                                  self.test.tractor['ra'].data.copy(), self.test.tractor['dec'].data.copy(),\
                                  1.0/3600.0)
        print("Matched: %d/%d objects" % (m1.size,len(self.ref.tractor)))
        miss1 = np.delete(np.arange(len(self.ref.tractor)), m1, axis=0)
        miss2 = np.delete(np.arange(len(self.test.tractor)), m2, axis=0)
        # Indices to bool array
        ref_rows= len(self.ref.tractor)
        test_rows= len(self.test.tractor)
        return dict(ref_match= self.indices2bool(m1,rows=ref_rows),\
                    ref_miss= self.indices2bool(miss1,rows=ref_rows),\
                    test_match= self.indices2bool(m2,rows=test_rows),\
                    test_miss= self.indices2bool(miss2,rows=test_rows))

    def indices2bool(self,indices,rows=None):
        '''convert indices for an array into boolean array of same length as array'''
        assert(rows is not None)
        b_arr= np.zeros(rows).astype(bool)
        b_arr[indices]= True
        return b_arr
 
    def use_matched(self):
        self.ref.keep['permanent']= self.keep['ref_match']
        self.test.keep['permanent']= self.keep['test_match']
    
    def use_missed(self):
        self.ref.keep['permanent']= self.keep['ref_miss']
        self.test.keep['permanent']= self.keep['test_miss']
          
    def select(self,keep_list):
        '''select keys in keep_list from ref and test data sets'''
        self.ref.select(keep_list)
        self.test.select(keep_list)


