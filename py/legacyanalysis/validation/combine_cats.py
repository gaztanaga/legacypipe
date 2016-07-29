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
    fns= read_lines(cat_list) 
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
        # Store tractor catalogue and extra data from it as two astropy Tables
        self.tractor, self.meta= combine_cats(cat_list, debug=debug)
        self.clean_up()
        self.extra= self.get_extra_data()
        # Dict of boolean arrays to filter on
        self.b_arrays= self.get_keep_dict()
        # Dict to hold cut astropy tables
        self.data= dict(tractor=None,extra=None)
        # Initialize without cuts, Creates self.data dict which has 'tractor','extra' cut astropy tables
        self.apply_cut(['all'])
        # Outdir for plots
        self.outdir= get_outdir(comparison)

    def clean_up(self):
        # "PSF" not 'PSF '
        self.tractor['type']= np.char.strip(self.tractor['type']) 
        
    def get_extra_data(self):
        '''computes additional data using tractor catalogue, returns as astropy table'''
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
        '''current -- the currently set boolean array'''
        # True where want to keep
        return dict(all= np.ones(len(self.tractor)).astype(bool),\
                    clean= self.clean_cut(),\
                    psf= self.tractor['type'].data == 'PSF',\
                    simp= self.tractor['type'].data == 'SIMP',\
                    exp= self.tractor['type'].data == 'EXP',\
                    dev= self.tractor['type'].data == 'DEV',\
                    comp= self.tractor['type'].data == 'COMP',\
                    extended= self.tractor['type'].data != 'PSF')
        #self.targets= self.get_TargetSelectin(data)

    def clean_cut(self): 
        return  np.any((self.tractor['decam_flux'].data[:,1] > 0,\
                        self.tractor['decam_flux'].data[:,2] > 0,\
                        self.tractor['decam_flux'].data[:,4] > 0, \
                        self.tractor['decam_anymask'].data[:,1] == 0,\
                        self.tractor['decam_anymask'].data[:,2] == 0,\
                        self.tractor['decam_anymask'].data[:,4] == 0,\
                        self.tractor['decam_fracflux'].data[:,1] <= 0.05,\
                        self.tractor['decam_fracflux'].data[:,2] <= 0.05,\
                        self.tractor['decam_fracflux'].data[:,4] <= 0.05,\
                        self.tractor['brick_primary'].data == True),axis=0)

    def get_cut(self,list_of_names):
        '''list_of_names -- keys from b_arrays dict use as a subset
        return boolean array to be used for cut'''
        # Check that keep_list is a list
        assert("<type 'list'>" == repr(type(list_of_names)))
        # Union of cuts
        keep= self.b_arrays[list_of_names[0]]
        if len(list_of_names) > 1:
            for i,name in enumerate(list_of_names[1:]): 
                keep= np.all((keep, self.b_arrays[name]), axis=0)
        return keep
    
    def apply_cut(self,list_of_names):
        '''cut astropy tables to keep_list bool arrays, store in 
        self.data'''
        keep= self.get_cut(list_of_names)
        self.data['tractor']= self.tractor[keep] 
        self.data['extra']= self.extra[keep] 

    def add_myown_cut(self, name=None,b_array=None):
        '''add "name" to self.b_arrays dictionary'''
        assert(name is not None and b_array is not None)
        if name in self.b_arrays.keys(): 
            print("choose a different name for %s, already in self.b_arrays.keys()" % name)
            raise ValueError
        self.b_arrays['name'] = b_array
    
    def n_in_cut(self):
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
        # Add bool arrays for matched, unmatched sources
        m_dict= self.do_matching()
        self.ref.add_myown_cut(name='match',b_array=m_dict['ref_match'])
        self.ref.add_myown_cut(name='unmatch',b_array=m_dict['ref_miss'])
        self.test.add_myown_cut(name='match',b_array=m_dict['test_match'])
        self.test.add_myown_cut(name='unmatch',b_array=m_dict['test_miss'])
        # Use all sources by default
        self.apply_cut(['all'])   
        # Get meta data, like output dir
        self.outdir= self.ref.outdir 

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
 
    def apply_cut(self,keep_list):
        '''call the apply_cut atrribute of self.ref and self.test'''
        if 'match' in keep_list:
            # Ref gflux > 0, test gflux > 0 could through out diff number objects in Ref,Test prevent this
            keep= np.all((self.ref.get_cut(keep_list),\
                          self.test.get_cut(keep_list)), axis=0)
            self.ref.add_myown_cut(name='both',b_array=keep)
            self.test.add_myown_cut(name='both',b_array=keep)
            keep_list= ['both']
        # Apply cut to for unique b_arrays of ref and test
        self.ref.apply_cut(keep_list)
        self.test.apply_cut(keep_list)


