
import numpy as np
import fitsio
import glob
import os
from astropy import units
from astropy.coordinates import SkyCoord

from astrometry.util.fits import fits_table, merge_tables
from tractor.brightness import NanoMaggies

def read_lines(fn_list):
    fin=open(fn_list,'r')
    lines=fin.readlines()
    fin.close()
    return np.sort(np.array( list(np.char.strip(lines)) ))

class CatalogueFuncs(object):
    '''funcs for using Dustins fits_table objects'''
    def stack(self,fn_list):
        '''concatenates fits tables'''
        cats= []
        fns=read_lines(fn_list)
        if len(fns) < 1: raise ValueError('Error: fns=',fns)
        for fn in fns:
            cats.append( fits_table(fn) )
        return merge_tables(cats, columns='fillzero')
    
    def set_extra_data(self,cat):
        '''adds columns to fits_table cat'''
        # Remove white spaces
        cat.set('type', np.char.strip(cat.get('type')))
        # AB Mags
        shp=cat.get('decam_flux').shape
        mag,mag_sigma= np.zeros(shp),np.zeros(shp)
        for iband in range(shp[1]):
            mag[:,iband],mag_sigma[:,iband]=NanoMaggies.fluxErrorsToMagErrors(\
                        cat.get('decam_flux')[:,iband], cat.get('decam_flux_ivar')[:,iband])
        cat.set('decam_mag',mag)
        cat.set('decam_mag_ivar',1./np.power(mag_sigma,2))
        
class Cuts4MatchedCats(object):
    '''Cuts for MATCHED cats only'''
    def __init__(self,matched1,matched2):
        self.psf1 = (matched1.get('type') == 'PSF')
        self.psf2 = (matched2.get('type') == 'PSF')
        # Band dependent
        bands='ugrizY'
        self.good={}
        for band,iband in zip(bands,range(6)):
            self.good[band]= ((matched1.decam_flux_ivar[:,iband] > 0) *\
                              (matched2.decam_flux_ivar[:,iband] > 0))      
        

#     def get_mags(self,cat):
#         mag= cat.get('decam_flux')/cat.get('decam_mw_transmission')
#         return 22.5 -2.5*np.log10(mag)

#     def get_mags_ivar(self,cat):
#         return np.power(np.log(10.)/2.5*cat.get('decam_flux'), 2)* \
#                                     cat.get('decam_flux_ivar')

class Matcher(object):
    '''sphere matches two ra,dec lists,
    ref,obs are astronometry.net "fits_tables" objects
    '''
    def __init__(self):
        pass
    def match_within(self,ref,obs,dist=1./3600):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        idx, d2d, d3d = cat1.match_to_catalog_3d(cat2)
        b= np.array(d2d) <= dist
        # Return 4 index arrays for indices where matches and where missing
        imatch=dict(ref=[],obs=[])
        imiss=dict(ref=[],obs=[])
        # 
        imatch['ref']= np.arange(len(ref))[b]
        imatch['obs']= np.array(idx)[b]
        print("Matched: %d/%d objects" % (imatch['ref'].size,len(ref)))
        imiss['ref'] = np.delete(np.arange(len(ref)), imatch['ref'], axis=0)
        imiss['obs'] = np.delete(np.arange(len(obs)), imatch['obs'], axis=0)
        return imatch,imiss,d2d

    def nearest_neighbors_within(self,ref,obs,within=1./3600,min_nn=1,max_nn=5):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        ref_nn,obs_nn,dist={},{},{}
        for nn in range(min_nn,max_nn+1):
            idx, d2d, d3d = cat1.match_to_catalog_3d(cat2,nthneighbor=nn)
            b= np.array(d2d) <= within
            ref_nn[str(nn)]= np.arange(len(ref))[b]
            obs_nn[str(nn)]= np.array(idx)[b]
            dist[str(nn)]= np.array(d2d)[b]
            print("within 1arcsec, nn=%d, %d/%d" % (nn,ref_nn[str(nn)].size,len(ref)))
        return ref_nn,obs_nn,dist
 
