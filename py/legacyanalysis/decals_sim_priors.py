if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.io import fits
#from astropy.table import vstack, Table
from astrometry.util.fits import fits_table, merge_tables
import os
import sys
from scipy.optimize import newton

class _GaussianMixtureModel(object):
    """Read and sample from a pre-defined Gaussian mixture model.
    This assumes 'sklearn.mixture.GMM' has already been used to determine MoGs 
    """
    def __init__(self, weights_, means_, covars_, covtype):
        self.weights_ = weights_
        self.means_ = means_
        self.covars_ = covars_
        self.covtype = covtype
        self.n_components, self.n_dimensions = self.means_.shape
    
    @staticmethod
    def save(model, filename,index=None):
        '''index: optional nddex array for subset of compenents to save'''
        hdus = fits.HDUList()
        hdr = fits.Header()
        hdr['covtype'] = model.covariance_type
        if index is None:
            index=np.arange(len(model.weights_))
        hdus.append(fits.ImageHDU(model.weights_[index], name='weights_', header=hdr))
        hdus.append(fits.ImageHDU(model.means_[index,...], name='means_'))
        hdus.append(fits.ImageHDU(model.covars_[index,...], name='covars_'))
        hdus.writeto(filename, clobber=True)
        
    @staticmethod
    def load(filename):
        hdus = fits.open(filename, memmap=False)
        hdr = hdus[0].header
        covtype = hdr['covtype']
        model = _GaussianMixtureModel(
            hdus['weights_'].data, hdus['means_'].data, hdus['covars_'].data, covtype)
        hdus.close()
        return model
    
    def sample(self, n_samples=1, random_state=None):
        
        if self.covtype != 'full':
            return NotImplementedError(
                'covariance type "{0}" not implemented yet.'.format(self.covtype))
        
        # Code adapted from sklearn's GMM.sample()
        if random_state is None:
            random_state = np.random.RandomState()

        weight_cdf = np.cumsum(self.weights_)
        X = np.empty((n_samples, self.n_dimensions))
        rand = random_state.rand(n_samples)
        # decide which component to use for each sample
        comps = weight_cdf.searchsorted(rand)
        # for each component, generate all needed samples
        for comp in range(self.n_components):
            # occurrences of current component in X
            comp_in_X = (comp == comps)
            # number of those occurrences
            num_comp_in_X = comp_in_X.sum()
            if num_comp_in_X > 0:
                X[comp_in_X] = random_state.multivariate_normal(
                    self.means_[comp], self.covars_[comp], num_comp_in_X)
        return X
    
    def sample_full_pdf(self, n_samples=1, random_state=None):
        ''''sample() uses the component datum x is closest too, sample_full_pdf() uses sum of components at datum x
        this is more time consuming than sample() and difference is negligible'''
        if self.covtype != 'full':
            return NotImplementedError(
                'covariance type "{0}" not implemented yet.'.format(self.covtype))

        from scipy.stats import multivariate_normal
        def get_mv(means_,covars_):
            mv=[]
            for mu, C in zip(means_, covars_):
                mv+= [ multivariate_normal(mean=mu, cov=C) ]
            return mv
                
        def prob_map(means_,covars_,weights_,\
                     xrng=(0.,1.),yrng=(0.,1.),npts=2**10):
            '''returns 
            -pmap: 2d probability map, with requirement that integral(pdf d2x) within 1% of 1
            -xvec,yvec: vectors where x[ix] and y[ix] data points have probability pmap[ix,iy]'''
            assert(xrng[1] > xrng[0] and yrng[1] > yrng[0])
            xstep= (xrng[1]-xrng[0])/float(npts-1)
            ystep= (yrng[1]-yrng[0])/float(npts-1)
            x,y  = np.mgrid[xrng[0]:xrng[1]+xstep:xstep, yrng[0]:yrng[1]+ystep:ystep]
            pos = np.empty(x.shape + (2,)) #npts x npts x 2
            pos[:, :, 0] = x; pos[:, :, 1] = y
            maps= np.zeros(x.shape+ (len(weights_),)) # x n_components
            # Multi-Variate function
            mv= get_mv(means_,covars_)
            # probability map for each component
            for dist,W,cnt in zip(mv,weights_, range(len(weights_))):
                maps[:,:,cnt]= dist.pdf(pos) * W
                print "map %d, dist.pdf max=%.2f, wt=%.3f" % (cnt,dist.pdf(pos).max(),W)
            # summed prob map
            pmap= np.sum(maps, axis=2) #some over components*weights
            xvec= x[:,0]
            yvec= y[0,:]
            # intregal of pdf over 2d map = 1
        #     assert( abs(1.- pmap.sum()*xstep*ystep) <= 0.01 )
            assert( np.diff(xvec).min() > 0. and np.diff(yvec).min() > 0.)
            return pmap, xvec,yvec
            
        # 2D probability map
        grrange = (-0.2, 2.0)
        rzrange = (-0.4, 2.5)
        pmap,xvec,yvec= prob_map(self.means_,self.covars_,self.weights_,\
                                 xrng=rzrange,yrng=grrange)
        # Sample self.n_dimensions using map
        if random_state is None:
            r = np.random.RandomState()
        # Make max pdf = 1 so always pick that cell
        pmap/= pmap.max()
        # Store x,y values to use
        X = np.empty((n_samples, self.n_dimensions))+np.nan
        # Get samples
        cnt=0
        # pick a random cell
        ix,iy=(r.rand(2)*len(xvec)).astype(int)
        # Random [0,1)
        likely= r.rand(1)
        while cnt < n_samples:
            if likely <= pmap[ix,iy]: # Sample it!
                X[cnt,:]= xvec[ix],yvec[iy]
                cnt+=1
            # Pick new cell in either case
            ix,iy=(r.rand(2)*len(xvec)).astype(int) 
            likely= r.rand(1)
        assert( np.where(~np.isfinite(X))[0].size == 0)
        return X


def add_MoG_curves(ax, means_, covars_, weights_):
    '''plot 2-sigma ellipses for each multivariate component'''
    ax.scatter(means_[:, 0], means_[:, 1], c='w')
    scale=2.
    for cnt, mu, C, w in zip(range(means_.shape[0]),means_, covars_, weights_):
    #     draw_ellipse(mu, C, scales=[1.5], ax=ax, fc='none', ec='k')
        # Draw MoG outlines
        sigma_x2 = C[0, 0]
        sigma_y2 = C[1, 1]
        sigma_xy = C[0, 1]

        alpha = 0.5 * np.arctan2(2 * sigma_xy,
                             (sigma_x2 - sigma_y2))
        tmp1 = 0.5 * (sigma_x2 + sigma_y2)
        tmp2 = np.sqrt(0.25 * (sigma_x2 - sigma_y2) ** 2 + sigma_xy ** 2)

        sigma1 = np.sqrt(tmp1 + tmp2)
        sigma2 = np.sqrt(tmp1 - tmp2)
        print('comp=%d, sigma1=%f,sigma2=%f' % (cnt+1,sigma1,sigma2))

        ax.text(mu[0],mu[1],str(cnt+1),color='blue')
        ax.add_patch(Ellipse((mu[0], mu[1]),
                     2 * scale * sigma1, 2 * scale * sigma2,
                     alpha * 180. / np.pi,\
                     fc='none', ec='k'))

def get_rgb_cols():
    return [(255,0,255),(102,255,255),(0,153,153),\
            (255,0,0),(0,255,0),(0,0,255),\
            (0,0,0)]

def flux2mag(nanoflux):
    return 22.5-2.5*np.log10(nanoflux)

# Globals
xyrange=dict(x_star=[-0.5,2.2],\
             y_star=[-0.3,2.],\
             x_elg=[-0.5,2.2],\
             y_elg=[-0.3,2.],\
             x_lrg= [0, 2.5],\
             y_lrg= [-2, 6])
def rm_last_ticklabel(ax):
    '''for multiplot'''
    labels=ax.get_xticks().tolist()
    labels=np.array(labels).astype(float) #prevent from making float
    labels=list(labels)
    labels[-1]=''
    ax.set_xticklabels(labels)



class TSBox(object):
    '''functions to add Target Selection box to ELG, LRG, etc plot
    add_ts_box -- main functino to call'''
    def __init__(self,src='ELG'):
        self.src=src

    def add_ts_box(self, ax, xlim=None,ylim=None):
        '''draw color selection box'''
        assert(xlim is not None and ylim is not None)
        if self.src == 'ELG':
            #g-r vs. r-z
            xint= newton(self.ts_root,np.array([1.]),args=('y1-y2',))
            
            x=np.linspace(xlim[0],xlim[1],num=1000)
            y=np.linspace(ylim[0],ylim[1],num=1000)
            x1,y1= x,self.ts_box(x,'y1')
            x2,y2= x,self.ts_box(x,'y2')
            x3,y3= np.array([0.3]*len(x)),y
            x4,y4= np.array([0.6]*len(x)),y
            b= np.all((x >= 0.3,x <= xint),axis=0)
            x1,y1= x1[b],y1[b]
            b= np.all((x >= xint,x <= 1.6),axis=0)
            x2,y2= x2[b],y2[b]
            b= y3 <= np.min(y1)
            x3,y3= x3[b],y3[b]
            b= y4 <= np.min(y2)
            x4,y4= x4[b],y4[b]
            ax.plot(x1,y1,'k--',lw=2)
            ax.plot(x2,y2,'k--',lw=2)
            ax.plot(x3,y3,'k--',lw=2)
            ax.plot(x4,y4,'k--',lw=2) 
        elif self.src == 'LRG':
            #r-w1 vs. r-z
            x=np.linspace(xlim[0],xlim[1],num=1000)
            y=np.linspace(ylim[0],ylim[1],num=1000)
            x1,y1= x,self.ts_box(x,'y1')
            x2,y2= np.array([1.5]*len(x)),y
            b= x >= 1.5
            x1,y1= x1[b],y1[b]
            b= y2 >= np.min(y1)
            x2,y2= x2[b],y2[b]
            ax.plot(x1,y1,'k--',lw=2)
            ax.plot(x2,y2,'k--',lw=2)
        else: raise ValueError('src=%s not supported' % src)

    def ts_box(self, x,name):
        if self.src == 'ELG':
            if name == 'y1': return 1.15*x-0.15
            elif name == 'y2': return -1.2*x+1.6
            else: raise ValueError
        elif self.src == 'LRG':
            if name == 'y1': return 1.8*x-1.
            else: raise ValueError
        else: raise ValueError('src=%s not supported' % self.src)

    def ts_root(self,x,name):
        if self.src == 'ELG':
            if name == 'y1-y2': return self.ts_box(x,'y1')-self.ts_box(x,'y2')
            else: raise ValueError
        else: raise ValueError('non ELG not supported')



def elg_data():
    '''Use DEEP2 ELGs whose SEDs have been modeled.'''
    elgs = fits.getdata('/project/projectdirs/desi/spectro/templates/basis_templates/v2.2/elg_templates_v2.0.fits', 1)
    # Colors
    gg = elgs['DECAM_G']
    rr = elgs['DECAM_R']
    zz = elgs['DECAM_Z']
    gr = gg - rr
    rz = rr - zz
    Xall = np.array([rz, gr]).T
    # Cuts
    has_morph = elgs['radius_halflight'] > 0
    cuts= dict(has_morph=has_morph) 
    print('%d/%d of Fit Template Spectra ELGs have morphologies' % (len(elgs[has_morph]), len(elgs)))
    # Morphology
    morph= {}
    morph['rz'] = rz[has_morph]
    morph['gr'] = gr[has_morph]
    morph['r50'] = elgs['RADIUS_HALFLIGHT'][has_morph] #arcsec
    morph['n'] = elgs['SERSICN'][has_morph]                            
    morph['ba'] = elgs['AXIS_RATIO'][has_morph] #minor/major
    return Xall,cuts,morph                            


class ReadWrite(object):
    def read_fits(self,fn):
        #return Table(fits.getdata(fn, 1))
        return fits_table(fn)

#DR=2,savefig=False,alpha=1.,\
#brick_primary=True,anymask=False,allmask=False,fracflux=False):
class CommonInit(ReadWrite):
    def __init__(self,**kwargs):
        # Syntax is self.val2 = kwargs.get('val2',"default value")
        self.DR= kwargs.get('DR',2)
        self.savefig= kwargs.get('savefig',False)
        self.alpha= kwargs.get('alpha',1.)
        self.brick_primary= kwargs.get('brick_primary',True)
        self.anymask= kwargs.get('anymask',False)
        self.allmask= kwargs.get('allmask',False)
        self.fracflux= kwargs.get('fracflux',False)
        print('self.fracflux=',self.fracflux,'kwargs= ',kwargs)
        if self.DR == 2:
            self.truth_dir= '/project/projectdirs/desi/target/analysis/truth'
        elif self.DR == 3:
            self.truth_dir= '/project/projectdirs/desi/users/burleigh/desi/target/analysis/truth'
        else: raise valueerror()

    def imaging_cut(self,data):
        '''data is a fits_table object with Tractor Catalogue columns'''
        cut=np.ones(len(data)).astype(bool)
        if self.brick_primary:
            if data.get('brick_primary').dtype == 'bool':
                cut*= data.get('brick_primary') == True
            elif data.get('brick_primary').dtype == 'S1':
                cut*= data.get('brick_primary') == 'T'
            else: 
                raise ValueError('brick_primary has type=',data.get('brick_primary').dtype)
        if self.anymask:
            cut*= np.all((data.get('decam_anymask')[:, [1,2,4]] == 0),axis=1)
        if self.allmask:
            cut*= np.all((data.get('decam_allmask')[:, [1,2,4]] == 0),axis=1)
        if self.fracflux:
            cut*= np.all((data.get('decam_fracflux')[:, [1,2,4]] < 0.05),axis=1)
        return cut

    def std_star_cut(self,data):
        '''See: https://desi.lbl.gov/trac/wiki/TargetSelectionWG/TargetSelection#SpectrophotometricStandardStarsFSTD
           data is a fits_table object with Tractor Catalogue columns
        '''
        # Remove whitespaces 'PSF ' --> 'PSF'
        if not 'decam_mag' in data.get_columns():
            from theValidator.catalogues import CatalogueFuncs 
            CatalogueFuncs().set_extra_data(data)
        RFLUX_obs = data.get('decam_flux')[:,2]
        GFLUX = data.get('decam_flux')[:,1] / data.get('decam_mw_transmission')[:,1]
        RFLUX = data.get('decam_flux')[:,2] / data.get('decam_mw_transmission')[:,2]
        ZFLUX = data.get('decam_flux')[:,4] / data.get('decam_mw_transmission')[:,4]
        GRZSN = data.get('decam_flux')[:,[1,2,4]] * np.sqrt(data.get('decam_flux_ivar')[:,[1,2,4]])
        GRCOLOR = 2.5 * np.log10(RFLUX / GFLUX)
        RZCOLOR = 2.5 * np.log10(ZFLUX / RFLUX)
        cut= np.all((data.get('brick_primary') == True,\
                     data.get('type') == 'PSF',\
                     np.all((data.get('decam_fracflux')[:, [1,2,4]] < 0.04),axis=1),\
                     np.all((GRZSN > 10),axis=1),\
                     RFLUX_obs < 10**((22.5-16.0)/2.5)),axis=0)
                     #np.power(GRCOLOR - 0.32,2) + np.power(RZCOLOR - 0.13,2) < 0.06**2,\
                     #RFLUX_obs > 10**((22.5-19.0)/2.5)),axis=0)
        return cut 


class ELG(CommonInit):
    def __init__(self,**kwargs):
        super(ELG, self).__init__(**kwargs)
        self.rlimit= 23.4
        if 'rlimit' in kwargs:
            self.rlimit= kwargs['rlimit']
        print('ELGs, self.rlimit= ',self.rlimit)
    
    def get_elgs_FDR_cuts(self):
        '''version 3.0 of data discussed in
        https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=912'''
        if self.DR == 2:
            zcat = self.read_fits(os.path.join(self.truth_dir,'../deep2/v3.0/','deep2-field1-oii.fits.gz'))
            G_MAG= zcat.get('cfhtls_g')
            R_MAG= zcat.get('cfhtls_r')
            Z_MAG= zcat.get('cfhtls_z')
            rmag_cut= R_MAG<self.rlimit 
        elif self.DR == 3:
            zcat = self.read_fits(os.path.join(self.truth_dir,'deep2f234-dr3matched.fits'))
            decals = self.read_fits(os.path.join(self.truth_dir,'dr3-deep2f234matched.fits'))
            # Add mag data 
            G_MAG= decals.get('decam_mag')[:,1]
            R_MAG= decals.get('decam_mag')[:,2]
            Z_MAG= decals.get('decam_mag')[:,4]
            rflux= decals.get('decam_flux')[:,2]/decals.get('decam_mw_transmission')[:,2]
            rmag_cut= rflux > 10**((22.5-self.rlimit)/2.5)
            rmag_cut*= self.imaging_cut(decals) 
        # Cuts
        oiicut1 = 8E-17 # [erg/s/cm2]
        zmin = 0.6
        loz = np.all((zcat.get('zhelio')<zmin,\
                      rmag_cut),axis=0)
        oiifaint = np.all((zcat.get('zhelio')>zmin,\
                           rmag_cut,\
                           zcat.get('oii_3727_err')!=-2.0,\
                           zcat.get('oii_3727')<oiicut1),axis=0)
        oiibright_loz = np.all((zcat.get('zhelio')>zmin,\
                                zcat.get('zhelio')<1.0,\
                                rmag_cut,\
                                zcat.get('oii_3727_err')!=-2.0,\
                                zcat.get('oii_3727')>oiicut1),axis=0)
        oiibright_hiz = np.all((zcat.get('zhelio')>1.0,\
                                rmag_cut,\
                                zcat.get('oii_3727_err')!=-2.0,\
                                zcat.get('oii_3727')>oiicut1),axis=0)
        med2hiz_oiibright= np.all((zcat.get('zhelio')>0.6,\
                                   rmag_cut,\
                                   zcat.get('oii_3727_err')!=-2.0,\
                                   zcat.get('oii_3727')>oiicut1),axis=0)

        # color data
        rz,gr,r={},{},{}
        for key,cut in zip(['loz','oiifaint','oiibright_loz','oiibright_hiz','med2hiz_oiibright'],\
                           [loz,oiifaint,oiibright_loz,oiibright_hiz,med2hiz_oiibright]):
            rz[key]= (R_MAG - Z_MAG)[cut]
            gr[key]= (G_MAG - R_MAG)[cut]
            r[key]= R_MAG[cut]
        return rz,gr,r

    def plot_FDR(self):
        rz,gr,r= self.get_elgs_FDR_cuts()
        # Object to add target selection box
        ts= TSBox(src='ELG')
        fig, ax = plt.subplots()
        # Add box
        xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
        ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
        # Add points
        b= 'loz'
        ax.scatter(rz[b],gr[b], marker='^', color='magenta', label=r'$z<0.6$')

        b='oiifaint'
        ax.scatter(rz[b],gr[b], marker='s', color='tan',
                        label=r'$z>0.6, [OII]<8\times10^{-17}$')

        b= 'oiibright_loz'
        ax.scatter(rz[b],gr[b], marker='o', color='powderblue',
                        label=r'$z>0.6, [OII]>8\times10^{-17}$')

        b='oiibright_hiz'
        ax.scatter(rz[b],gr[b], marker='o', color='blue',
                        label=r'$z>1.0, [OII]>8\times10^{-17}$')
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        xlab= ax.set_xlabel('r-z')
        ylab= ax.set_ylabel('g-r')
        leg=ax.legend(loc=(0,1.05), ncol=2,prop={'size': 14}, labelspacing=0.2,
                  markerscale=1.5)
        name='dr%d_FDR_ELG.png' % self.DR
        kwargs= dict(bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
        if self.savefig:
            plt.savefig(name, **kwargs)
            plt.close()
            print('Wrote {}'.format(name))
 
    def plot_FDR_multipanel(self):
        rz,gr,r= self.get_elgs_FDR_cuts()
        ts= TSBox(src='ELG')
        xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
        # Plot
        fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(18,4))
        plt.subplots_adjust(wspace=0.1,hspace=0)
        for cnt,key,col,marker,ti in zip(range(4),\
                               ['loz','oiifaint','oiibright_loz','oiibright_hiz'],\
                               ['magenta','tan','powderblue','blue'],\
                               ['^','s','o','o'],\
                               [r'$z<0.6$',r'$z>0.6, [OII]<8\times10^{-17}$',r'$z>0.6, [OII]>8\times10^{-17}$',r'$z>1.0, [OII]>8\times10^{-17}$']):
            # Add box
            ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
            # Add points
            b= key
            ax[cnt].scatter(rz[b],gr[b], marker=marker,color=col)
            ti_loc=ax[cnt].set_title(ti)
            ax[cnt].set_xlim(xrange)
            ax[cnt].set_ylim(yrange)
            xlab= ax[cnt].set_xlabel('r-z')
            ylab= ax[cnt].set_ylabel('g-r')
            name='dr%d_FDR_ELG_multipanel.png' % self.DR
        kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
        if self.savefig:
            plt.savefig(name, **kwargs)
            plt.close()
            print('Wrote {}'.format(name))

    def plot(self):
        self.plot_FDR()
        self.plot_FDR_multipanel()
        #plot_FDR(self.Xall,self.cuts,src='ELG')
        #b= self.cuts['any_elg']
        #color_color_plot(self.Xall[b,:],src='ELG',append='_FDR') #,extra=True)
        #Xall,cuts, morph= elg_data()
        #color_color_plot(Xall, src='ELG',append='_synth') #,extra=True)
        #b= cuts['has_morph']
        #color_color_plot(Xall[b,:],src='ELG',append='_synth+morph') #,extra=True)

class LRG(CommonInit):
    def __init__(self,**kwargs):
        super(LRG, self).__init__(**kwargs)
        self.zlimit= 20.46
        if 'zlimit' in  kwargs:
            self.zlimit= kwargs['zlimit']
        print('LRGs, self.zlimit= ',self.zlimit)
    
    def get_lrgs_FDR_cuts(self):
        if self.DR == 2:
            decals=self.read_fits( os.path.join(self.truth_dir,'decals-dr2-cosmos-zphot.fits.gz') )
            spec=self.read_fits( os.path.join(self.truth_dir,'cosmos-zphot.fits.gz') )
        elif self.DR == 3:
            decals=self.read_fits( os.path.join(self.truth_dir,'dr3-cosmoszphotmatched.fits') )
            spec=self.read_fits( os.path.join(self.truth_dir,'cosmos-zphot-dr3matched.fits') )
        # DECaLS
        Z_FLUX = decals.get('decam_flux'.lower())[:,4] / decals.get('decam_mw_transmission'.lower())[:,4]
        W1_FLUX = decals.get('wise_flux'.lower())[:,0] / decals.get('wise_mw_transmission'.lower())[:,0]
        index={}
        # BUG!!
        #index['decals']= np.all((Z_FLUX < 10**((22.5-20.46)/2.5),\
        index['decals']= np.all((Z_FLUX > 10**((22.5-self.zlimit)/2.5),\
                                 W1_FLUX > 0.),axis=0)
        index['decals']*= self.imaging_cut(decals)
        # Cosmos
        # http://irsa.ipac.caltech.edu/data/COSMOS/gator_docs/cosmos_zphot_mag25_colDescriptions.html
        # http://irsa.ipac.caltech.edu/data/COSMOS/tables/redshift/cosmos_zphot_mag25.README
        for key in ['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy']:
            if key == 'star': 
                index[key]= np.all((index['decals'],\
                                    spec.get('type') == 1),axis=0)
            elif key == 'blue_galaxy': 
                index[key]= np.all((index['decals'],\
                                    spec.get('type') == 0,\
                                    spec.get('mod_gal') > 8),axis=0)
            elif key == 'red_galaxy_lowz': 
                index[key]= np.all((index['decals'],\
                                    spec.get('type') == 0,\
                                    spec.get('mod_gal') <= 8,\
                                    spec.get('zp_gal') <= 0.6),axis=0)
            elif key == 'red_galaxy_hiz': 
                index[key]= np.all((index['decals'],\
                                    spec.get('type') == 0,\
                                    spec.get('mod_gal') <= 8,\
                                    spec.get('zp_gal') > 0.6),axis=0)
        # return Mags
        rz,rW1,r={},{},{}
        R_FLUX = decals.get('decam_flux')[:,2] / decals.get('decam_mw_transmission')[:,2]
        for key in ['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy']:
            cut= index[key]
            rz[key]= flux2mag(R_FLUX[cut]) - flux2mag(Z_FLUX[cut])
            rW1[key]= flux2mag(R_FLUX[cut]) - flux2mag(W1_FLUX[cut])
            r[key]= flux2mag(R_FLUX[cut])
        return rz,rW1,r
                                 
    def plot_FDR(self):
        rz,rW1,r= self.get_lrgs_FDR_cuts()
        fig,ax = plt.subplots()
        rgb_cols=get_rgb_cols()
        xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
        # Add box
        ts= TSBox(src='LRG')
        ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
        # Data
        for cnt,key,rgb in zip(range(4),rz.keys(),rgb_cols):
            rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
            ax.scatter(rz[key],rW1[key],c=[rgb],edgecolors='none',marker='o',s=10.,rasterized=True,label=key)
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        xlab=ax.set_xlabel('r-z')
        ylab=ax.set_ylabel('r-W1')
        leg=ax.legend(loc=(0,1.05), ncol=2,prop={'size': 14}, labelspacing=0.2,\
                      markerscale=2,scatterpoints=1)
        #handles,labels = ax.get_legend_handles_labels()
        #index=[0,1,2,3]
        #handles,labels= np.array(handles)[index],np.array(labels)[index]
        #leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
        name='dr%d_FDR_LRG.png' % self.DR
        if self.savefig:
            plt.savefig(name,\
                        bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote {}'.format(name))

    def plot_FDR_multipanel(self):
        rz,rW1,r= self.get_lrgs_FDR_cuts()
        fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(16,4))
        plt.subplots_adjust(wspace=0.1,hspace=0)
        rgb_cols=get_rgb_cols()
        for cnt,key,rgb in zip(range(4),rz.keys(),rgb_cols):
            rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
            ax[cnt].scatter(rz[key],rW1[key],c=[rgb],edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
            ti=ax[cnt].set_title(key)
            ax[cnt].set_xlim([0,2.5])
            ax[cnt].set_ylim([-2,6])
            xlab=ax[cnt].set_xlabel('r-z')
            # Add box
            ts= TSBox(src='LRG')
            xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
            ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
            ylab=ax[cnt].set_ylabel('r-W1')
        #handles,labels = ax.get_legend_handles_labels()
        #index=[0,1,2,3]
        #handles,labels= np.array(handles)[index],np.array(labels)[index]
        #leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
        name='dr%d_FDR_LRG_multi.png' % self.DR
        if self.savefig:
            plt.savefig(name,\
                        bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote {}'.format(name))

    def get_vipers(self):
        '''LRGs from VIPERS in CFHTLS W4 field (12 deg2)'''
        # DR2 matched
        if self.DR == 2:
            decals = self.read_fits(os.path.join(self.truth_dir,'decals-dr2-vipers-w4.fits.gz'))
            vip = self.read_fits(os.path.join(self.truth_dir,'vipers-w4.fits.gz'))
        elif self.DR == 3:
            decals = self.read_fits(os.path.join(self.truth_dir,'dr3-vipersw1w4matched.fits'))
            vip = self.read_fits(os.path.join(self.truth_dir,'vipersw1w4-dr3matched.fits'))
        Z_FLUX = decals.get('decam_flux')[:,4] / decals.get('decam_mw_transmission')[:,4]
        W1_FLUX = decals.get('wise_flux')[:,0] / decals.get('wise_mw_transmission')[:,0]
        index={}
        index['decals']= np.all((Z_FLUX > 10**((22.5-self.zlimit)/2.5),\
                                 W1_FLUX > 0.),axis=0)
        index['decals']*= self.imaging_cut(decals)
        # VIPERS
        # https://arxiv.org/abs/1310.1008
        # https://arxiv.org/abs/1303.2623
        flag= vip.get('zflg').astype(int)
        index['good_z']= np.all((flag >= 2,\
                                 flag <= 9,\
                                 vip.get('zspec') < 9.9),axis=0) 
        # return Mags
        rz,rW1={},{}
        R_FLUX = decals.get('decam_flux')[:,2] / decals.get('decam_mw_transmission')[:,2]
        cut= np.all((index['decals'],\
                     index['good_z']),axis=0)
        rz= flux2mag(R_FLUX[cut]) - flux2mag(Z_FLUX[cut])
        rW1= flux2mag(R_FLUX[cut]) - flux2mag(W1_FLUX[cut])
        return rz,rW1
    
    def plot_vipers(self):
        rz,rW1= self.get_vipers()
        fig,ax = plt.subplots(figsize=(5,4))
        plt.subplots_adjust(wspace=0,hspace=0)
        rgb=get_rgb_cols()[0]
        rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
        ax.scatter(rz,rW1,c=[rgb],edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
        ax.set_xlim([0,2.5])
        ax.set_ylim([-2,6])
        xlab=ax.set_xlabel('r-z')
        ylab=ax.set_ylabel('r-W1')
        # Add box
        ts= TSBox(src='LRG')
        xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
        ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
        #handles,labels = ax.get_legend_handles_labels()
        #index=[0,1,2,3]
        #handles,labels= np.array(handles)[index],np.array(labels)[index]
        #leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
        name='dr%d_LRG_vipers.png' % self.DR
        if self.savefig:
            plt.savefig(name,\
                        bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote {}'.format(name))


    def plot(self):
        self.plot_FDR()
        self.plot_FDR_multipanel()
        self.plot_vipers()
        #plot_FDR(self.Xall,self.cuts,src='LRG')
        #color_color_plot(self.Xall,src='LRG',append='cc') #,extra=True)
        #b= self.cuts['lrg')
        #color_color_plot(self.Xall[b,:],src='LRG') #,extra=True)

class STAR(CommonInit):
    def __init__(self,**kwargs):
        super(STAR,self).__init__(**kwargs)
    
    def get_sweepstars(self):
        '''Model the g-r, r-z color-color sequence for stars'''
        # Build a sample of stars with good photometry from a single sweep.
        rbright = 18
        rfaint = 19.5
        swp_dir='/global/project/projectdirs/cosmo/data/legacysurvey/dr2/sweep/2.0'
        sweep = self.read_fits(os.path.join(swp_dir,'sweep-340p000-350p005.fits'))
        keep = np.where((sweep.get('type') == 'PSF ')*
                        (np.sum((sweep.get('decam_flux')[:, [1,2,4]] > 0)*1, axis=1)==3)*
                        (np.sum((sweep.get('DECAM_ANYMASK')[:, [1,2,4]] > 0)*1, axis=1)==0)*
                        (np.sum((sweep.get('DECAM_FRACFLUX')[:, [1,2,4]] < 0.05)*1, axis=1)==3)*
                        (sweep.get('decam_flux')[:,2]<(10**(0.4*(22.5-rbright))))*
                        (sweep.get('decam_flux')[:,2]>(10**(0.4*(22.5-rfaint)))))[0]
        stars = sweep[keep]
        print('dr2stars sample: {}'.format(len(stars)))
        gg = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 1])
        rr = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 2])
        zz = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 4])
        gr = gg - rr
        rz = rr - zz
        return np.array(rz),np.array(gr)

    def plot_sweepstars(self): 
        rz,gr= self.get_sweepstars()
        fig,ax = plt.subplots(figsize=(5,4))
        plt.subplots_adjust(wspace=0,hspace=0)
        rgb=get_rgb_cols()[0]
        rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
        ax.scatter(rz,gr,c=[rgb],edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
        ax.set_xlim([-0.5,2.2])
        ax.set_ylim([-0.3,2.])
        xlab=ax.set_xlabel('r-z')
        ylab=ax.set_ylabel('g-r')
        #handles,labels = ax.get_legend_handles_labels()
        #index=[0,1,2,3]
        #handles,labels= np.array(handles)[index],np.array(labels)[index]
        #leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
        name='STAR_dr2.png'
        if self.savefig:
            plt.savefig(name,\
                        bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote {}'.format(name))
      
    def get_purestars(self):
        # https://desi.lbl.gov/trac/wiki/TargetSelectionWG/TargetSelection#SpectrophotometricStandardStarsFSTD
        if self.DR == 2:
            stars=fits_table(os.path.join(self.truth_dir,'Stars_str82_355_4.DECaLS.dr2.fits'))
            stars.cut( self.std_star_cut(stars) )
            return stars
        elif self.DR == 3:
            raise ValueError()
 
    def plot(self,savefig):
        self.plot_sweepstars() 
        #plt.plot(self.cat.get('ra'),self.cat.get('dec'))
        #plt.savefig('test.png')
        #plt.close()

class QSO(CommonInit):
    def __init__(self,**kwargs):
        super(QSO,self).__init__(**kwargs)
        self.rlimit= 22.7
        if 'rlimit' in  kwargs:
            self.rlimit= kwargs['rlimit']
        print('QSOs, self.rlimit= ',self.rlimit)
  
    def get_qsos(self):
        from theValidator.catalogues import CatalogueFuncs 
        if self.DR == 2:        
            qsos= self.read_fits( os.path.join(self.truth_dir,'AllQSO.DECaLS.dr2.fits') )
            # Add AB mags
            CatalogueFuncs().set_extra_data(qsos)
            qsos.cut( self.imaging_cut(qsos) )
            # r < 22.7, grz > 17
            GFLUX = qsos.get('decam_flux')[:,1] / qsos.get('decam_mw_transmission')[:,1]
            RFLUX = qsos.get('decam_flux')[:,2] / qsos.get('decam_mw_transmission')[:,2]
            ZFLUX = qsos.get('decam_flux')[:,4] / qsos.get('decam_mw_transmission')[:,4]
            GRZFLUX = (GFLUX + 0.8* RFLUX + 0.5* ZFLUX ) / 2.3
            cut= np.all((RFLUX > 10**((22.5-self.rlimit)/2.5),\
                         GRZFLUX < 10**((22.5-17.0)/2.5)),axis=0)
            qsos.cut(cut)
            return qsos
        elif self.DR == 3:
            qsos=self.read_fits( os.path.join(self.truth_dir,'qso-dr3sweepmatched.fits') )
            decals=self.read_fits( os.path.join(self.truth_dir,'dr3-qsosweepmatched.fits') )
            decals.set('z',qsos.get('z'))
            CatalogueFuncs().set_extra_data(decals)
            decals.cuts( self.imaging_cut(decals) )
            return decals
        #G_FLUX= decals.get('decam_flux')[:,1]/decals.get('decam_mw_transmission')[:,1]
        #R_FLUX= decals.get('decam_flux')[:,2]/decals.get('decam_mw_transmission')[:,2]
        #Z_FLUX= decals.get('decam_flux')[:,4]/decals.get('decam_mw_transmission')[:,4]
        # DECaLS
        #index={}
        #rfaint=22.7
        #grzbright=17.
        #index['decals']= np.all((R_FLUX > 10**((22.5-rfaint)/2.5),\
        #                         G_FLUX < 10**((22.5-grzbright)/2.5),\
        #                         R_FLUX < 10**((22.5-grzbright)/2.5),\
        #                         Z_FLUX < 10**((22.5-grzbright)/2.5),\
        #                         decals.get('brick_primary') == True),axis=0)
        # QSO
        # Return
     
    def plot_FDR(self):
        # Data
        qsos= self.get_qsos()
        star_obj= STAR(DR=self.DR,savefig=False)
        stars= star_obj.get_purestars()
        hiz=2.1
        index={}
        index['hiz']= qsos.get('z') > hiz
        index['loz']= qsos.get('z') <= hiz
        # Plot
        fig,ax = plt.subplots(1,2,figsize=(10,4))
        plt.subplots_adjust(wspace=0.1,hspace=0)
        # Stars
        ax[0].scatter(stars.get('decam_mag')[:,2]-stars.get('decam_mag')[:,4],\
                      stars.get('decam_mag')[:,1]-stars.get('decam_mag')[:,2],\
                      c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
        W= 0.75*stars.get('wise_mag')[:,0]+ 0.25*stars.get('wise_mag')[:,1]
        ax[1].scatter(stars.get('decam_mag')[:,1]-stars.get('decam_mag')[:,4],\
                      stars.get('decam_mag')[:,2]-W,\
                      c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
        # QSOs
        for key,lab,col in zip(['loz','hiz'],['(z < 2.1)','(z > 2.1)'],['magenta','red']):
            i= index[key]
            ax[0].scatter(qsos.get('decam_mag')[:,2][i]-qsos.get('decam_mag')[:,4][i],\
                          qsos.get('decam_mag')[:,1][i]-qsos.get('decam_mag')[:,2][i],\
                          c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
            W= 0.75*qsos.get('wise_mag')[:,0]+ 0.25*qsos.get('wise_mag')[:,1]
            ax[1].scatter(qsos.get('decam_mag')[:,1][i]-qsos.get('decam_mag')[:,4][i],\
                          qsos.get('decam_mag')[:,2][i]-W[i],\
                          c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
        
        #for xlim,ylim,x_lab,y_lab in ax[0].set_xlim([-0.5,3.])
        ax[0].set_xlim([-0.5,3.])
        ax[1].set_xlim([-0.5,4.5])
        ax[0].set_ylim([-0.5,2.5])
        ax[1].set_ylim([-2.5,3.5])
        xlab=ax[0].set_xlabel('r-z')
        xlab=ax[1].set_xlabel('g-z')
        ylab=ax[0].set_ylabel('g-r')
        ylab=ax[1].set_ylabel('r-W')
        leg=ax[0].legend(loc=(0,1.02),scatterpoints=1,ncol=3,markerscale=2)
        ## Add box
        #ts= TSBox(src='LRG')
        #xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
        #ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
        name='dr%d_FDR_QSO.png' % self.DR
        if self.savefig:
            plt.savefig(name,\
                        bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote {}'.format(name))

 
    def plot(self):
        self.plot_FDR()


class GalaxyPrior(object):
    def __init__(self):
        #self.qso= QSO()
        self.star= STAR()
        self.lrg= LRG()
        self.elg= ELG()
    def plot_all(self):
        #self.qso.plot()
        self.star.plot()
        self.lrg.plot()
        self.elg.plot()

if __name__ == '__main__':
    gals=GalaxyPrior()
    gals.plot_all()
    print "gals.__dict__= ",gals.__dict__