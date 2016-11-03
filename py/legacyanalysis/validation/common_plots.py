'''plotting functions for Validation
input -- one or more Single_TractorCat() objects
      -- ref is reference Single_TractorCat()
         test is test ...
'''
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg') #display backend
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from scipy import stats as sp_stats

# Globals
class PlotKwargs(object):
    def __init__(self):
        self.ax= dict(fontweight='bold',fontsize='medium')
        self.text=dict(fontweight='bold',fontsize='medium',va='top',ha='left')
        self.leg=dict(frameon=True,fontsize='x-small')
        self.save= dict(bbox_inches='tight',dpi=150)

kwargs= PlotKwargs()
##########

# Helpful functions for making plots
def bin_up(data_bin_by,data_for_percentile, bin_minmax=(18.,26.),nbins=20):
    '''bins "data_for_percentile" into "nbins" using "data_bin_by" to decide how indices are assigned to bins
    returns bin center,N,q25,50,75 for each bin
    '''
    bin_edges= np.linspace(bin_minmax[0],bin_minmax[1],num= nbins+1)
    vals={}
    for key in ['q50','q25','q75','n']: vals[key]=np.zeros(nbins)+np.nan
    vals['binc']= (bin_edges[1:]+bin_edges[:-1])/2.
    for i,low,hi in zip(range(nbins), bin_edges[:-1],bin_edges[1:]):
        keep= np.all((low < data_bin_by,data_bin_by <= hi),axis=0)
        if np.where(keep)[0].size > 0:
            vals['n'][i]= np.where(keep)[0].size
            vals['q25'][i]= np.percentile(data_for_percentile[keep],q=25)
            vals['q50'][i]= np.percentile(data_for_percentile[keep],q=50)
            vals['q75'][i]= np.percentile(data_for_percentile[keep],q=75)
        else:
            vals['n'][i]=0 
    return vals

# Main plotting functions 
def nobs(tractor, outname='test.png',savefig=False):
    '''make histograms of nobs so can compare depths of g,r,z between the two catalogues
    tractor -- Tractor catalogue in a table'''   
    hi= np.max(tractor.get('decam_nobs')[:,[1,2,4]])
    fig,ax= plt.subplots(3,1)
    for i, band,iband in zip(range(3),['g','r','z'],[1,2,4]):
        ax[i].hist(tractor.get('decam_nobs')[:,iband],\
                   bins=hi+1,normed=True,cumulative=True,align='mid')
        xlab=ax[i].set_xlabel('nobs %s' % band, **kwargs.ax)
        ylab=ax[i].set_ylabel('CDF', **kwargs.ax)
    if savefig == True:
        plt.savefig(outname, bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()

def radec(tractor,outname='test.png',savefig=False): 
    '''ra,dec distribution of objects
    obj -- Single_TractorCat()'''
    plt.scatter(tractor.get('ra'), tractor.get('dec'), \
                edgecolor='b',c='none',lw=1.)
    xlab=plt.xlabel('RA', **kwargs.ax)
    ylab=plt.ylabel('DEC', **kwargs.ax)
    if savefig == True:
        plt.savefig(outname,bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()


def hist_types(obj, name='',savefig=False):
    '''number of psf,exp,dev,comp, etc
    obj -- Single_TractorCat()'''
    types= ['PSF','SIMP','EXP','DEV','COMP']
    # the x locations for the groups
    ind = np.arange(len(types))  
    # the width of the bars
    width = 0.35       
    ht= np.zeros(len(types),dtype=int)
    for cnt,typ in enumerate(types):
        # Mask to type desired
        ht[cnt]= obj.number_not_masked(['current',typ.lower()])
    # Plot
    fig, ax = plt.subplots()
    rects = ax.bar(ind, ht, width, color='b')
    ylab= ax.set_ylabel("counts")
    ax.set_xticks(ind + width)
    ax.set_xticklabels(types)
    if savefig == True:
        plt.savefig(os.path.join(obj.outdir,'hist_types_%s.png' % name), \
                    bbox_extra_artists=[ylab], **kwargs.save)
        plt.close()


def sn_vs_mag(tractor, mag_minmax=(18.,26.),name='',savefig=False):
    '''plots Signal to Noise vs. mag for each band'''
    min,max= mag_minmax
    # Bin up SN values
    bin_SN={}
    for band,iband in zip(['g','r','z'],[1,2,4]):
        bin_SN[band]= bin_up(tractor.get('decam_mag')[:,iband], \
                       tractor.get('decam_flux')[:,iband]*np.sqrt(tractor.get('decam_flux_ivar')[:,iband]),\
                       bin_minmax=mag_minmax)
    #setup plot
    fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
    plt.subplots_adjust(wspace=0.25)
    #plot SN
    for cnt,band,color in zip(range(3),['g','r','z'],['g','r','m']):
        #horiz line at SN = 5
        ax[cnt].plot([mag_minmax[0],mag_minmax[1]],[5,5],'k--',lw=2)
        #data
        ax[cnt].plot(bin_SN[band]['binc'], bin_SN[band]['q50'],c=color,ls='-',lw=2)
        ax[cnt].fill_between(bin_SN[band]['binc'],bin_SN[band]['q25'],bin_SN[band]['q75'],color=color,alpha=0.25)
    #labels
    for cnt,band in zip(range(3),['g','r','z']):
        ax[cnt].set_yscale('log')
        xlab=ax[cnt].set_xlabel('%s' % band, **kwargs.ax)
        ax[cnt].set_ylim(1,100)
        ax[cnt].set_xlim(mag_minmax)
    ylab=ax[0].set_ylabel('S/N', **kwargs.ax)
    ax[2].text(26,5,'S/N = 5  ',**kwargs.text)
    if savefig == True:
        plt.savefig(os.path.join(obj.outdir,'sn_%s.png' % name), \
                    bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()

def create_confusion_matrix(ref_tractor,test_tractor):
    '''compares MATCHED reference (truth) to test (prediction)
    ref_obj,test_obj -- reference,test Single_TractorCat()
    return 5x5 confusion matrix and colum/row names'''
    cm=np.zeros((5,5))-1
    types=['PSF','SIMP','EXP','DEV','COMP']
    for i_ref,ref_type in enumerate(types):
        cut_ref= np.where(ref_tractor.get('type') == ref_type)[0]
        #n_ref= ref_obj.number_not_masked(['current',ref_type.lower()])
        for i_test,test_type in enumerate(types):
            n_test= np.where(test_tractor.get('type')[ cut_ref ] == test_type)[0].size
            if cut_ref.size > 0: cm[i_ref,i_test]= float(n_test)/cut_ref.size
            else: cm[i_ref,i_test]= np.nan
    return cm,types


def confusion_matrix(ref_tractor,test_tractor, outname='test.png',savefig=False,\
                     ref_name='ref',test_name='test'):
    '''plot confusion matrix
    ref_obj,test_obj -- reference,test Single_TractorCat()'''
    cm,ticknames= create_confusion_matrix(ref_tractor,test_tractor)
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues,vmin=0,vmax=1)
    cbar=plt.colorbar()
    plt.xticks(range(len(ticknames)), ticknames)
    plt.yticks(range(len(ticknames)), ticknames)
    ylab=plt.ylabel('True (%s)' % ref_name, **kwargs.ax)
    xlab=plt.xlabel('Predicted (%s)' % test_name, **kwargs.ax)
    for row in range(len(ticknames)):
        for col in range(len(ticknames)):
            if np.isnan(cm[row,col]):
                plt.text(col,row,'n/a',va='center',ha='center')
            elif cm[row,col] > 0.5:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
            else:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
    if savefig == True:
        plt.savefig(outname,bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()

def create_stack(answer_type,predict_type, types=['PSF','SIMP','EXP','DEV','COMP'],slim=True):
    '''compares classifications of matched objects, returns 2D array which is conf matrix and xylabels
    return 5x5 confusion matrix and colum/row names
   answer_type,predict_type -- arrays of same length with reference and prediction types'''
    for typ in set(answer_type): assert(typ in types)
    for typ in set(predict_type): assert(typ in types)
    # if a type was not in answer (training) list then don't put in cm
    if slim: ans_types= set(answer_type)
    # put in cm regardless
    else: ans_types= set(types)
    cm=np.zeros((len(ans_types),len(types)))-1
    for i_ans,ans_type in enumerate(ans_types):
        ind= np.where(answer_type == ans_type)[0]
        for i_pred,pred_type in enumerate(types):
            n_pred= np.where(predict_type[ind] == pred_type)[0].size
            if ind.size > 0: cm[i_ans,i_pred]= float(n_pred)/ind.size # ind.size is constant for loop over pred_types
            else: cm[i_ans,i_pred]= np.nan
    if slim: return cm,ans_types,types #size ans_types != types
    else: return cm,types

def plot_stack(cm_stack,stack_names,all_names, \
               ref_name='ref',test_name='test',\
               outname='test.png',savefig=False):
    '''cm_stack -- list of single row confusion matrices
    stack_names -- list of same len as cm_stack, names for each row of cm_stack'''
    # combine list into single cm
    cm=np.zeros((len(cm_stack),len(all_names)))+np.nan
    for i in range(cm.shape[0]): cm[i,:]= cm_stack[i]
    # make usual cm, but labels repositioned
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    cbar=plt.colorbar()
    plt.xticks(range(len(all_names)), all_names)
    plt.yticks(range(len(stack_names)), stack_names)
    ylab=plt.ylabel('True (%s) = PSF' % ref_name)
    xlab=plt.xlabel('Predicted (%s)' % test_name)
    for row in range(len(stack_names)):
        for col in range(len(all_names)):
            if np.isnan(cm[row,col]):
                plt.text(col,row,'n/a',va='center',ha='center')
            elif cm[row,col] > 0.5:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
            else:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
            #if np.isnan(cm[row,col]): 
            #    plt.text(col,row,'n/a',va='center',ha='center')
            #else: plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center')
    if savefig == True:
        plt.savefig(outname, bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()

def stacked_confusion_matrix(ref_tractor,test_tractor,\
                             ref_name='ref',test_name='test', \
                             outname='test.png',savefig=False):
    cm_stack,stack_names=[],[]
    rbins= np.array([18.,20.,22.,23.,24.])
    for rmin,rmax in zip(rbins[:-1],rbins[1:]):
        # master cut
        br_cut= np.all((ref_tractor.get('decam_mag')[:,2] > rmin,ref_tractor.get('decam_mag')[:,2] <= rmax),axis=0)
        stack_names+= ["%d < r <= %d" % (int(rmin),int(rmax))]
        cm,ans_names,all_names= create_stack(np.array(['PSF']*len(ref_tractor)),
                                             test_tractor.get('type'))
        cm_stack+= [cm]
    plot_stack(cm_stack, stack_names,all_names, \
               ref_name=ref_name,test_name=test_name,outname=outname,savefig=savefig)


def matched_dist(obj,dist, name=''):
    '''dist -- array of distances in degress between matched objects'''
    pixscale=dict(decam=0.25,bass=0.45)
    # Plot
    fig,ax=plt.subplots()
    ax.hist(dist*3600,bins=50,color='b',align='mid')
    ax2 = ax.twiny()
    ax2.hist(dist*3600./pixscale['bass'],bins=50,color='g',align='mid',\
                visible=False)
    xlab= ax.set_xlabel("arcsec")
    xlab= ax2.set_xlabel("pixels [BASS]")
    ylab= ax.set_ylabel("Matched")
    plt.savefig(os.path.join(obj.outdir,"separation_hist_%s.png" % name), \
                bbox_extra_artists=[xlab,ylab], **kwargs.save)
    plt.close()

def chi_v_gaussian(ref_tractor,test_tractor,\
                   low=-8.,hi=8., outname='test.png',savefig=False):
    # Compute Chi
    chi={} 
    for band,iband in zip(['g','r','z'],[1,2,4]):
        chi[band]= (ref_tractor.get('decam_flux')[:,iband]-test_tractor.get('decam_flux')[:,iband])/\
                   np.sqrt( np.power(ref_tractor.get('decam_flux_ivar')[:,iband],-1)+\
                            np.power(test_tractor.get('decam_flux_ivar')[:,iband],-1))
    for b_low,b_hi in zip([18,19,20,21,22,23],[19,20,21,22,23,24]):
        #loop over mag bins, one 3 panel for each mag bin
        hist= dict(g=0,r=0,z=0)
        binc= dict(g=0,r=0,z=0)
        stats=dict(g=0,r=0,z=0)
        # Counts per bin
        for band,iband in zip(['g','r','z'],[1,2,4]):
            imag= np.all((b_low <= ref_tractor.get('decam_mag')[:,iband],\
                          ref_tractor.get('decam_mag')[:,iband] < b_hi),axis=0)
            hist[band],bins= np.histogram(chi[band][imag],\
                                    range=(low,hi),bins=50,normed=True)
            db= (bins[1:]-bins[:-1])/2
            binc[band]= bins[:-1]+db
        # Unit gaussian N(0,1)
        G= sp_stats.norm(0,1)
        xvals= np.linspace(low,hi)
        # Plot for each mag range
        fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
        plt.subplots_adjust(wspace=0.25)
        for cnt,band in zip(range(3),['g','r','z']):
            ax[cnt].step(binc[band],hist[band], where='mid',c='b',lw=2)
            ax[cnt].plot(xvals,G.pdf(xvals),c='g',label=r'$N(0,1)$')
        #labels
        for cnt,band in zip(range(3),['g','r','z']):
            if band == 'r': xlab=ax[cnt].set_xlabel(r'%s  $(F_{d}-F_{bm})/\sqrt{\sigma^2_{d}+\sigma^2_{bm}}$' % band, **kwargs.ax)
            else: xlab=ax[cnt].set_xlabel('%s' % band, **kwargs.ax)
            ax[cnt].set_ylim(0,0.6)
            ax[cnt].set_xlim(low,hi)
            ti=ax[cnt].set_title("%.1f < %s < %.1f" % (b_low,band,b_hi),**kwargs.ax)
        ylab=ax[0].set_ylabel('PDF', **kwargs.ax)
        ax[0].legend(loc='upper left',fontsize='medium')
        # Need unique name
        name=os.path.basename(outname).replace('.ipynb','')+'_%d-%d.png' % (b_low,b_hi) 
        if savefig == True:
            plt.savefig(os.path.join(os.path.dirname(outname),name), \
                        bbox_extra_artists=[ti,xlab,ylab], **kwargs.save)
            plt.close()

def delta_mag_vs_mag(ref_tractor,test_tractor, ref_name='ref',test_name='test',\
                     ylim=None,outname='test.png',savefig=False):
    fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
    plt.subplots_adjust(wspace=0.25)
    for cnt,iband in zip(range(3),[1,2,4]):
        delta= test_tractor.get('decam_mag')[:,iband]- ref_tractor.get('decam_mag')[:,iband]
        ax[cnt].scatter(ref_tractor.get('decam_mag')[:,iband],delta/ref_tractor.get('decam_mag')[:,iband],\
                        c='b',edgecolor='b',s=5) #,c='none',lw=2.)
    for cnt,band in zip(range(3),['g','r','z']):
        xlab=ax[cnt].set_xlabel('%s AB' % band, **kwargs.ax)
        if ylim is None:
            ax[cnt].set_ylim(-0.1,0.1)
        else: ax[cnt].set_ylim(ylim)
        ax[cnt].set_xlim(18,26)
    ylab=ax[0].set_ylabel('mag (%s) - mag(%s)' % (test_name,ref_name), **kwargs.ax)
    if savefig == True:
        plt.savefig(outname, bbox_extra_artists=[xlab,ylab], **kwargs.save)
        plt.close()


#def n_per_deg2(obj,deg2=1., req_mags=[24.,23.4,22.5],name=''):
#    '''compute number density in each bin for each band mag [18,requirement]
#    deg2 -- square degrees spanned by sources in obj table
#    req_mags -- image requirements grz<=24,23.4,22.5'''
#    bin_nd={}
#    for band,iband,req in zip(['g','r','z'],[1,2,4],req_mags):
#        bin_nd[band]={}
#        bins= np.linspace(18.,req,num=15)
#        bin_nd[band]['cnt'],junk= np.histogram(obj.t['decam_mag'][:,iband], bins=bins)
#        bin_nd[band]['binc']= (bins[1:]+bins[:-1])/2.
#        # bins and junk should be identical arrays
#        assert( np.all(np.all((bins,junk),axis=0)) )
#    # Plot
#    fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
#    plt.subplots_adjust(wspace=0.25)
#    for cnt,band in zip(range(3),['g','r','z']):
#        ax[cnt].step(bin_nd[band]['binc'],bin_nd[band]['cnt']/deg2, \
#                     where='mid',c='b',lw=2)
#    #labels
#    for cnt,band in zip(range(3),['g','r','z']):
#        xlab=ax[cnt].set_xlabel('%s' % band, **kwargs.ax) 
#    ylab=ax[0].set_ylabel('counts/deg2', **kwargs.ax)
#    # Make space for and rotate the x-axis tick labels
#    fig.autofmt_xdate()
#    plt.savefig(os.path.join(obj.outdir,'n_per_deg2_%s.png' % name), \
#                bbox_extra_artists=[xlab,ylab], **kwargs.save)
#    plt.close()



