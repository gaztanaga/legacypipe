
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as st
from scipy.optimize import curve_fit as cu

gfun = lambda x, m0, s0 : st.norm.pdf(x,loc=m0,scale=s0)
f=open("depth-comparison-new.txt","w")

def noJunk(cat):
    '''cat is a fits_table object
    keep only S/N>5 '''
    i = (cat.get('brick_primary')) & (cat.get('decam_anymask')[:,1]==0) & \
        (cat.get('decam_anymask')[:,2]==0) & (cat.get('decam_anymask')[:,4]==0) & \
        (cat.get('decam_flux')[:,1]>0)  &  (cat.get('decam_flux')[:,2]>0) &  \
        (cat.get('decam_flux')[:,4]>0) & \
        (cat.get('decam_flux_ivar')[:,1]*(cat.get('decam_flux')[:,1])**2>25)  & \
        (cat.get('decam_flux_ivar')[:,2]*(cat.get('decam_flux')[:,2])**2>25)  & \
        (cat.get('decam_flux_ivar')[:,4]*(cat.get('decam_flux')[:,4])**2>25)
    return i

def areaCat(cat):
    '''cat is a fits_table object'''
    area = ( np.max(cat.get('ra')) - np.min(cat.get('ra')) ) * \
           ( np.max(cat.get('dec')) - np.min(cat.get('dec')) )*\
           np.cos( np.pi * np.mean(cat.get('dec')) / 180. )
    return area

def compareMags(dr3,dr2,\
               N1='ref',N2='obs',savefig=False):
    '''Enriques cosmso script used cosmo40/ as reference and cosmos41 as observed
    so N1=40 and N2=41
    catalague matching to the ref,
    he also called these dr3 for cosmos40, dr2 for cosmos41'''
    i_junk= noJunk(dr3)
    i_junk*= noJunk(dr2)
    dr2= dr2[i_junk]
    dr3= dr3[i_junk]
    areaDR3 = areaCat(dr3)
    areaDR2 = areaCat(dr2)

    #maglimits = dict(g=(21., 24.), r=(21., 23.5), z=(21., 22.5))
    # full-depth
    #maglimits = dict(g=(21.5, 24.), r=(21.5, 24.0), z=(21.5, 23))
    maglimits = dict(g=(23, 24.), r=(22.4, 23.4), z=(21.5, 22.5))

    # Check whether fitting the histogram introduces bias (binning is sinning)...
    # not much.
    G = np.random.normal(size=10000)
    nnn,bbb, ppp=plt.hist(G, bins=np.arange(-6,6,0.2), histtype='step', normed=True)
    out = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
    print('Fit to Gaussian random draws:', out)
    plt.clf()

    f.write(N1 + " " + N2 + '\n')

    bands= {"g","r","z"}

    for band in bands:
        wb ={"g":1,"r":2,"z":4}
        w = int(wb[band])
        print "band,w=", band,w,N1,N2
        mag_dr2 = dr2.get('decam_mag')
        mag_dr3 = dr3.get('decam_mag')
        mag_mean = 22.5 - 2.5 * np.log10((dr2.get('decam_flux')[:,w] + dr3.get('decam_flux')[:,w])/2.)

        iv_dr2 = dr2.get('decam_flux_ivar')[:,w]
        iv_dr3 = dr3.get('decam_flux_ivar')[:,w]

        df = dr3.get('decam_flux')[:,w] - dr2.get('decam_flux')[:,w]
        sigma = (1./dr2.get('decam_flux_ivar')[:,w] + 1./dr3.get('decam_flux_ivar')[:,w])**(0.5)

        plt.figure(2,(5,5))
        plt.axes([0.17,0.15,0.75,0.75])
        plt.plot(np.arange(-6,6,0.2), st.norm.pdf(np.arange(-6,6,0.2),loc=0,scale=1.), 'k--', lw=0.5, label='N(0,1)')
        maglo,maghi = maglimits[band]

        ok = ((mag_mean > maglo) & (mag_mean < maghi) &
        (iv_dr2 > 0) & (iv_dr3 > 0))

        wb ={"g":1,"r":2,"z":4}
        w = int(wb[band])
        xwb= wb
        del xwb[band]
        xbands=xwb.keys()
        for xband in xbands:   # apply magnitude cuts in all bands
            wx= int(xwb[xband])
            mag_mean2 = 22.5 - 2.5 * np.log10((dr2.get('decam_flux')[:,wx] + \
                                               dr3.get('decam_flux')[:,wx])/2.)
            maglo2,maghi2 = maglimits[xband]
            ok = (ok) & (mag_mean2 > maglo2) & (mag_mean2 < maghi2)
        
        sig = df[ok] / sigma[ok]
        mean_df=sum(sig)/len(sig)
        sigma_df=(sum((sig-mean_df)**2)/len(sig))**(0.5)

        xcor = np.zeros(len(xwb))
        ib= 0
        for xband in xbands:  # estimate covariance in xcor
            wx= int(xwb[xband])          
            dfx = dr3.get('decam_flux')[:,wx] - dr2.get('decam_flux')[:,wx]
            sigmax = (1./dr2.get('decam_flux_ivar')[:,wx] + 1./dr3.get('decam_flux_ivar')[:,wx])**(0.5)
            sigx = dfx[ok] / sigmax[ok]
            meanx_df=sum(sigx)/len(sigx)
            sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
            xcor[ib]=sum((sig-mean_df)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df
            ib=ib+1

        xbandst=" Cov["+band+"-"+str(xbands[0])+","+band+"-"+str(xbands[1])+"]"
        txtr1=xbandst+"=["+str(np.round(xcor[0],2))+","+str(np.round(xcor[1],2))+"]"
        print txtr1
        out1= (mean_df,sigma_df)
        out3 = (np.median(sig), (np.percentile(sig,84) - np.percentile(sig,16))/2.)
        txt1="$\sigma$="+str(np.round(out1[1],2))+" , "+str(np.round(out3[1],2))

        nnn,bbb, ppp=plt.hist(sig, bins=np.arange(-6,6,0.2),histtype='step',label="All "+txt1+" "+txtr1, normed=True)
        out = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
        print  str(out[0]),txt1
        #gfit=st.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out[0][0],scale=out[0][1])
        #gfit1=gfit[(abs(bbb)<2)]
        #nnn1=nnn[(abs(bbb)<2)]
        #chi2=sum((nnn1-gfit1)**2/nnn1)*len(sig) # error is sqrt(N)/len(sig)
        #print "chi2=",chi2,chi2/(len(nnn1)-3)
        kst=sm.stats.lillifors(sig) # https://en.wikipedia.org/wiki/Lilliefors_test

        plt.plot(np.arange(-6,6,0.2), st.norm.pdf(np.arange(-6,6,0.2),loc=out[0][0],scale=out[0][1]), 'b--', lw=2, label='All N fit [mean,$\sigma$]='+str(np.round(out[0],2))+' logP='+str(np.round(np.log10(kst[1]),1)))

        #ok2 = (ok) & (dr2.get('type')=="PSF") & (dr3.get('type')=='PSF')
        ok2 = (ok) & (dr2.get('type')=="SIMP") & (dr3.get('type')=='SIMP')

        sig2 = df[ok2] / sigma[ok2]
        mean_df2=sum(sig2)/len(sig2)
        sigma_df2=(sum((sig2-mean_df2)**2)/len(sig2))**(0.5)
        #sigma2_df=(np.var(sig))**(0.5)
        ib= 0
        for xband in xbands:
            wx= int(xwb[xband])
            dfx = dr3.get('decam_flux')[:,wx] - dr2.get('decam_flux')[:,wx]
            sigmax = (1./dr2.get('decam_flux_ivar')[:,wx] + 1./dr3.get('decam_flux_ivar')[:,wx])**(0.5)
            sigx = dfx[ok2] / sigmax[ok2]
            meanx_df=sum(sigx)/len(sigx)
            sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
            xcor[ib]=sum((sig2-mean_df2)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df2
            ib=ib+1
        txtr2=xbandst+"=["+str(np.round(xcor[0],2))+","+str(np.round(xcor[1],2))+"]"
        print txtr2

        out4= (mean_df2,sigma_df2)
        out5 = (np.median(sig2), (np.percentile(sig2,84) - np.percentile(sig2,16))/2.)
        #txt2="$\sigma$="+str(np.round(out4[1],2))+" s68="+str(np.round(out5[1],2))
        txt2="$\sigma$="+str(np.round(out4[1],2))+" , "+str(np.round(out5[1],2))
        
        nnn,bbb,ppp = plt.hist(sig2, bins=np.arange(-6,6,0.2), histtype='step', label='SIMP '+txt2+" "+txtr2, normed=True)
        out2 = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
        print str(out2[0]),txt2

        #gfit=st.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out2[0][0],scale=out2[0][1])
        #gfit2=gfit[(abs(bbb)<2)]
        #nnn2=nnn[(abs(bbb)<2)]
        #chi22=sum((nnn2-gfit2)**2/nnn2)*len(sig2) # error is 1/sqrt(N)
        #print "chi22=",chi22,chi22/(len(nnn2)-3)     
        kst2=sm.stats.lillifors(sig2) #https://en.wikipedia.org/wiki/Lilliefors_test

        plt.plot(np.arange(-6,6,0.2), st.norm.pdf(np.arange(-6,6,0.2),\
                 loc=out2[0][0],scale=out2[0][1]), 'g--', lw=2, \
                 label=r'SIMP N fit [mean,$\sigma$]='+\
                       str(np.round(out2[0],2))+\
                       ' logP='+str(np.round(np.log10(kst2[1]),1)))

        print  str(out2[0]),out4,out5

        f.write(band + str(out[0])+str(out1)+str(out3)+str(out2[0])+str(out4)+str(out5))
        f.write('\n')

        #ok = (dr2.get('type')=="PSF")&(dr3['decam_psfsize'].T[1]<1.5)
        #plt.hist(df[ok] / sigma[ok], bins=np.arange(-6,6,0.1), weights = np.ones_like(df[ok])/areaDR3, histtype='step', label='type PSF, PSFsize<1.5', normed=True)
        plt.xlabel("[ "+band+"("+N1+")-"+band+"("+N2+") ] / sqrt[ Var("+N2+")+Var("+N1+") ]")
        plt.ylabel('Normed counts %g < %s < %g' %(maglo,band,maghi))
        plt.xlim((-4,4))
        plt.ylim((0,0.7))
        gp = plt.legend(loc=2, fontsize=10)
        gp.set_frame_on(False)
        plt.title('Cosmos '+band+N1+"-"+N2+" #= "+str(len(sig))+'(All) '+str(len(sig2)))
        plt.grid()
        plt.savefig("enrique-cosmos-"+band+"-v"+N1+"-"+N2+".png")
        plt.clf()
