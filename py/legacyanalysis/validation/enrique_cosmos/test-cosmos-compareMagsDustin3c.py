import astropy.io.fits as fits
import numpy as n
import glob
import os
from os.path import join
import sys
import matplotlib.pyplot as p
import statsmodels.api as sm
from scipy.spatial import KDTree
import scipy.stats as st
from scipy.optimize import curve_fit as cu

import numpy as np


gfun = lambda x, m0, s0 : st.norm.pdf(x,loc=m0,scale=s0)
f=open("depth-comparison.txt","w")

def compareMags(N1 = "40", N2 = "41"):
    hdu = fits.open("catalog-"+N1+".fits")
    dat = hdu[1].data
    noJunk3 = (dat['brick_primary']) & (dat['decam_anymask'].T[1]==0) & (dat['decam_anymask'].T[2]==0) & (dat['decam_anymask'].T[4]==0) 
    dr3 = dat[noJunk3]
    areaDR3 = ( n.max(dat['ra']) - n.min(dat['ra']) ) * ( n.max(dat['dec']) - n.min(dat['dec']) )*n.cos( n.pi * n.mean(dat['dec']) / 180. )

    hdu = fits.open("catalog-"+N2+".fits")
    dat = hdu[1].data
    noJunk2 = (dat['brick_primary']) & (dat['decam_anymask'].T[1]==0) & (dat['decam_anymask'].T[2]==0) & (dat['decam_anymask'].T[4]==0) 
    dr2 = dat[noJunk2]
    areaDR2 = ( n.max(dat['ra']) - n.min(dat['ra']) ) * ( n.max(dat['dec']) - n.min(dat['dec']) )*n.cos( n.pi * n.mean(dat['dec']) / 180. )

    dr2T = KDTree(n.transpose([dr2['ra'],dr2['dec']]))
    dr3T = KDTree(n.transpose([dr3['ra'],dr3['dec']]))

    # look for the nearest neighbour from dr3 sources in dr2
    ds, ids = dr2T.query(dr3T.data)

    dsmax = 1/3600.

    sel = (ds<dsmax)
    print len(dr3[sel]) # 116,791 withins 1 arcsec.
    print len(dr2[ids[sel]])

    #maglimits = dict(g=(21., 24.), r=(21., 23.5), z=(21., 22.5))
    # full-depth
    maglimits = dict(g=(21.5, 24.), r=(21.5, 24.0), z=(21.5, 23))

    p.clf()
    # Check whether fitting the histogram introduces bias (binning is sinning)...
    # not much.
    G = np.random.normal(size=10000)
    nnn,bbb, ppp=p.hist(G, bins=n.arange(-6,6,0.2), histtype='step', normed=True)
    out = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
    print('Fit to Gaussian random draws:', out)
    p.clf()

    f.write(N1 + " " + N2 + '\n')

    bands= {"g","r","z"}

    for band in bands:
        wb ={"g":1,"r":2,"z":4}    
        w = int(wb[band])
        print "band,w=", band,w,N1,N2
        mag_dr2 = 22.5 - 2.5 * n.log10(dr2[ids[sel]]['decam_flux'].T[w])
        mag_dr3 = 22.5 - 2.5 * n.log10(dr3[sel]['decam_flux'].T[w])
        mag_mean = 22.5 - 2.5 * n.log10((dr2[ids[sel]]['decam_flux'].T[w] + dr3[sel]['decam_flux'].T[w])/2.)
    
        iv_dr2 = dr2[ids[sel]]['decam_flux_ivar'].T[w]
        iv_dr3 = dr3[sel]['decam_flux_ivar'].T[w]
    
        df = dr3[sel]['decam_flux'].T[w] - dr2[ids[sel]]['decam_flux'].T[w]
        sigma = (1./dr2[ids[sel]]['decam_flux_ivar'].T[w] + 1./dr3[sel]['decam_flux_ivar'].T[w])**(0.5)

        p.figure(2,(5,5))
        p.axes([0.17,0.15,0.75,0.75])
        p.plot(n.arange(-6,6,0.2), st.norm.pdf(n.arange(-6,6,0.2),loc=0,scale=1.), 'k--', lw=0.5, label='N(0,1)')
        maglo,maghi = maglimits[band]

        ok = ((mag_mean > maglo) & (mag_mean < maghi) &
        (iv_dr2 > 0) & (iv_dr3 > 0))

        wb ={"g":1,"r":2,"z":4}    
        w = int(wb[band])
        xwb= wb
        del xwb[band]
        xbands=xwb.keys()
        print band,xbands
        for xband in xbands:
            wx= int(xwb[xband])
            mag_mean2 = 22.5 - 2.5 * n.log10((dr2[ids[sel]]['decam_flux'].T[wx] + dr3[sel]['decam_flux'].T[wx])/2.)
            maglo2,maghi2 = maglimits[xband]
            ok = (ok) & (mag_mean2 > maglo2) & (mag_mean2 < maghi2)

        
        sig = df[ok] / sigma[ok]
        mean_df=sum(sig)/len(sig)
        sigma_df=(sum((sig-mean_df)**2)/len(sig))**(0.5)

        xcor = n.zeros(len(xwb))        
        ib= 0
        for xband in xbands[0]:  # do only one color
            wx= int(xwb[xband])
            dfx = dr3[sel]['decam_flux'].T[wx] - dr2[ids[sel]]['decam_flux'].T[wx]
            sigmax = (1./dr2[ids[sel]]['decam_flux_ivar'].T[wx] + 1./dr3[sel]['decam_flux_ivar'].T[wx])**(0.5)
            sigc= (dfx[ok]-df[ok])/(sigmax[ok]**2+sigma[ok]**2)**0.5
            sigx = dfx[ok] / sigmax[ok]   # current value
            meanx_df=sum(sigx)/len(sigx)
            sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
            xcor[ib]=sum((sig-mean_df)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df
            ib=ib+1

        sig = sigc
        mean_df=sum(sig)/len(sig)
        sigma_df=(sum((sig-mean_df)**2)/len(sig))**(0.5)
            
        xbandst=" Cov["+band+"-"+str(xbands[0])+","+band+"-"+str(xbands[1])+"]"
        txtr1=xbandst+"=["+str(n.round(xcor[0],2))+","+str(n.round(xcor[1],2))+"]"
        print txtr1
        out1= (mean_df,sigma_df)
        out3 = (n.median(sig), (n.percentile(sig,84) - n.percentile(sig,16))/2.)
        txt1="$\sigma$="+str(n.round(out1[1],2))+" , "+str(n.round(out3[1],2))

        ##nnn,bbb, ppp=p.hist(sig, bins=n.arange(-6,6,0.2),histtype='step',label="All "+txt1+" "+txtr1, normed=True)
        nnn,bbb, ppp=p.hist(sig, bins=n.arange(-6,6,0.2),histtype='step',label="All "+txt1+" "+txtr1, normed=True)
        out = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
        print  str(out[0]),txt1
        #gfit=st.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out[0][0],scale=out[0][1])
        #gfit1=gfit[(abs(bbb)<2)]
        #nnn1=nnn[(abs(bbb)<2)]
        #chi2=sum((nnn1-gfit1)**2/nnn1)*len(sig) # error is sqrt(N)/len(sig)
        #print "chi2=",chi2,chi2/(len(nnn1)-3)
        kst=sm.stats.lillifors(sig) # https://en.wikipedia.org/wiki/Lilliefors_test
        
        p.plot(n.arange(-6,6,0.2), st.norm.pdf(n.arange(-6,6,0.2),loc=out[0][0],scale=out[0][1]), 'b--', lw=2, label='All N fit [mean,$\sigma$]='+str(n.round(out[0],2))+' logP='+str(n.round(n.log10(kst[1]),1)))

        #ok2 = (ok) & (dr2[ids[sel]]['type']=="PSF") & (dr3[sel]['type']=='PSF')
        ok2 = (ok) & (dr2[ids[sel]]['type']=="SIMP") & (dr3[sel]['type']=='SIMP')

        sig2 = df[ok2] / sigma[ok2]
        mean_df2=sum(sig2)/len(sig2)
        sigma_df2=(sum((sig2-mean_df2)**2)/len(sig2))**(0.5)
        #sigma2_df=(n.var(sig))**(0.5)
        ib= 0
        for xband in xbands[0]:  # do only one color            
            wx= int(xwb[xband])
            dfx = dr3[sel]['decam_flux'].T[wx] - dr2[ids[sel]]['decam_flux'].T[wx]
            sigmax = (1./dr2[ids[sel]]['decam_flux_ivar'].T[wx] + 1./dr3[sel]['decam_flux_ivar'].T[wx])**(0.5)
            sigx = dfx[ok2] / sigmax[ok2]
            sig2c= (dfx[ok2]-df[ok2])/(sigmax[ok2]**2+sigma[ok2]**2)**0.5
            meanx_df=sum(sigx)/len(sigx)
            sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
            xcor[ib]=sum((sig2-mean_df2)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df2
            ib=ib+1
        txtr2=xbandst+"=["+str(n.round(xcor[0],2))+","+str(n.round(xcor[1],2))+"]"
        print txtr2

        sig2 = sig2c
        mean_df2=sum(sig2)/len(sig2)
        sigma_df2=(sum((sig2-mean_df2)**2)/len(sig2))**(0.5)

        out4= (mean_df2,sigma_df2)
        out5 = (n.median(sig2c), (n.percentile(sig2c,84) - n.percentile(sig2,16))/2.)
        #txt2="$\sigma$="+str(n.round(out4[1],2))+" s68="+str(n.round(out5[1],2))
        txt2="$\sigma$="+str(n.round(out4[1],2))+" , "+str(n.round(out5[1],2))

        
        nnn,bbb,ppp = p.hist(sig2, bins=n.arange(-6,6,0.2), histtype='step', label='SIMP '+txt2+" "+txtr2, normed=True)
        out2 = cu(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
        print str(out2[0]),txt2

        #gfit=st.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out2[0][0],scale=out2[0][1])
        #gfit2=gfit[(abs(bbb)<2)]
        #nnn2=nnn[(abs(bbb)<2)]
        #chi22=sum((nnn2-gfit2)**2/nnn2)*len(sig2) # error is 1/sqrt(N)
        #print "chi22=",chi22,chi22/(len(nnn2)-3)     
        kst2=sm.stats.lillifors(sig2) #https://en.wikipedia.org/wiki/Lilliefors_test

        
        p.plot(n.arange(-6,6,0.2), st.norm.pdf(n.arange(-6,6,0.2),loc=out2[0][0],scale=out2[0][1]), 'g--', lw=2, label='SIMP N fit [mean,$\sigma$]='+str(n.round(out2[0],2))+' logP='+str(n.round(n.log10(kst2[1]),1)))

        print  str(out2[0]),out4,out5

        f.write(band + str(out[0])+str(out1)+str(out3)+str(out2[0])+str(out4)+str(out5))
        f.write('\n')
        
        #ok = (dr2[ids[sel]]['type']=="PSF")&(dr3[sel]['decam_psfsize'].T[1]<1.5)
        #p.hist(df[ok] / sigma[ok], bins=n.arange(-6,6,0.1), weights = n.ones_like(df[ok])/areaDR3, histtype='step', label='type PSF, PSFsize<1.5', normed=True)
        p.xlabel("[ ("+band+"-"+str(xbands[0])+") "+"("+N1+")-"+"("+N2+") ] / sqrt[ Var("+N2+")+Var("+N1+") ]")
        p.ylabel('Normed counts %g < %s < %g' %(maglo,band,maghi))
        p.xlim((-4,4))
        p.ylim((0,0.7))
        gp = p.legend(loc=2, fontsize=10)
        gp.set_frame_on(False)
        p.title('Cosmos Colors ('+band+"-"+str(xbands[0])+") "+N1+"-"+N2)
        p.grid()
        p.savefig(join("plots","cosmos-comparison-depth-normed-"+band+"-v"+N1+"-"+N2+"c.png"))
        p.clf()
        

        

# compareMags(N1 = "30", N2 = "31")
# compareMags(N1 = "31", N2 = "32")
# compareMags(N1 = "32", N2 = "30")

compareMags(N1 = "40", N2 = "41")
compareMags(N1 = "41", N2 = "42")
compareMags(N1 = "42", N2 = "40")

# compareMags(N1 = "10", N2 = "11")
# compareMags(N1 = "11", N2 = "12")
# compareMags(N1 = "12", N2 = "10")
# 
# compareMags(N1 = "20", N2 = "21")
# compareMags(N1 = "21", N2 = "22")
# compareMags(N1 = "22", N2 = "20")
