
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from tractor.brightness import NanoMaggies
from catalogues import Cuts4MatchedCats

class Dustins(object):
    def __init__(self,ref_cat,obs_cat,imatch,d2d,\
                 ref_name='ref',obs_name='obs',savefig=False):
        #plot_all(self,ref_cat,obs_cat,
        #        ref_name='ref',obs_name='obs')
        self.outdir='./'
        self.cuts= Cuts4MatchedCats(ref_cat[ imatch['ref']], obs_cat[imatch['obs']])
        # Plot
        self.match_distance(d2d,savefig=savefig)
        self.fluxflux(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.fluxdeltaflux(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.grz_ratio_fluxerr(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.magmag(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.magdeltamag(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.deltamag_err(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.deltamag_errbars(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.stephist(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)
        self.curvehist(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig)

    def match_distance(self,dist,prefix='',savefig=False):
        plt.hist(dist * 3600., 100,range=(0,1))
        plt.xlabel('Match distance (arcsec)')
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smatch_dist.png' % prefix))
            plt.close()

    def fluxflux(self,matched1,matched2,\
                 ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            K = np.flatnonzero(self.cuts.good[band]) 

            print('Median mw_trans', band, 'is',
                  np.median(matched1.decam_mw_transmission[:,iband]))
            ax[cnt].errorbar(matched1.decam_flux[K,iband],
                         matched2.decam_flux[K,iband],
                         fmt='.', color=cc,
                         xerr=1./np.sqrt(matched1.decam_flux_ivar[K,iband]),
                         yerr=1./np.sqrt(matched2.decam_flux_ivar[K,iband]),
                         alpha=0.1,
                         )
            ax[cnt].set_xlabel('%s flux: %s' % (ref_name, band))
            ax[cnt].set_ylabel('%s flux: %s' % (obs_name, band))
            ax[cnt].plot([-1e6, 1e6], [-1e6,1e6], 'k-', alpha=1.)
            ax[cnt].axis([-100, 1000, -100, 1000])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sfluxflux_%s.png' % (prefix,band)))
            plt.close()

    def fluxdeltaflux(self,matched1,matched2,\
                     ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)  

            mag1, magerr1 = matched1.decam_mag[:,iband],1./np.sqrt(matched1.decam_mag_ivar[:,iband])

            iv1 = matched1.decam_flux_ivar[:, iband]
            iv2 = matched2.decam_flux_ivar[:, iband]
            std = np.sqrt(1./iv1 + 1./iv2)

            ax[cnt].plot(mag1[K],
                     (matched2.decam_flux[K,iband] - matched1.decam_flux[K,iband]) / std[K],
                     '.', alpha=0.1, color=cc)
            ax[cnt].plot(mag1[P],
                     (matched2.decam_flux[P,iband] - matched1.decam_flux[P,iband]) / std[P],
                     '.', alpha=0.1, color='k')
            ax[cnt].set_ylabel('(%s - %s) flux / flux errors (sigma): %s' % \
                               (obs_name, ref_name, band))
            ax[cnt].set_xlabel('%s mag: %s' % (ref_name, band))
            ax[cnt].axhline(0, color='k', alpha=0.5)
            ax[cnt].axis([24, 16, -10, 10])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sfluxdeltaflux_%s.png' % (prefix,band)))
            plt.close()

    def grz_ratio_fluxerr(self,matched1,matched2,\
                          ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        lp,lt = [],[]
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:        
            mag1, magerr1 = matched1.decam_mag[:,iband],1./np.sqrt(matched1.decam_mag_ivar[:,iband])

            iv1 = matched1.decam_flux_ivar[:, iband]
            iv2 = matched2.decam_flux_ivar[:, iband]
            std = np.sqrt(1./iv1 + 1./iv2)
            #std = np.hypot(std, 0.01)
            G= np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2 *\
                              np.isfinite(mag1) *\
                              (mag1 >= 20) * (mag1 < dict(g=24, r=23.5, z=22.5)[band]))

            n,b,p = ax[cnt].hist((matched2.decam_flux[G,iband] -
                              matched1.decam_flux[G,iband]) / std[G],
                     range=(-4, 4), bins=50, histtype='step', color=cc,
                     normed=True)

            sig = (matched2.decam_flux[G,iband] -
                   matched1.decam_flux[G,iband]) / std[G]
            print('Raw mean and std of points:', np.mean(sig), np.std(sig))
            med = np.median(sig)
            print('sig= ',sig,'len(sig)=',len(sig))
            if len(sig) > 0:
                rsigma = (np.percentile(sig, 84) - np.percentile(sig, 16)) / 2.
            else: rsigma=-1
            print('Median and percentile-based sigma:', med, rsigma)
            lp.append(p[0])
            lt.append('%s: %.2f +- %.2f' % (band, med, rsigma))

        bins = []
        gaussint = []
        for blo,bhi in zip(b, b[1:]):
            c = scipy.stats.norm.cdf(bhi) - scipy.stats.norm.cdf(blo)
            c /= (bhi - blo)
            #bins.extend([blo,bhi])
            #gaussint.extend([c,c])
            bins.append((blo+bhi)/2.)
            gaussint.append(c)
        ax[cnt].plot(bins, gaussint, 'k-', lw=2, alpha=0.5)
        ax[cnt].set_xlabel('Flux difference / error (sigma)')
        ax[cnt].axvline(0, color='k', alpha=0.1)
        ax[cnt].set_ylim(0, 0.45)
        ax[cnt].legend(lp, lt, loc='upper right')
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sgrz_ratio_fluxerr.png' % prefix))
            plt.close()

    def magmag(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * cuts.psf1 * cuts.psf2)

            ax[cnt].errorbar(mag1[K], mag2[K], fmt='.', color=cc,
                         xerr=magerr1[K], yerr=magerr2[K], alpha=0.1)
            ax[cnt].plot(mag1[P], mag2[P], 'k.', alpha=0.5)
            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('%s %s (mag)' % (obs_name, band))
            ax[cnt].plot([-1e6, 1e6], [-1e6,1e6], 'k-', alpha=1.)
            ax[cnt].axis([24, 16, 24, 16])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smagmag_%s.png' % (prefix,band)))
            plt.close()

    def magdeltamag(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * cuts.psf1 * cuts.psf2)

            ax[cnt].errorbar(mag1[K], mag2[K] - mag1[K], fmt='.', color=cc,
                         xerr=magerr1[K], yerr=magerr2[K], alpha=0.1)
            ax[cnt].plot(mag1[P], mag2[P] - mag1[P], 'k.', alpha=0.5)
            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('%s %s - %s %s (mag)' % (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)
            ax[cnt].axis([24, 16, -1, 1])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smagdeltamag_%s.png' % (prefix,band)))
            plt.close()

    def deltamag_err(self,matched1,matched2,\
                     ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            magbins = np.arange(16, 24.001, 0.5)
            ax[cnt].plot(mag1[K], (mag2[K]-mag1[K]) / np.hypot(magerr1[K], magerr2[K]),
                         '.', color=cc, alpha=0.1)
            ax[cnt].plot(mag1[P], (mag2[P]-mag1[P]) / np.hypot(magerr1[P], magerr2[P]),
                         'k.', alpha=0.5)

            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('(%s %s - %s %s) / errors (sigma)' %
                       (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)
            ax[cnt].axis([24, 16, -10, 10])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sdeltamag_err_%s.png' % (prefix,band)))
            plt.close()

    def deltamag_errbars(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * cuts.psf1 * cuts.psf2)

            magbins = np.arange(16, 24.001, 0.5)
            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
            ax[cnt].plot(meanmag[P], y[P], 'k.', alpha=0.1)

            midmag = []
            vals = np.zeros((len(magbins)-1, 5))
            median_err1 = []

            iqd_gauss = scipy.stats.norm.ppf(0.75) - scipy.stats.norm.ppf(0.25)

            # FIXME -- should we do some stats after taking off the mean difference?

            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
                I = P[(meanmag[P] >= mlo) * (meanmag[P] < mhi)]
                midmag.append((mlo+mhi)/2.)
                median_err1.append(np.median(magerr1[I]))
                if len(I) == 0:
                    continue
                # median and +- 1 sigma quantiles
                ybin = y[I]
                vals[bini,0] = np.percentile(ybin, 16)
                vals[bini,1] = np.median(ybin)
                vals[bini,2] = np.percentile(ybin, 84)
                # +- 2 sigma quantiles
                vals[bini,3] = np.percentile(ybin, 2.3)
                vals[bini,4] = np.percentile(ybin, 97.7)

                iqd = np.percentile(ybin, 75) - np.percentile(ybin, 25)

                print('Mag bin', midmag[-1], ': IQD is factor', iqd / iqd_gauss,
                      'vs expected for Gaussian;', len(ybin), 'points')

                # if iqd > iqd_gauss:
                #     # What error adding in quadrature would you need to make the IQD match?
                #     err = median_err1[-1]
                #     target_err = err * (iqd / iqd_gauss)
                #     sys_err = np.sqrt(target_err**2 - err**2)
                #     print('--> add systematic error', sys_err)

            # ~ Johan's cuts
            mlo = 21.
            mhi = dict(g=24., r=23.5, z=22.5)[band]
            I = P[(meanmag[P] >= mlo) * (meanmag[P] < mhi)]
            print('y=',y)
            print('I=',I)
            ybin = y[I]
            print('ybin =',ybin)
            if len(ybin) > 0:
                iqd = np.percentile(ybin, 75) - np.percentile(ybin, 25)
            else: iqd=-1
            print('Mag bin', mlo, mhi, 'band', band, ': IQD is factor',
                  iqd / iqd_gauss, 'vs expected for Gaussian;', len(ybin), 'points')
            if iqd > iqd_gauss:
                # What error adding in quadrature would you need to make
                # the IQD match?
                err = np.median(np.hypot(magerr1[I], magerr2[I]))
                print('Median error (hypot):', err)
                target_err = err * (iqd / iqd_gauss)
                print('Target:', target_err)
                sys_err = np.sqrt((target_err**2 - err**2) / 2.)
                print('--> add systematic error', sys_err)

                # check...
                err_sys = np.hypot(np.hypot(magerr1, sys_err),
                                   np.hypot(magerr2, sys_err))
                ysys = (mag2 - mag1) / err_sys
                ysys = ysys[I]
                print('Resulting median error:', np.median(err_sys[I]))
                iqd_sys = np.percentile(ysys, 75) - np.percentile(ysys, 25)
                print('--> IQD', iqd_sys / iqd_gauss, 'vs Gaussian')
                # Hmmm, this doesn't work... totally overshoots.


            ax[cnt].errorbar(midmag, vals[:,1], fmt='o', color='b',
                         yerr=(vals[:,1]-vals[:,0], vals[:,2]-vals[:,1]),
                         capthick=3, zorder=20)
            ax[cnt].errorbar(midmag, vals[:,1], fmt='o', color='b',
                         yerr=(vals[:,1]-vals[:,3], vals[:,4]-vals[:,1]),
                         capthick=2, zorder=20)
            ax[cnt].axhline( 1., color='b', alpha=0.2)
            ax[cnt].axhline(-1., color='b', alpha=0.2)
            ax[cnt].axhline( 2., color='b', alpha=0.2)
            ax[cnt].axhline(-2., color='b', alpha=0.2)

            for mag,err,y in zip(midmag, median_err1, vals[:,3]):
                if not np.isfinite(err):
                    continue
                if y < -6:
                    continue
                plt.text(mag, y-0.1, '%.3f' % err, va='top', ha='center', color='k',
                         fontsize=10)

            ax[cnt].set_xlabel('(%s + %s)/2 %s (mag), PSFs' % (ref_name, obs_name, band))
            ax[cnt].set_ylabel('(%s %s - %s %s) / errors (sigma)' %
                       (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)

            ax[cnt].axvline(21, color='k', alpha=0.3)
            ax[cnt].axvline(dict(g=24, r=23.5, z=22.5)[band], color='k', alpha=0.3)

            ax[cnt].axis([24.1, 16, -6, 6])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sdeltamag_errbars_%s.png' % (prefix,band)))
            plt.close()

    def stephist(self,matched1,matched2,\
                 ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            #magbins = np.append([16, 18], np.arange(20, 24.001, 0.5))
            if band == 'g':
                magbins = [20, 24]
            elif band == 'r':
                magbins = [20, 23.5]
            elif band == 'z':
                magbins = [20, 22.5]

            slo,shi = -5,5
            ha = dict(bins=25, range=(slo,shi), histtype='step', normed=True)
            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
            midmag = []
            nn = []
            rgbs = []
            lt,lp = [],[]
            I_empty=True
            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
                I = P[(mag1[P] >= mlo) * (mag1[P] < mhi)]
                if len(I) == 0:
                    continue
                I_empty=False
                ybin = y[I]
                rgb = [0.,0.,0.]
                rgb[0] = float(bini) / (len(magbins)-1)
                rgb[2] = 1. - rgb[0]
                n,b,p = plt.hist(ybin, color=rgb, **ha)
                lt.append('mag %g to %g' % (mlo,mhi))
                lp.append(p[0])
                midmag.append((mlo+mhi)/2.)
                nn.append(n)
                rgbs.append(rgb)

            bins = []
            gaussint = []
            if I_empty == False:
                for blo,bhi in zip(b, b[1:]):
                    #midbin.append((blo+bhi)/2.)
                    #gaussint.append(scipy.stats.norm.cdf(bhi) -
                    #                scipy.stats.norm.cdf(blo))
                    c = scipy.stats.norm.cdf(bhi) - scipy.stats.norm.cdf(blo)
                    c /= (bhi - blo)
                    bins.extend([blo,bhi])
                    gaussint.extend([c,c])
            ax[cnt].plot(bins, gaussint, 'k-', lw=2, alpha=0.5)

            ax[cnt].legend(lp, lt)
            ax[cnt].set_xlim(slo,shi)
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sstephist_%s.png' % (prefix,band)))
            plt.close()

    def curvehist(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            #magbins = np.append([16, 18], np.arange(20, 24.001, 0.5))
            if band == 'g':
                magbins = [20, 24]
            elif band == 'r':
                magbins = [20, 23.5]
            elif band == 'z':
                magbins = [20, 22.5]

            slo,shi = -5,5
            ha = dict(bins=25, range=(slo,shi), histtype='step', normed=True)
            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
            midmag = []
            nn = []
            rgbs = []
            lt,lp = [],[]
            I_empty=True
            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
                I = P[(mag1[P] >= mlo) * (mag1[P] < mhi)]
                if len(I) == 0:
                    continue
                I_empty=False
                ybin = y[I]
                rgb = [0.,0.,0.]
                rgb[0] = float(bini) / (len(magbins)-1)
                rgb[2] = 1. - rgb[0]
                n,b,p = plt.hist(ybin, color=rgb, **ha)
                lt.append('mag %g to %g' % (mlo,mhi))
                lp.append(p[0])
                midmag.append((mlo+mhi)/2.)
                nn.append(n)
                rgbs.append(rgb)

            if I_empty == False:
                bincenters = b[:-1] + (b[1]-b[0])/2.
                lp = []
                for n,rgb,mlo,mhi in zip(nn, rgbs, magbins, magbins[1:]):
                    p = ax[cnt].plot(bincenters, n, '-', color=rgb)
                    lp.append(p[0])
                ax[cnt].plot(bincenters, gaussint[::2], 'k-', alpha=0.5, lw=2)
                ax[cnt].legend(lp, lt)
                ax[cnt].set_xlim(slo,shi)
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%scurvehist_%s.png' % (prefix,band)))
            plt.close()
