from __future__ import print_function
import numpy as np
import fitsio
from astrometry.util.fits import fits_table
from legacypipe.cpimage import CP_DQ_BITS
from legacypipe.survey import tim_get_resamp

def make_coadds(tims, bands, targetwcs,
                mods=None, xy=None, apertures=None, apxy=None,
                ngood=False, detmaps=False, psfsize=False,
                callback=None, callback_args=[],
                plots=False, ps=None,
                lanczos=True, mp=None):
    from astrometry.util.ttime import Time
    t0 = Time()
    
    class Duck(object):
        pass
    C = Duck()

    W = int(targetwcs.get_width())
    H = int(targetwcs.get_height())

    # always, for patching SATUR, etc pixels?
    unweighted=True

    if not xy:
        psfsize = False
    
    C.coimgs = []
    if detmaps:
        C.galdetivs = []
        C.detivs = []
    if mods is not None:
        C.comods = []
        C.coresids = []
        
    if apertures is not None:
        unweighted = True
        C.AP = fits_table()

    if xy:
        ix,iy = xy
        C.T = fits_table()
        C.T.nobs    = np.zeros((len(ix), len(bands)), np.uint8)
        C.T.anymask = np.zeros((len(ix), len(bands)), np.int16)
        C.T.allmask = np.zeros((len(ix), len(bands)), np.int16)
        if psfsize:
            C.T.psfsize = np.zeros((len(ix), len(bands)), np.float32)
        if detmaps:
            C.T.depth    = np.zeros((len(ix), len(bands)), np.float32)
            C.T.galdepth = np.zeros((len(ix), len(bands)), np.float32)

    if lanczos:
        print('Doing Lanczos resampling')

    for tim in tims:
        # surface-brightness correction
        tim.sbscale = (targetwcs.pixel_scale() / tim.subwcs.pixel_scale())**2

    # We create one iterator per band to do the tim resampling.  These all run in
    # parallel when multi-processing.
    imaps = []
    for band in bands:
        args = []
        for itim,tim in enumerate(tims):
            if tim.band != band:
                continue
            if mods is None:
                mo = None
            else:
                mo = mods[itim]
            args.append((itim,tim,mo,lanczos,targetwcs))
        if mp is not None:
            imaps.append(mp.imap_unordered(_resample_one, args))
        else:
            import itertools
            imaps.append(itertools.imap(_resample_one, args))

    # Args for aperture photometry
    apargs = []
            
    tinyw = 1e-30
    for iband,(band,timiter) in enumerate(zip(bands, imaps)):
        print('Computing coadd for band', band)
        
        # coadded weight map (moo)
        cow    = np.zeros((H,W), np.float32)
        # coadded weighted image map
        cowimg = np.zeros((H,W), np.float32)

        kwargs = dict(cowimg=cowimg, cow=cow)

        if detmaps:
            # detection map inverse-variance (depth map)
            detiv = np.zeros((H,W), np.float32)
            C.detivs.append(detiv)
            kwargs.update(detiv=detiv)
            # galaxy detection map inverse-variance (galdepth map)
            galdetiv = np.zeros((H,W), np.float32)
            C.galdetivs.append(galdetiv)
            kwargs.update(galdetiv=galdetiv)

        if mods is not None:
            # model image
            cowmod = np.zeros((H,W), np.float32)
            # chi-squared image
            cochi2 = np.zeros((H,W), np.float32)
            kwargs.update(cowmod=cowmod, cochi2=cochi2)

        if unweighted:
            # unweighted image
            coimg  = np.zeros((H,W), np.float32)
            if mods is not None:
                # unweighted model
                comod  = np.zeros((H,W), np.float32)
            # number of exposures
            con    = np.zeros((H,W), np.uint8)
            # inverse-variance
            coiv   = np.zeros((H,W), np.float32)
            kwargs.update(coimg=coimg, coiv=coiv)

        # Note that we have 'congood' as well as 'nobs':
        # * 'congood' is used for the 'nexp' *image*.
        # * 'nobs' is used for the per-source measurements
        #
        # (you want to know the number of observations within the
        # source footprint, not just the peak pixel which may be
        # saturated, etc.)

        if ngood:
            congood = np.zeros((H,W), np.uint8)
            kwargs.update(congood=congood)

        if xy:
            # These match the type of the "DQ" images.
            # "any" mask
            ormask  = np.zeros((H,W), np.int16)
            # "all" mask
            andmask = np.empty((H,W), np.int16)
            allbits = reduce(np.bitwise_or, CP_DQ_BITS.values())
            andmask[:,:] = allbits
            # number of observations
            nobs  = np.zeros((H,W), np.uint8)
            kwargs.update(ormask=ormask, andmask=andmask, nobs=nobs)

        if psfsize:
            psfsizemap = np.zeros((H,W), np.float32)

        for R in timiter:
            if R is None:
                continue

            itim,Yo,Xo,iv,im,mo,dq = R
            #print('timiter Yo,Xo,im.shape=',Yo,Xo,im.shape)

            tim = tims[itim]

            # invvar-weighted image
            cowimg[Yo,Xo] += iv * im
            cow   [Yo,Xo] += iv

            if unweighted:
                if dq is None:
                    goodpix = 1
                else:
                    # include BLEED, SATUR, INTERP pixels if no other
                    # pixels exists (do this by eliminating all other CP
                    # flags)
                    badbits = 0
                    for bitname in ['badpix', 'cr', 'trans', 'edge', 'edge2']:
                        badbits |= CP_DQ_BITS[bitname]
                    goodpix = ((dq & badbits) == 0)
                    
                coimg[Yo,Xo] += goodpix * im
                con  [Yo,Xo] += goodpix
                coiv [Yo,Xo] += goodpix * 1./(tim.sig1 * tim.sbscale)**2  # ...ish
                
            if xy:
                if dq is not None:
                    ormask [Yo,Xo] |= dq
                    andmask[Yo,Xo] &= dq
                # raw exposure count
                nobs[Yo,Xo] += 1

            if psfsize:
                # psfnorm is in units of 1/pixels.
                # (eg, psfnorm for a gaussian is ~ 1/psf_sigma)
                # Neff is in pixels**2
                neff = 1./tim.psfnorm**2
                # Narcsec is in arcsec**2
                narcsec = neff * tim.wcs.pixel_scale()**2
                psfsizemap[Yo,Xo] += iv * (1. / narcsec)
                
            if detmaps:
                # point-source depth
                detsig1 = tim.sig1 / tim.psfnorm
                detiv[Yo,Xo] += (iv > 0) * (1. / detsig1**2)

                # Galaxy detection map
                gdetsig1 = tim.sig1 / tim.galnorm
                galdetiv[Yo,Xo] += (iv > 0) * (1. / gdetsig1**2)

            if ngood:
                congood[Yo,Xo] += (iv > 0)

            if mods is not None:
                # straight-up
                comod[Yo,Xo] += goodpix * mo
                # invvar-weighted
                cowmod[Yo,Xo] += iv * mo
                # chi-squared
                cochi2[Yo,Xo] += iv * (im - mo)**2
                del mo
                del goodpix

            del Yo,Xo,im,iv
            # END of loop over tims
        # Per-band:
        cowimg /= np.maximum(cow, tinyw)
        C.coimgs.append(cowimg)
        if mods is not None:
            cowmod  /= np.maximum(cow, tinyw)
            C.comods.append(cowmod)
            coresid = cowimg - cowmod
            coresid[cow == 0] = 0.
            C.coresids.append(coresid)

        if unweighted:
            coimg  /= np.maximum(con, 1)
            del con
            cowimg[cow == 0] = coimg[cow == 0]
            if mods is not None:
                cowmod[cow == 0] = comod[cow == 0]

        if xy:
            C.T.nobs [:,iband] = nobs[iy,ix]
            C.T.anymask[:,iband] =  ormask [iy,ix]
            C.T.allmask[:,iband] =  andmask[iy,ix]
            # unless there were no images there...
            C.T.allmask[nobs[iy,ix] == 0, iband] = 0

            if detmaps:
                C.T.depth   [:,iband] =    detiv[iy, ix]
                C.T.galdepth[:,iband] = galdetiv[iy, ix]
        
        if psfsize:
            wt = cow[iy,ix]
            # psfsizemap is in units of iv * (1 / arcsec**2)
            sz = psfsizemap[iy,ix]
            sz /= np.maximum(wt, tinyw)
            sz[wt == 0] = 0.
            # Back to units of linear arcsec.
            sz = 1. / np.sqrt(sz)
            sz[wt == 0] = 0.
            # Correction factor to get back to equivalent of Gaussian sigma
            sz /= (2. * np.sqrt(np.pi))
            # Conversion factor to FWHM (2.35)
            sz *= 2. * np.sqrt(2. * np.log(2.))
            C.T.psfsize[:,iband] = sz
            del psfsizemap

        if apertures is not None:
            # Aperture photometry, using the unweighted "coimg" and
            # "coiv" arrays.
            with np.errstate(divide='ignore'):
                imsigma = 1.0/np.sqrt(coiv)
                imsigma[coiv == 0] = 0

            for irad,rad in enumerate(apertures):
                apargs.append((irad, band, rad, coimg, imsigma, True, apxy))
                if mods is not None:
                    apargs.append((irad, band, rad, coresid, None, False, apxy))

        if callback is not None:
            callback(band, *callback_args, **kwargs)
        # END of loop over bands

    t2 = Time()
    print('coadds: images:', t2-t0)

    if apertures is not None:
        # Aperture phot, in parallel
        if mp is not None:
            apresults = mp.map(_apphot_one, apargs)
        else:
            apresults = map(_apphot_one, apargs)
        del apargs
        apresults = iter(apresults)
        
        for iband,band in enumerate(bands):
            apimg = []
            apimgerr = []
            if mods is not None:
                apres = []
            for irad,rad in enumerate(apertures):
                (airad, aband, isimg, ap_img, ap_err) = apresults.next()
                assert(airad == irad)
                assert(aband == band)
                assert(isimg)
                apimg.append(ap_img)
                apimgerr.append(ap_err)
    
                if mods is not None:
                    (airad, aband, isimg, ap_img, ap_err) = apresults.next()
                    assert(airad == irad)
                    assert(aband == band)
                    assert(not isimg)
                    apres.append(ap_img)
                    assert(ap_err is None)
                
            ap = np.vstack(apimg).T
            ap[np.logical_not(np.isfinite(ap))] = 0.
            C.AP.set('apflux_img_%s' % band, ap)
            ap = 1./(np.vstack(apimgerr).T)**2
            ap[np.logical_not(np.isfinite(ap))] = 0.
            C.AP.set('apflux_img_ivar_%s' % band, ap)
            if mods is not None:
                ap = np.vstack(apres).T
                ap[np.logical_not(np.isfinite(ap))] = 0.
                C.AP.set('apflux_resid_%s' % band, ap)

        t3 = Time()
        print('coadds apphot:', t3-t2)

    return C

def _resample_one((itim,tim,mod,lanczos,targetwcs)):
    from astrometry.util.resample import resample_with_wcs, OverlapError
    if lanczos:
        from astrometry.util.miscutils import patch_image
        patched = tim.getImage().copy()
        okpix = (tim.getInvError() > 0)
        patch_image(patched, okpix)
        del okpix
        imgs = [patched]
        if mod is not None:
            imgs.append(mod)
    else:
        imgs = []

    try:
        Yo,Xo,Yi,Xi,rimgs = resample_with_wcs(
            targetwcs, tim.subwcs, imgs, 3)
    except OverlapError:
        return None
    if len(Yo) == 0:
        return None
    mo = None
    if lanczos:
        im = rimgs[0]
        if mod is not None:
            mo = rimgs[1]
        del patched,imgs,rimgs
    else:
        im = tim.getImage ()[Yi,Xi]
        if mod is not None:
            mo = mods[itim][Yi,Xi]
    iv = tim.getInvvar()[Yi,Xi]
    fscale = tim.sbscale
    print('Applying surface-brightness scaling of %.3f to' % fscale, tim.name)
    im *=  fscale
    iv /= (fscale**2)
    if mod is not None:
        mo *= fscale
    if tim.dq is None:
        dq = None
    else:
        dq = tim.dq[Yi,Xi]
    return itim,Yo,Xo,iv,im,mo,dq

def _apphot_one((irad, band, rad, img, sigma, isimage, apxy)):
    import photutils
    result = [irad, band, isimage]
    aper = photutils.CircularAperture(apxy, rad)
    p = photutils.aperture_photometry(img, aper, error=sigma)
    result.append(p.field('aperture_sum'))
    if sigma is not None:
        result.append(p.field('aperture_sum_err'))
    else:
        result.append(None)
    return result

def write_coadd_images(band,
                       survey, brickname, version_header, tims, targetwcs,
                       cowimg=None, cow=None, cowmod=None, cochi2=None,
                       detiv=None, galdetiv=None, congood=None, **kwargs):

    # copy version_header before modifying...
    hdr = fitsio.FITSHDR()
    for r in version_header.records():
        hdr.add_record(r)
    # Grab these keywords from all input files for this band...
    keys = ['TELESCOP','OBSERVAT','OBS-LAT','OBS-LONG','OBS-ELEV',
            'INSTRUME','FILTER']
    vals = set()
    for tim in tims:
        if tim.band != band:
            continue
        v = []
        for key in keys:
            v.append(tim.primhdr.get(key,''))
        vals.add(tuple(v))
    for i,v in enumerate(vals):
        for ik,key in enumerate(keys):
            if i == 0:
                kk = key
            else:
                kk = key[:7] + '%i'%i
            hdr.add_record(dict(name=kk, value=v[ik]))
    hdr.add_record(dict(name='FILTERX', value=band))

    # DATE-OBS converted to TAI.
    # print('Times:', [tim.time for tim in tims if tim.band == band])
    mjds = [tim.time.toMjd() for tim in tims if tim.band == band]
    minmjd = min(mjds)
    maxmjd = max(mjds)
    #print('MJDs', mjds, 'range', minmjd, maxmjd)
    # back to date string in UTC...
    import astropy.time
    tt = [astropy.time.Time(mjd, format='mjd', scale='tai').utc.isot
          for mjd in [minmjd, maxmjd]]
    hdr.add_record(dict(
        name='DATEOBS1', value=tt[0],
        comment='DATE-OBS for the first image in the stack (UTC)'))
    hdr.add_record(dict(
        name='DATEOBS2', value=tt[1],
        comment='DATE-OBS for the last  image in the stack (UTC)'))

    # Plug the WCS header cards into these images
    targetwcs.add_to_header(hdr)
    hdr.delete('IMAGEW')
    hdr.delete('IMAGEH')
    hdr.add_record(dict(name='EQUINOX', value=2000.))

    imgs = [
        ('image',  'image', cowimg),
        ('invvar', 'wtmap', cow   ),
        ]
    if congood is not None:
        imgs.append(
            ('nexp',   'expmap',   congood),
            )
    if detiv is not None:
        imgs.extend([
                ('depth',    'psfdepth', detiv   ),
                ])
    if galdetiv is not None:
        imgs.extend([
                ('galdepth', 'galdepth', galdetiv),
                ])
    if cowmod is not None:
        imgs.extend([
                ('model',    'model',    cowmod  ),
                ('chi2',     'chi2',     cochi2  ),
                ])
    for name,prodtype,img in imgs:
        from legacypipe.survey import MyFITSHDR
        hdr2 = MyFITSHDR()
        # Make a copy, because each image has different values for
        # these headers...
        #hdr2 = fitsio.FITSHDR()
        for r in hdr.records():
            hdr2.add_record(r)
        hdr2.add_record(dict(name='IMTYPE', value=name,
                             comment='LegacySurvey image type'))
        hdr2.add_record(dict(name='PRODTYPE', value=prodtype,
                             comment='NOAO image type'))
        if name in ['image', 'model']:
            hdr2.add_record(dict(name='MAGZERO', value=22.5,
                                 comment='Magnitude zeropoint'))
            hdr2.add_record(dict(name='BUNIT', value='nanomaggy',
                                 comment='AB mag = 22.5 - 2.5*log10(nanomaggy)'))
        if name in ['invvar', 'depth']:
            hdr2.add_record(dict(name='BUNIT', value='1/nanomaggy^2',
                                 comment='Ivar of ABmag=22.5-2.5*log10(nmgy)'))

        with survey.write_output(name, brick=brickname, band=band) as out:
            fitsio.write(out.fn, img, clobber=True, header=hdr2)
            print('Wrote', out.fn)

# Pretty much only used for plots; the real deal is make_coadds()
def quick_coadds(tims, bands, targetwcs, images=None,
                 get_cow=False, get_n2=False, fill_holes=True):

    W = int(targetwcs.get_width())
    H = int(targetwcs.get_height())

    coimgs = []
    cons = []
    if get_n2:
        cons2 = []
    if get_cow:
        # moo
        cowimgs = []
        wimgs = []

    for ib,band in enumerate(bands):
        coimg = np.zeros((H,W), np.float32)
        coimg2 = np.zeros((H,W), np.float32)
        con   = np.zeros((H,W), np.uint8)
        con2  = np.zeros((H,W), np.uint8)
        if get_cow:
            cowimg = np.zeros((H,W), np.float32)
            wimg  = np.zeros((H,W), np.float32)
        for itim,tim in enumerate(tims):
            if tim.band != band:
                continue
            R = tim_get_resamp(tim, targetwcs)
            if R is None:
                continue
            (Yo,Xo,Yi,Xi) = R
            nn = (tim.getInvError()[Yi,Xi] > 0)
            if images is None:
                coimg [Yo,Xo] += tim.getImage()[Yi,Xi] * nn
                coimg2[Yo,Xo] += tim.getImage()[Yi,Xi]
            else:
                coimg [Yo,Xo] += images[itim][Yi,Xi] * nn
                coimg2[Yo,Xo] += images[itim][Yi,Xi]
            con   [Yo,Xo] += nn
            if get_cow:
                cowimg[Yo,Xo] += tim.getInvvar()[Yi,Xi] * tim.getImage()[Yi,Xi]
                wimg  [Yo,Xo] += tim.getInvvar()[Yi,Xi]
            con2  [Yo,Xo] += 1
        coimg /= np.maximum(con,1)
        if fill_holes:
            coimg[con == 0] = coimg2[con == 0] / np.maximum(1, con2[con == 0])
        if get_cow:
            cowimg /= np.maximum(wimg, 1e-16)
            cowimg[wimg == 0] = coimg[wimg == 0]
            cowimgs.append(cowimg)
            wimgs.append(wimg)
        coimgs.append(coimg)
        cons.append(con)
        if get_n2:
            cons2.append(con2)

    rtn = [coimgs,cons]
    if get_cow:
        rtn.extend([cowimgs, wimgs])
    if get_n2:
        rtn.append(cons2)
    return rtn

