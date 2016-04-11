from __future__ import print_function

import numpy as np
from astrometry.util.ttime import Time
from astrometry.util.resample import resample_with_wcs, OverlapError
from astrometry.util.fits import fits_table
from tractor import Tractor, PointSource, Image, NanoMaggies, Catalog, Patch
from tractor.galaxy import DevGalaxy, ExpGalaxy, FixedCompositeGalaxy, SoftenedFracDev, FracDev, disable_galaxy_cache, enable_galaxy_cache
from legacypipe.common import SimpleGalaxy, LegacyEllipseWithPriors


def one_blob(X):
    '''
    Fits sources contained within a "blob" of pixels.
    '''
    (nblob, iblob, Isrcs, brickwcs, bx0, by0, blobw, blobh, blobmask, timargs,
     srcs, bands, plots, ps, simul_opt, use_ceres, hastycho) = X

    print('Fitting blob number', nblob, 'val', iblob, ':', len(Isrcs),
          'sources, size', blobw, 'x', blobh, len(timargs), 'images')

    if len(timargs) == 0:
        return None

    plots2 = False
    tlast = Time()
    alphas = [0.1, 0.3, 1.0]
    optargs = dict(priors=True, shared_params=False, alphas=alphas)
    bigblob = (blobw * blobh) > 100*100
    trargs = dict()
    
    if use_ceres:
        from tractor.ceres_optimizer import CeresOptimizer
        ceres_optimizer = CeresOptimizer()
        optargs.update(scale_columns=False,
                       scaled=False,
                       dynamic_scale=False)
        trargs.update(optimizer=ceres_optimizer)
        
    # 50 CCDs is over 90th percentile of bricks in DR2.
    many_exposures = len(timargs) >= 50
    #PTF special handling len(timargs) >= 1000

    blobwcs = brickwcs.get_subimage(bx0, by0, blobw, blobh)
    ok,x0,y0 = blobwcs.radec2pixelxy(
        np.array([src.getPosition().ra  for src in srcs]),
        np.array([src.getPosition().dec for src in srcs]))
    started_in_blob = blobmask[np.clip(np.round(y0-1).astype(int), 0,blobh-1),
                               np.clip(np.round(x0-1).astype(int), 0,blobw-1)]

    tims = []
    for (img, inverr, twcs, wcs, pcal, sky, psf, name, sx0, sx1, sy0, sy1,
         band, sig1, modelMinval, imobj) in timargs:

        # Mask out inverr for pixels that are not within the blob.
        subwcs = wcs.get_subimage(int(sx0), int(sy0),
                                  int(sx1-sx0), int(sy1-sy0))
        try:
            Yo,Xo,Yi,Xi,rims = resample_with_wcs(subwcs, blobwcs, [], 2)
        except OverlapError:
            #print('No overlap')
            continue
        if len(Yo) == 0:
            continue
        inverr2 = np.zeros_like(inverr)
        I = np.flatnonzero(blobmask[Yi,Xi])
        inverr2[Yo[I],Xo[I]] = inverr[Yo[I],Xo[I]]
        inverr = inverr2

        # If the subimage (blob) is small enough, instantiate a
        # constant PSF model in the center.
        if sy1-sy0 < 400 and sx1-sx0 < 400:
            subpsf = psf.constantPsfAt((sx0+sx1)/2., (sy0+sy1)/2.)
        else:
            # Otherwise, instantiate a (shifted) spatially-varying
            # PsfEx model.
            subpsf = psf.getShifted(sx0, sy0)

        tim = Image(data=img, inverr=inverr, wcs=twcs,
                    psf=subpsf, photocal=pcal, sky=sky, name=name)
        tim.band = band
        tim.sig1 = sig1
        tim.modelMinval = modelMinval
        tim.subwcs = subwcs
        tim.meta = imobj
        tims.append(tim)

        if plots:
            try:
                Yo,Xo,Yi,Xi,rims = resample_with_wcs(blobwcs, subwcs,[],2)
            except OverlapError:
                continue
            tim.resamp = (Yo, Xo, Yi, Xi)
            if False:
                plt.clf()
                plt.subplot(1,2,1)
                dimshow(img, vmin=-2.*sig1, vmax=5.*sig1)
                plt.subplot(1,2,2)
                dimshow(inverr, vmin=0, vmax=1.1/sig1)
                plt.suptitle('Subimage: ' + name)
                ps.savefig()

    if plots:
        print('Plotting blob image for blob', nblob, 'blob id', iblob)
        coimgs,cons = compute_coadds(tims, bands, blobwcs, fill_holes=False)
        plt.clf()
        dimshow(get_rgb(coimgs, bands))
        ps.savefig()

        plt.plot(x0, y0, 'r.')
        plt.title('initial sources')
        ps.savefig()

        plt.clf()
        ccmap = dict(g='g', r='r', z='m')
        for tim in tims:
            chi = (tim.data * tim.inverr)[tim.inverr > 0]
            plt.hist(chi.ravel(), range=(-5,10), bins=100, histtype='step',
                     color=ccmap[tim.band])
        plt.xlabel('signal/noise per pixel')
        ps.savefig()

    cat = Catalog(*srcs)
    tr = Tractor(tims, cat, **trargs)
    tr.freezeParam('images')

    _fit_fluxes(cat, tims, bands, use_ceres, alphas)
    cat.thawAllRecursive()

    if plots:
        bslc = (slice(by0, by0+blobh), slice(bx0, bx0+blobw))
        plotmods = []
        plotmodnames = []
        plotmods.append(list(tr.getModelImages()))
        plotmodnames.append('Initial models')
        _plot_mods(tims, plotmods, plotmodnames, bands, None, None,
                   bslc, blobw, blobh, ps, chi_plots=False)

    # Optimize individual sources in order of flux
    fluxes = []
    for src in cat:
        # HACK -- here we just *sum* the nanomaggies in each band.  Bogus!
        br = src.getBrightness()
        flux = sum([br.getFlux(band) for band in bands])
        fluxes.append(flux)
    Ibright = np.argsort(-np.array(fluxes))

    if len(cat) > 1:
        # -Remember the original images
        # -Compute initial models for each source (in each tim)
        # -Subtract initial models from images
        # -During fitting, for each source:
        #   -add back in the source's initial model (to each tim)
        #   -fit, with Catalog([src])
        #   -subtract final model (from each tim)
        # -Replace original images

        models = SourceModels()
        # Remember original tim images
        models.save_images(tims)
        # Create & subtract initial models for each tim x each source
        models.create(tims, srcs, subtract=True)

        # For sources, in decreasing order of brightness
        for numi,i in enumerate(Ibright):
            #tsrc = Time()
            print('Fitting source', i, '(%i of %i in blob)' %
                  (numi, len(Ibright)))
            src = cat[i]

            # Add this source's initial model back in.
            models.add(i, tims)

            if bigblob:
                # Create super-local sub-sub-tims around this source

                # Make the subimages the same size as the modelMasks.
                #tbb0 = Time()
                mods = [mod[i] for mod in models.models]
                srctims,modelMasks = _get_subimages(tims, mods, src)
                #print('Creating srctims:', Time()-tbb0)

                if plots and (numi < 3 or numi >= len(Ibright)-3):
                    bx1 = bx0 + blobw
                    by1 = by0 + blobh
                    plt.clf()
                    coimgs,cons = compute_coadds(tims, bands, blobwcs,
                                                 fill_holes=False)
                    dimshow(get_rgb(coimgs, bands), extent=(bx0,bx1,by0,by1))
                    # plt.plot([bx0,bx0,bx1,bx1,bx0],[by0,by1,by1,by0,by0],
                    #          'r-')
                    # for tim in srctims:
                    #     h,w = tim.shape
                    #     tx,ty = [0,0,w,w,0], [0,h,h,0,0]
                    #     rd = [tim.getWcs().pixelToPosition(xi,yi)
                    #           for xi,yi in zip(tx,ty)]
                    #     ra  = [p.ra  for p in rd]
                    #     dec = [p.dec for p in rd]
                    #     ok,x,y = brickwcs.radec2pixelxy(ra, dec)
                    #     plt.plot(x, y, 'g-')
                    #
                    #     ra,dec = tim.subwcs.pixelxy2radec(tx, ty)
                    #     ok,x,y = brickwcs.radec2pixelxy(ra, dec)
                    #     plt.plot(x, y, 'm-')
                    for tim in tims:
                        h,w = tim.shape
                        tx,ty = [0,0,w,w,0], [0,h,h,0,0]
                        rd = [tim.getWcs().pixelToPosition(xi,yi)
                              for xi,yi in zip(tx,ty)]
                        ra  = [p.ra  for p in rd]
                        dec = [p.dec for p in rd]
                        ok,x,y = brickwcs.radec2pixelxy(ra, dec)
                        plt.plot(x, y, 'b-')

                        ra,dec = tim.subwcs.pixelxy2radec(tx, ty)
                        ok,x,y = brickwcs.radec2pixelxy(ra, dec)
                        plt.plot(x, y, 'c-')
                    plt.title('source %i of %i' % (numi, len(Ibright)))
                    ps.savefig()

            else:
                srctims = tims
                modelMasks = models.model_masks(i, src)

            srctractor = Tractor(srctims, [src], **trargs)
            srctractor.freezeParams('images')
            srctractor.setModelMasks(modelMasks)

            if plots and False:
                spmods,spnames = [],[]
                spallmods,spallnames = [],[]
                if numi == 0:
                    spallmods.append(list(tr.getModelImages()))
                    spallnames.append('Initial (all)')
                spmods.append(list(srctractor.getModelImages()))
                spnames.append('Initial')

            # First-round optimization
            print('First-round initial log-prob:', srctractor.getLogProb())
            for step in range(50):
                #R = srctractor.optimize(**optargs)
                #print('Optimize result:', R)
                dlnp,X,alpha = srctractor.optimize(**optargs)
                
                print('dlnp:', dlnp, 'src', src)
                if dlnp < 0.1:
                    break
            print('First-round final log-prob:', srctractor.getLogProb())

            if plots and False:
                spmods.append(list(srctractor.getModelImages()))
                spnames.append('Fit')
                spallmods.append(list(tr.getModelImages()))
                spallnames.append('Fit (all)')

            if plots and False:
                tims_compute_resamp(None, srctractor.getImages(), brickwcs)
                tims_compute_resamp(None, tims, brickwcs)
                plt.figure(1, figsize=(8,6))
                plt.subplots_adjust(left=0.01, right=0.99, top=0.95,
                                    bottom=0.01, hspace=0.1, wspace=0.05)
                #plt.figure(2, figsize=(3,3))
                #plt.subplots_adjust(left=0.005, right=0.995,
                #                    top=0.995,bottom=0.005)
                #_plot_mods(tims, spmods, spnames, bands, None, None, bslc,
                #           blobw, blobh, ps, chi_plots=plots2)
                plt.figure(2, figsize=(3,3.5))
                plt.subplots_adjust(left=0.005, right=0.995,
                                    top=0.88, bottom=0.005)
                plt.suptitle('Blob %i' % iblob)
                tempims = [tim.getImage() for tim in tims]

                _plot_mods(list(srctractor.getImages()), spmods, spnames,
                           bands, None, None, bslc, blobw, blobh, ps,
                           chi_plots=plots2, rgb_plots=True, main_plot=False,
                           rgb_format=('spmods Blob %i, src %i: %%s' %
                                       (iblob, i)))
                _plot_mods(tims, spallmods, spallnames, bands, None, None,
                           bslc, blobw, blobh, ps,
                           chi_plots=plots2, rgb_plots=True, main_plot=False,
                           rgb_format=('spallmods Blob %i, src %i: %%s' %
                                       (iblob, i)))

                models.restore_images(tims)
                _plot_mods(tims, spallmods, spallnames, bands, None, None,
                           bslc, blobw, blobh, ps,
                           chi_plots=plots2, rgb_plots=True, main_plot=False,
                           rgb_format='Blob %i, src %i: %%s' % (iblob, i))
                for tim,im in zip(tims, tempims):
                    tim.data = im

            # Re-remove the final fit model for this source
            models.update_and_subtract(i, src, tims)

            srctractor.setModelMasks(None)
            disable_galaxy_cache()

            #print('Fitting source took', Time()-tsrc)
            #print(src)

        models.restore_images(tims)
        del models

    else:
        # Single source (though this is coded to handle multiple sources)
        # Fit sources one at a time, but don't subtract other models
        cat.freezeAllParams()

        models = SourceModels()
        models.create(tims, [src])
        modelMasks = models.model_masks(0, src)

        tr.setModelMasks(modelMasks)
        enable_galaxy_cache()

        for numi,i in enumerate(Ibright):
            #tsrc = Time()
            #print('Fitting source', i, '(%i of %i in blob)' %
            #  (numi, len(Ibright)))
            cat.freezeAllBut(i)

            for step in range(50):
                dlnp,X,alpha = tr.optimize(**optargs)
                # print('dlnp:', dlnp)
                if dlnp < 0.1:
                    break
            #print('Fitting source took', Time()-tsrc)
            # print(cat[i])

        tr.setModelMasks(None)
        disable_galaxy_cache()

    if plots and False:
        plotmods.append(list(tr.getModelImages()))
        plotmodnames.append('Per Source')

    if len(srcs) > 1 and len(srcs) <= 10:
        #tfit = Time()
        # Optimize all at once?
        cat.thawAllParams()
        #print('Optimizing:', tr)
        # tr.printThawedParams()
        for step in range(20):
            dlnp,X,alpha = tr.optimize(**optargs)
            if dlnp < 0.1:
                break
        #print('Simultaneous fit took:', Time()-tfit)

        if plots and False:
            plotmods.append(list(tr.getModelImages()))
            plotmodnames.append('All Sources')

    if plots and False:
        _plot_mods(tims, plotmods, plotmodnames, bands, None, None,
                   bslc, blobw, blobh, ps)

    # FIXME -- for large blobs, fit strata of sources simultaneously?

    print('Blob finished fitting:', Time()-tlast)
    tlast = Time()

    if plots:
        plt.clf()
        dimshow(get_rgb(coimgs, bands))
        ok,sx,sy = blobwcs.radec2pixelxy(
            np.array([src.getPosition().ra  for src in srcs]),
            np.array([src.getPosition().dec for src in srcs]))
        plt.plot(sx, sy, 'r.')
        plt.title('after source fitting')
        ps.savefig()


    # Next, model selections: point source vs dev/exp vs composite.

    # FIXME -- render initial models and find significant flux overlap
    # (product)??  (Could use the same logic above!)  This would give
    # families of sources to fit simultaneously.  (The
    # not-friends-of-friends version of blobs!)

    # We repeat the "compute & subtract initial models" logic from above.
    # -Remember the original images
    # -Compute initial models for each source (in each tim)
    # -Subtract initial models from images
    # -During fitting, for each source:
    #   -add back in the source's initial model (to each tim)
    #   -fit, with Catalog([src])
    #   -subtract final model (from each tim)
    # -Replace original images

    models = SourceModels()
    # Remember original tim images
    models.save_images(tims)
    # Create initial models for each tim x each source
    # tt = Time()
    models.create(tims, cat, subtract=True)
    # print('Subtracting initial models:', Time()-tt)

    # table of per-source measurements for this blob.
    B = fits_table()
    B.flags = np.zeros(len(Isrcs), np.uint16)
    B.dchisqs = np.zeros((len(Isrcs), 5), np.float32)
    B.sources = srcs
    B.Isrcs = Isrcs
    B.started_in_blob = started_in_blob
    B.all_models        = np.array([{} for i in range(len(Isrcs))])
    B.all_model_fluxivs = np.array([{} for i in range(len(Isrcs))])
    B.all_model_flags   = np.array([{} for i in range(len(Isrcs))])

    del srcs
    del Isrcs
    del started_in_blob

    # Model selection for sources, in decreasing order of brightness
    for numi,i in enumerate(Ibright):

        src = cat[i]
        print('Model selection for source %i of %i in blob' %
              (numi, len(Ibright)))
        #tsel = Time()

        # Add this source's initial model back in.
        models.add(i, tims)

        if bigblob:

            if plots:
                plt.clf()
                for j,tim in enumerate(tims):
                    plt.subplot(len(tims), 2, j+1)
                    dimshow(tim.getImage(), vmin=-2*tim.sig1, vmax=5*tim.sig1)
                    ax = plt.axis()
                    x,y = tim.wcs.positionToPixel(src.getPosition())
                    plt.plot(x, y, 'r.')
                ps.savefig()

            mods = [mod[i] for mod in models.models]
            srctims,modelMasks = _get_subimages(tims, mods, src)

            if plots:
                for j,tim in enumerate(srctims):
                    plt.subplot(len(tims), 2, len(tims)+j+1)
                    dimshow(tim.getImage(), vmin=-2*tim.sig1, vmax=5*tim.sig1)
                    ax = plt.axis()
                    x,y = tim.wcs.positionToPixel(src.getPosition())
                    plt.plot(x, y, 'r.')
                ps.savefig()
            # Create a little local WCS subregion for this source, by
            # resampling non-zero inverrs from the srctims into blobwcs
            insrc = np.zeros((blobh,blobw), bool)
            for tim in srctims:
                try:
                    Yo,Xo,Yi,Xi,nil = resample_with_wcs(blobwcs, tim.subwcs,
                                                        [],2)
                except:
                    continue
                insrc[Yo,Xo] |= (tim.inverr[Yi,Xi] > 0)

            if np.sum(insrc) == 0:
                # No source pixels touching blob... this can happen when a source
                # scatters outside the blob in the fitting stage.
                # Drop the source here.
                B.sources[i] = cat[i] = None
                continue
            yin = np.max(insrc, axis=1)
            xin = np.max(insrc, axis=0)
            yl,yh = np.flatnonzero(yin)[np.array([0,-1])]
            xl,xh = np.flatnonzero(xin)[np.array([0,-1])]
            srcwcs = blobwcs.get_subimage(xl, yl, 1+xh-xl, 1+yh-yl)
            # A mask for which pixels in the 'srcwcs' square are occupied.
            srcpix = insrc[yl:yh+1, xl:xh+1]
            from scipy.ndimage.morphology import binary_erosion
            srcpix2 = binary_erosion(srcpix)
        else:
            modelMasks = models.model_masks(i, src)
            srctims = tims
            srcwcs = blobwcs
            srcpix = None

        srctractor = Tractor(srctims, [src], **trargs)
        srctractor.freezeParams('images')
        srctractor.setModelMasks(modelMasks)
        enable_galaxy_cache()

        if plots:
            tims_compute_resamp(None, srctims, blobwcs)
            plt.clf()
            coimgs,cons = compute_coadds(srctims, bands, blobwcs,
                                         fill_holes=False)
            dimshow(get_rgb(coimgs, bands))#, extent=(bx0,bx1,by0,by1))
            plt.title('Model selection: stage1 data')
            ps.savefig()

            for tim in srctims:
                del tim.resamp
            tims_compute_resamp(None, srctims, srcwcs)
            plt.clf()
            coimgs,cons = compute_coadds(srctims, bands, srcwcs,
                                         fill_holes=False)
            dimshow(get_rgb(coimgs, bands))
            plt.title('Model selection: stage1 data')
            ps.savefig()

            # srch,srcw = srcwcs.shape
            # _plot_mods(srctims, [list(srctractor.getModelImages())],
            #            ['Model selection init'], bands, None, None,
            #            None, srch,srcw, ps, chi_plots=False)

        # use log likelihood rather than log prior because we use
        # priors on, eg, the ellipticity to avoid crazy fits.  Here we
        # sort of want to just evaluate the fit quality regardless of
        # priors on parameters...?

        srccat = srctractor.getCatalog()

        # Compute the log-likehood without a source here.
        srccat[0] = None
        chisqs_none = _per_band_chisqs(srctractor, bands)

        nparams = dict(ptsrc=2, simple=2, exp=5, dev=5, comp=9)
        # This is our "upgrade" threshold: how much better a galaxy
        # fit has to be versus ptsrc, and comp versus galaxy.
        galaxy_margin = 3.**2 + (nparams['exp'] - nparams['ptsrc'])

        # *chisqs* is actually chi-squared improvement vs no source;
        # larger is a better fit.
        chisqs = dict(none=0)

        oldmodel, ptsrc, simple, dev, exp, comp = _initialize_models(src)

        trymodels = [('ptsrc', ptsrc), ('simple', simple)]

        if oldmodel == 'ptsrc':
            # Try galaxy models if simple > ptsrc, or if bright.
            # The 'gals' model is just a marker
            trymodels.extend([('gals', None)])
        else:
            trymodels.extend([('dev', dev), ('exp', exp), ('comp', comp)])

        # If lots of exposures, cut to a subset that reach the DECaLS
        # depth goals and use those in an initial round?
        if many_exposures:
            timsubset = set()

            for band in bands:
                # Order to try them: first, DECaLS data (our propid),
                # then in point-source depth (* npix?) order.
                otims = []
                value = []
                for tim in srctims:
                    if tim.band != band:
                        continue
                    otims.append(tim)

                    detsig1 = tim.sig1 / tim.meta.galnorm
                    tim.detiv1 = 1./detsig1**2
                    h,w = tim.shape
                    propid = tim.meta.propid
                    value.append((propid == DECALS_PROPID) * 1e12 +
                                 tim.detiv1 * h*w)

                t1,t2,t3 = dict(g=(24.0, 23.7, 23.4),
                                r=(23.4, 23.1, 22.8),
                                z=(22.5, 22.2, 21.9))[band]
                Nsigma = 5.
                sig = NanoMaggies.magToNanomaggies(t1) / Nsigma
                target1 = 1./sig**2
                sig = NanoMaggies.magToNanomaggies(t2) / Nsigma
                target2 = 1./sig**2
                sig = NanoMaggies.magToNanomaggies(t3) / Nsigma
                target3 = 1./sig**2

                detiv = np.zeros(srcwcs.shape, np.float32)
                I = np.argsort(-np.array(value))
                for cnt in I:
                    tim = otims[cnt]
                    try:
                        Yo,Xo,Yi,Xi,nil = resample_with_wcs(
                            srcwcs, tim.subwcs, [], 2)
                    except:
                        continue
                    if len(Yo) == 0:
                        continue
                    detiv[Yo,Xo] += tim.detiv1
                    timsubset.add(tim.name)

                    print('Tim:', tim.name, 'exptime', tim.meta.exptime,
                          'FWHM', tim.meta.fwhm)
                    print('Tim: detiv', tim.detiv1, 'depth mag',
                          NanoMaggies.nanomaggiesToMag(np.sqrt(1./tim.detiv1)
                                                       * Nsigma))
                    print('overlapping pixels:', len(Yo))

                    # Hit DECaLS depth targets?
                    pctiles = [100-90, 100-95, 100-98]
                    if srcpix is None:
                        p1,p2,p3 = np.percentile(detiv, pctiles)
                    else:
                        p1,p2,p3 = np.percentile(detiv[srcpix2], pctiles)

                    m1 = NanoMaggies.nanomaggiesToMag(np.sqrt(1./p1) * Nsigma)
                    m2 = NanoMaggies.nanomaggiesToMag(np.sqrt(1./p2) * Nsigma)
                    m3 = NanoMaggies.nanomaggiesToMag(np.sqrt(1./p3) * Nsigma)
                    print('Added image', tim.name, 'and got depths (mag)',
                          '%.2f, %.2f, %.2f' % (m1,m2,m3),
                          'vs target mags', t1, t2, t3)
                    print('  detivs', p1, p2, p3, 'vs targets',
                          target1, target2, target3)

                    if p1 >= target1 and p2 >= target2 and p3 >= target3:
                        # Got enough depth, thank you!
                        print('Reached target depth!')
                        break

            if plots:
                dtims = [tim for tim in srctims if tim.name in timsubset]
                plt.clf()
                coimgs,cons = compute_coadds(dtims, bands, srcwcs,
                                             fill_holes=False)
                dimshow(get_rgb(coimgs, bands))
                plt.title('To-depth data')
                ps.savefig()

        allflags = {}
        for name,newsrc in trymodels:

            if name == 'gals':
                # If 'simple' was better than 'ptsrc', or the source is
                # bright, try the galaxy models.
                if ((chisqs['simple'] > chisqs['ptsrc']) or
                    (chisqs['ptsrc'] > 400)):

                    if hastycho:
                        print('Not computing galaxy models: Tycho-2 star in '+
                              'blob')
                        continue

                    trymodels.extend([
                        ('dev', dev), ('exp', exp), ('comp', comp)])
                continue

            if name == 'comp' and newsrc is None:
                # Compute the comp model if exp or dev would be accepted
                if (max(chisqs['dev'], chisqs['exp']) <
                    (chisqs['ptsrc'] + galaxy_margin)):
                    #print('dev/exp not much better than ptsrc; not computing
                    # comp model.')
                    continue
                newsrc = comp = FixedCompositeGalaxy(
                    src.getPosition(), src.getBrightness(),
                    SoftenedFracDev(0.5), exp.getShape(),
                    dev.getShape()).copy()
            #print('New source:', newsrc)
            srccat[0] = newsrc

            # Use the same modelMask shapes as the original source ('src').
            # Need to create newsrc->mask mappings though:
            mm = []
            for mim in modelMasks:
                d = dict()
                mm.append(d)
                try:
                    d[newsrc] = mim[src]
                except KeyError:
                    pass
            srctractor.setModelMasks(mm)
            enable_galaxy_cache()

            # Save these modelMasks for later...
            newsrc_mm = mm

            #lnp = srctractor.getLogProb()
            #print('Initial log-prob:', lnp)
            #print('vs original src: ', lnp - lnp0)
            if plots and False:
                # Grid of derivatives.
                _plot_derivs(tims, newsrc, srctractor, ps)
            if plots:
                mods = list(srctractor.getModelImages())
                plt.clf()
                coimgs,cons = compute_coadds(srctims, bands, srcwcs,
                                             images=mods, fill_holes=False)
                dimshow(get_rgb(coimgs, bands))
                plt.title('Initial: ' + name)
                ps.savefig()

            if many_exposures:
                # Run a quick round of optimization with our to-depth subset
                dtims = []
                dmm = []
                for tim,mim in zip(srctims, modelMasks):
                    if not tim.name in timsubset:
                        continue
                    dtims.append(tim)
                    d = dict()
                    dmm.append(d)
                    try:
                        d[newsrc] = mim[src]
                    except KeyError:
                        pass

                dtractor = Tractor(dtims, [newsrc], **trargs)
                dtractor.freezeParams('images')
                dtractor.setModelMasks(dmm)
                enable_galaxy_cache()

                for step in range(50):
                    #print('optimizing:', newsrc)
                    dlnp,X,alpha = dtractor.optimize(**optargs)
                    #print('  dlnp:', dlnp, 'new src', newsrc)
                    if dlnp < 0.1:
                        break

                # print('Mod', name, 'round0 opt', Time()-t0)
                # print('New source (after to-depth round optimization):',
                #   newsrc)

                if plots:
                    plt.clf()
                    modimgs = list(dtractor.getModelImages())
                    comods,nil = compute_coadds(dtims, bands, srcwcs,
                                                images=modimgs)
                    dimshow(get_rgb(comods, bands))
                    plt.title('To-depth opt: ' + name)
                    ps.savefig()

            # First-round optimization (during model selection)
            thisflags = 0
            for step in range(50):
                #print('optimizing:', newsrc)
                dlnp,X,alpha = srctractor.optimize(**optargs)
                #print('  dlnp:', dlnp, 'new src', newsrc)
                if dlnp < 0.1:
                    break
            else:
                thisflags |= FLAG_STEPS_A

            # print('Mod', name, 'round1 opt', Time()-t0)
            #print('New source (after first round optimization):', newsrc)

            if plots:
                # _plot_mods(srctims, [list(srctractor.getModelImages())],
                #            ['Model selection: ' + name], bands, None, None,
                #            None, srch,srcw, ps, chi_plots=False)
                plt.clf()
                modimgs = list(srctractor.getModelImages())
                comods,nil = compute_coadds(srctims, bands, srcwcs,
                                            images=modimgs)
                dimshow(get_rgb(comods, bands))
                plt.title('First-round opt: ' + name)
                ps.savefig()

            srctractor.setModelMasks(None)
            disable_galaxy_cache()

            # Recompute modelMasks in the original tims

            ## FIXME -- avoid huge patches?  Clip to significant
            ## pixels in model?  Via minval, I guess; appealing to use
            ## same cut as determining the blobs, but that's in brick
            ## coadd space.

            tim = tims[0]
            from tractor.galaxy import ProfileGalaxy
            if isinstance(newsrc, ProfileGalaxy):
                px,py = tim.wcs.positionToPixel(newsrc.getPosition())
                h = newsrc._getUnitFluxPatchSize(tim, px, py, tim.modelMinval)
                MAXHALF = 128
                if h > MAXHALF:
                    print('halfsize', h,'for',newsrc,'-> setting to',MAXHALF)
                    newsrc.halfsize = MAXHALF

            if hastycho:
                modtims = None
            elif bigblob:
                mods = []
                for tim in tims:
                    mod = newsrc.getModelPatch(tim)
                    if mod is not None:
                        h,w = tim.shape
                        mod.clipTo(w,h)
                        if mod.patch is None:
                            mod = None
                    mods.append(mod)
                modtims,mm = _get_subimages(tims, mods, newsrc)
            else:
                mm = []
                modtims = []
                for tim in tims:
                    d = dict()
                    mod = newsrc.getModelPatch(tim)
                    if mod is None:
                        continue
                    #print('After first-round fit: model is', mod.shape)
                    mod = _clip_model_to_blob(mod, tim.shape,tim.getInvError())
                    if mod is None:
                        continue
                    d[newsrc] = Patch(mod.x0, mod.y0, mod.patch != 0)
                    modtims.append(tim)
                    mm.append(d)

            if modtims is not None:
                modtractor = Tractor(modtims, [newsrc], **trargs)
                modtractor.freezeParams('images')
                modtractor.setModelMasks(mm)
                enable_galaxy_cache()

                #t0 = Time()

                # Run another round of opt.
                for step in range(50):
                    dlnp,X,alpha = modtractor.optimize(**optargs)
                    #print('  dlnp:', dlnp, 'new src', newsrc)
                    if dlnp < 0.1:
                        break
                else:
                    thisflags |= FLAG_STEPS_B

                # print('Mod', name, 'round2 opt', Time()-t0)

                if plots:
                    plt.clf()
                    modimgs = list(modtractor.getModelImages())
                    tims_compute_resamp(None, modtims, srcwcs)
                    comods,nil = compute_coadds(modtims, bands, srcwcs,
                                                images=modimgs)
                    dimshow(get_rgb(comods, bands))
                    plt.title('Second-round opt: ' + name)
                    ps.savefig()
            else:
                # Tycho-2 star; set modtractor = srctractor for the ivars
                srctractor.setModelMasks(newsrc_mm)
                modtractor = srctractor

            # Compute FLUX inverse-variances for each source.
            # This uses the second-round modelMasks.
            newsrc.freezeAllBut('brightness')
            allderivs = modtractor.getDerivs()
            ivs = np.zeros(len(bands), np.float32)
            B.all_model_fluxivs[i][name] = ivs
            for iparam,derivs in enumerate(allderivs):
                chisq = 0
                for deriv,tim in derivs:
                    h,w = tim.shape
                    deriv.clipTo(w,h)
                    ie = tim.getInvError()
                    slc = deriv.getSlice(ie)
                    chi = deriv.patch * ie[slc]
                    chisq += (chi**2).sum()
                ivs[iparam] = chisq
            newsrc.thawAllParams()

            # Use the original 'srctractor' here so that the different
            # models are evaluated on the same pixels.
            # ---> AND with the same modelMasks as the original source...
            # FIXME -- it is not clear that this is what we want!!!
            srctractor.setModelMasks(newsrc_mm)
            ch = _per_band_chisqs(srctractor, bands)
            chisqs[name] = _chisq_improvement(newsrc, ch, chisqs_none)
            B.all_models[i][name] = newsrc.copy()
            B.all_model_flags[i][name] = thisflags

        # if plots:
        #    _plot_mods(tims, plotmods, plotmodnames, bands, None, None,
        #               bslc, blobw, blobh, ps)

        if plots:
            from collections import OrderedDict
            plt.clf()
            rows,cols = 2, 6
            mods = OrderedDict([('none',None), ('ptsrc',ptsrc),
                                ('simple',simple),
                                ('dev',dev), ('exp',exp), ('comp',comp)])
            for imod,modname in enumerate(mods.keys()):

                if mod != 'none' and not modname in chisqs:
                    continue

                srccat[0] = mods[modname]

                print('Plotting model for blob', iblob, 'source', i,
                      ':', modname)
                print(srccat[0])

                srctractor.setModelMasks(None)

                plt.subplot(rows, cols, imod+1)

                if modname != 'none':
                    modimgs = list(srctractor.getModelImages())
                    comods,nil = compute_coadds(srctims, bands, srcwcs,
                                                images=modimgs)
                    dimshow(get_rgb(comods, bands), ticks=False)
                    plt.title(modname)
                    chis = [((tim.getImage() - mod) * tim.getInvError())**2
                            for tim,mod in zip(srctims, modimgs)]
                    res = [(tim.getImage() - mod) for tim,mod in
                           zip(srctims, modimgs)]
                else:
                    coimgs, cons = compute_coadds(srctims, bands, srcwcs)
                    dimshow(get_rgb(coimgs, bands))
                    ax = plt.axis()
                    ok,x,y = blobwcs.radec2pixelxy(src.getPosition().ra,
                                                     src.getPosition().dec)
                    plt.plot(x-1, y-1, 'r+')
                    plt.axis(ax)
                    plt.title('Image')
                    chis = [((tim.getImage()) * tim.getInvError())**2
                              for tim in srctims]
                    res = [tim.getImage() for tim in srctims]

                if False:
                    cochisqs,nil = compute_coadds(tims, bands, blobwcs,
                                                  images=chis)
                    cochisq = reduce(np.add, cochisqs)
                    plt.subplot(rows, cols, imod+1+cols)
                    dimshow(cochisq, vmin=0, vmax=25)

                else:
                    # residuals
                    coresids,nil = compute_coadds(srctims, bands, srcwcs,
                                                  images=res)
                    plt.subplot(rows, cols, imod+1+cols)
                    dimshow(get_rgb(coresids, bands, **rgbkwargs_resid),
                            ticks=False)
                plt.title('chisq %.0f' % chisqs[modname], fontsize=8)
            plt.suptitle('Blob %i, source %i: was: \n%s' %
                         (iblob, i, str(src)), fontsize=10)
            ps.savefig()

        # This determines the order of the elements in the DCHISQ
        # column of the catalog.
        modnames = ['ptsrc', 'simple', 'dev', 'exp', 'comp']

        keepmod = _select_model(chisqs, nparams, galaxy_margin)

        keepsrc = dict(none=None, ptsrc=ptsrc, simple=simple,
                       dev=dev, exp=exp, comp=comp)[keepmod]

        B.dchisqs[i, :] = np.array([chisqs.get(k,0) for k in modnames])
        B.flags[i] = allflags.get(keepmod, 0)
        B.sources[i] = keepsrc
        cat[i] = keepsrc

        # Re-remove the final fit model for this source.
        models.update_and_subtract(i, keepsrc, tims)

        #print('Keeping model:', keepmod)
        #print('Keeping source:', keepsrc)
        #print(Time() - tsel)

    models.restore_images(tims)
    del models

    print('Blob finished model selection:', Time()-tlast)
    tlast = Time()

    if plots:
        plotmods, plotmodnames = [],[]
        plotmods.append(list(tr.getModelImages()))
        plotmodnames.append('All model selection')
        _plot_mods(tims, plotmods, plotmodnames, bands, None, None, bslc,
                   blobw, blobh, ps)

    I = np.array([i for i,s in enumerate(cat) if s is not None])
    B.cut(I)

    cat = Catalog(*B.sources)
    tr.catalog = cat

    # Do another quick round of flux-only fitting?
    # This does horribly -- fluffy galaxies go out of control because
    # they're only constrained by pixels within this blob.
    #_fit_fluxes(cat, tims, bands, use_ceres, alphas)

    # print('After cutting sources:')
    # for src,dchi in zip(B.sources, B.dchisqs):
    #     print('  source', src, 'max dchisq', max(dchi), 'dchisqs', dchi)

    ### Simultaneous re-opt?
    if simul_opt and len(cat) > 1 and len(cat) <= 10:
        #tfit = Time()
        cat.thawAllParams()
        #print('Optimizing:', tr)
        #tr.printThawedParams()
        flags |= FLAG_TRIED_C
        max_cpu = 300.
        cpu0 = time.clock()
        for step in range(50):
            dlnp,X,alpha = tr.optimize(**optargs)
            cpu = time.clock()
            if cpu-cpu0 > max_cpu:
                print('Warning: Exceeded maximum CPU time for source')
                flags |= FLAG_CPU_C
                break
            if dlnp < 0.1:
                break
        #print('Simultaneous fit took:', Time()-tfit)

    # Compute variances on all parameters for the kept model
    B.srcinvvars = [[] for i in range(len(B))]
    cat.thawAllRecursive()
    cat.freezeAllParams()
    for isub in range(len(B.sources)):
        cat.thawParam(isub)
        src = cat[isub]
        if src is None:
            cat.freezeParam(isub)
            continue
        # Convert to "vanilla" ellipse parameterization
        if isinstance(src, (DevGalaxy, ExpGalaxy)):
            src.shape = src.shape.toEllipseE()
        elif isinstance(src, FixedCompositeGalaxy):
            src.shapeExp = src.shapeExp.toEllipseE()
            src.shapeDev = src.shapeDev.toEllipseE()
            src.fracDev = FracDev(src.fracDev.clipped())

        #print('Computing variances for source:', src)
        #tr.printThawedParams()
        #print('halfsize:', getattr(src, 'halfsize', None))

        allderivs = tr.getDerivs()
        for iparam,derivs in enumerate(allderivs):
            chisq = 0
            for deriv,tim in derivs:
                h,w = tim.shape
                deriv.clipTo(w,h)
                ie = tim.getInvError()
                slc = deriv.getSlice(ie)
                chi = deriv.patch * ie[slc]
                chisq += (chi**2).sum()
            B.srcinvvars[isub].append(chisq)
        assert(len(B.srcinvvars[isub]) == cat[isub].numberOfParams())
        cat.freezeParam(isub)
    #print('Blob variances:', Time()-tlast)
    #tlast = Time()

    # Check for sources with zero inverse-variance -- I think these
    # can be generated during the "Simultaneous re-opt" stage above --
    # sources can get scattered outside the blob.
    # Arbitrarily look at the first element (RA)
    I = np.flatnonzero(np.array([iv[0] for iv in B.srcinvvars]))
    if len(I) < len(B):
        print('Keeping', len(I), 'of', len(B), 'sources with non-zero ivar')
        B.cut(I)
        cat = Catalog(*B.sources)
        tr.catalog = cat

    M = _compute_source_metrics(B.sources, tims, bands, tr)
    for k,v in M.items():
        B.set(k, v)

    ok,x1,y1 = blobwcs.radec2pixelxy(
        np.array([src.getPosition().ra  for src in B.sources]),
        np.array([src.getPosition().dec for src in B.sources]))
    B.finished_in_blob = blobmask[
        np.clip(np.round(y1-1).astype(int), 0, blobh-1),
        np.clip(np.round(x1-1).astype(int), 0, blobw-1)]
    assert(len(B.finished_in_blob) == len(B))
    assert(len(B.finished_in_blob) == len(B.started_in_blob))

    B.hastycho = np.zeros(len(B), bool)
    if hastycho:
        B.hastycho[:] = True

    print('Blob number', nblob, 'iblob', iblob, 'finished:', Time()-tlast)

    B.iblob = iblob

    return B

def _compute_source_metrics(srcs, tims, bands, tr):
    # rchi2 quality-of-fit metric
    rchi2_num    = np.zeros((len(srcs),len(bands)), np.float32)
    rchi2_den    = np.zeros((len(srcs),len(bands)), np.float32)

    # fracflux degree-of-blending metric
    fracflux_num = np.zeros((len(srcs),len(bands)), np.float32)
    fracflux_den = np.zeros((len(srcs),len(bands)), np.float32)

    # fracin flux-inside-blob metric
    fracin_num = np.zeros((len(srcs),len(bands)), np.float32)
    fracin_den = np.zeros((len(srcs),len(bands)), np.float32)

    # fracmasked: fraction of masked pixels metric
    fracmasked_num = np.zeros((len(srcs),len(bands)), np.float32)
    fracmasked_den = np.zeros((len(srcs),len(bands)), np.float32)

    for iband,band in enumerate(bands):
        for tim in tims:
            if tim.band != band:
                continue
            mod = np.zeros(tim.getModelShape(), tr.modtype)
            srcmods = [None for src in srcs]
            counts = np.zeros(len(srcs))
            pcal = tim.getPhotoCal()

            # For each source, compute its model and record its flux
            # in this image.  Also compute the full model *mod*.
            for isrc,src in enumerate(srcs):
                patch = tr.getModelPatch(tim, src, minsb=tim.modelMinval)
                if patch is None or patch.patch is None:
                    continue
                counts[isrc] = np.sum([np.abs(pcal.brightnessToCounts(b))
                                              for b in src.getBrightnesses()])
                if counts[isrc] == 0:
                    continue
                H,W = mod.shape
                patch.clipTo(W,H)
                srcmods[isrc] = patch
                patch.addTo(mod)

            # Now compute metrics for each source
            for isrc,patch in enumerate(srcmods):
                if patch is None:
                    continue
                if patch.patch is None:
                    continue
                if counts[isrc] == 0:
                    continue
                if np.sum(patch.patch**2) == 0:
                    continue
                slc = patch.getSlice(mod)
                patch = patch.patch

                # print('fracflux: band', band, 'isrc', isrc, 'tim', tim.name)
                # print('src:', srcs[isrc])
                # print('patch sum', np.sum(patch),'abs',np.sum(np.abs(patch)))
                # print('counts:', counts[isrc])
                # print('mod slice sum', np.sum(mod[slc]))
                # print('mod[slc] - patch:', np.sum(mod[slc] - patch))

                # (mod - patch) is flux from others
                # (mod - patch) / counts is normalized flux from others
                # We take that and weight it by this source's profile;
                #  patch / counts is unit profile
                # But this takes the dot product between the profiles,
                # so we have to normalize appropriately, ie by
                # (patch**2)/counts**2; counts**2 drops out of the
                # denom.  If you have an identical source with twice the flux,
                # this results in fracflux being 2.0

                # fraction of this source's flux that is inside this patch.
                # This can be < 1 when the source is near an edge, or if the
                # source is a huge diffuse galaxy in a small patch.
                fin = np.abs(np.sum(patch) / counts[isrc])

                # print('fin:', fin)
                # print('fracflux_num: fin *',
                #      np.sum((mod[slc] - patch) * np.abs(patch)) /
                #      np.sum(patch**2))

                fracflux_num[isrc,iband] += (fin *
                    np.sum((mod[slc] - patch) * np.abs(patch)) /
                    np.sum(patch**2))
                fracflux_den[isrc,iband] += fin
                
                fracmasked_num[isrc,iband] += (
                    np.sum((tim.getInvError()[slc] == 0) * np.abs(patch)) /
                    np.abs(counts[isrc]))
                    
                fracmasked_den[isrc,iband] += fin

                fracin_num[isrc,iband] += np.abs(np.sum(patch))
                fracin_den[isrc,iband] += np.abs(counts[isrc])

            tim.getSky().addTo(mod)
            chisq = ((tim.getImage() - mod) * tim.getInvError())**2

            for isrc,patch in enumerate(srcmods):
                if patch is None or patch.patch is None:
                    continue
                if counts[isrc] == 0:
                    continue
                slc = patch.getSlice(mod)
                # We compute numerator and denom separately to handle
                # edge objects, where sum(patch.patch) < counts.
                # Also, to normalize by the number of images.  (Being
                # on the edge of an image is like being in half an
                # image.)
                rchi2_num[isrc,iband] += (np.sum(chisq[slc] * patch.patch) / 
                                          counts[isrc])
                # If the source is not near an image edge,
                # sum(patch.patch) == counts[isrc].
                rchi2_den[isrc,iband] += np.sum(patch.patch) / counts[isrc]

    #print('Fracflux_num:', fracflux_num)
    #print('Fracflux_den:', fracflux_den)
                
    fracflux   = fracflux_num   / fracflux_den
    rchi2      = rchi2_num      / rchi2_den
    fracmasked = fracmasked_num / fracmasked_den

    #print('Fracflux:', B.fracflux)
    
    # fracin_{num,den} are in flux * nimages units
    tinyflux = 1e-9
    fracin     = fracin_num     / np.maximum(tinyflux, fracin_den)

    return dict(fracin=fracin, fracflux=fracflux, rchi2=rchi2,
                fracmasked=fracmasked)




def _initialize_models(src):
    if isinstance(src, PointSource):
        ptsrc = src.copy()
        simple = SimpleGalaxy(src.getPosition(), src.getBrightness()).copy()
        # logr, ee1, ee2
        shape = LegacyEllipseWithPriors(-1., 0., 0.)
        dev = DevGalaxy(src.getPosition(), src.getBrightness(), shape).copy()
        exp = ExpGalaxy(src.getPosition(), src.getBrightness(), shape).copy()
        comp = None
        oldmodel = 'ptsrc'

    elif isinstance(src, DevGalaxy):
        ptsrc = PointSource(src.getPosition(), src.getBrightness()).copy()
        simple = SimpleGalaxy(src.getPosition(), src.getBrightness()).copy()
        dev = src.copy()
        exp = ExpGalaxy(src.getPosition(), src.getBrightness(),
                        src.getShape()).copy()
        comp = None
        oldmodel = 'dev'

    elif isinstance(src, ExpGalaxy):
        ptsrc = PointSource(src.getPosition(), src.getBrightness()).copy()
        simple = SimpleGalaxy(src.getPosition(), src.getBrightness()).copy()
        dev = DevGalaxy(src.getPosition(), src.getBrightness(),
                        src.getShape()).copy()
        exp = src.copy()
        comp = None
        oldmodel = 'exp'

    elif isinstance(src, FixedCompositeGalaxy):
        ptsrc = PointSource(src.getPosition(), src.getBrightness()).copy()
        simple = SimpleGalaxy(src.getPosition(), src.getBrightness()).copy()
        frac = src.fracDev.clipped()
        if frac > 0:
            shape = src.shapeDev
        else:
            shape = src.shapeExp
        dev = DevGalaxy(src.getPosition(), src.getBrightness(), shape).copy()
        if frac < 1:
            shape = src.shapeExp
        else:
            shape = src.shapeDev
        exp = ExpGalaxy(src.getPosition(), src.getBrightness(), shape).copy()
        comp = src.copy()
        oldmodel = 'comp'

    return oldmodel, ptsrc, simple, dev, exp, comp

def _get_subimages(tims, mods, src):
    srctims = []
    modelMasks = []
    #print('Big blob: trimming:')
    for tim,mod in zip(tims, mods):
        if mod is None:
            continue
        mh,mw = mod.shape
        if mh == 0 or mw == 0:
            continue
        # for modelMasks
        d = { src: Patch(0, 0, mod.patch != 0) }
        modelMasks.append(d)

        x0,y0 = mod.x0 , mod.y0
        x1,y1 = x0 + mw, y0 + mh
        slc = slice(y0,y1), slice(x0, x1)

        subimg = tim.getImage()[slc]
        if subimg.shape != (mh,mw):
            print('Subimage shape:', subimg.shape, 'image shape',
                  tim.getImage().shape, 'slice y', y0,y1, 'x', x0,x1,
                  'mod shape', mh,mw)
        # print('  srctim: x0,y0', x0,y0, 'shape', (y1-y0,x1-x0))
        subpsf = tim.psf.constantPsfAt((x0+x1)/2., (y0+y1)/2.)
        srctim = Image(data=subimg,
                       inverr=tim.getInvError()[slc],
                       wcs=tim.wcs.shifted(x0, y0),
                       psf=subpsf,
                       photocal=tim.getPhotoCal(),
                       sky=tim.sky.shifted(x0, y0),
                       name=tim.name)
        sh,sw = srctim.shape
        srctim.subwcs = tim.subwcs.get_subimage(x0, y0, sw, sh)
        srctim.band = tim.band
        srctim.sig1 = tim.sig1
        srctim.modelMinval = tim.modelMinval
        srctim.x0 = x0
        srctim.y0 = y0
        srctim.meta = tim.meta
        srctims.append(srctim)
        #print('  ', tim.shape, 'to', srctim.shape)
    return srctims, modelMasks

class SourceModels(object):
    '''
    This class maintains a list of the model patches for a set of sources
    in a set of images.
    '''
    def save_images(self, tims):
        self.orig_images = [tim.getImage() for tim in tims]
        for tim,img in zip(tims, self.orig_images):
            tim.data = img.copy()

    def restore_images(self, tims):
        for tim,img in zip(tims, self.orig_images):
            tim.data = img

    def create(self, tims, srcs, subtract=False):
        '''
        Note that this modifies the *tims* if subtract=True.
        '''
        self.models = []
        for tim in tims:
            mods = []
            sh = tim.shape
            ie = tim.getInvError()
            for src in srcs:
                mod = src.getModelPatch(tim)
                if mod is not None and mod.patch is not None:
                    if not np.all(np.isfinite(mod.patch)):
                        print('Non-finite mod patch')
                        print('source:', src)
                        print('tim:', tim)
                        print('PSF:', tim.getPsf())
                    assert(np.all(np.isfinite(mod.patch)))
                    mod = _clip_model_to_blob(mod, sh, ie)
                    if subtract and mod is not None:
                        mod.addTo(tim.getImage(), scale=-1)
                mods.append(mod)
            self.models.append(mods)

    def add(self, i, tims):
        '''
        Adds the models for source *i* back into the tims.
        '''
        for tim,mods in zip(tims, self.models):
            mod = mods[i]
            if mod is not None:
                mod.addTo(tim.getImage())

    def update_and_subtract(self, i, src, tims):
        for tim,mods in zip(tims, self.models):
            #mod = srctractor.getModelPatch(tim, src)
            if src is None:
                mod = None
            else:
                mod = src.getModelPatch(tim)
            if mod is not None:
                mod.addTo(tim.getImage(), scale=-1)
            mods[i] = mod

    def model_masks(self, i, src):
        modelMasks = []
        for mods in self.models:
            d = dict()
            modelMasks.append(d)
            mod = mods[i]
            if mod is not None:
                d[src] = Patch(mod.x0, mod.y0, mod.patch != 0)
        return modelMasks



def _clip_model_to_blob(mod, sh, ie):
    '''
    mod: Patch
    sh: tim shape
    ie: tim invError
    Returns: new Patch
    '''
    mslc,islc = mod.getSlices(sh)
    sy,sx = mslc
    patch = mod.patch[mslc] * (ie[islc]>0)
    if patch.shape == (0,0):
        return None
    mod = Patch(mod.x0 + sx.start, mod.y0 + sy.start, patch)

    # Check
    mh,mw = mod.shape
    assert(mod.x0 >= 0)
    assert(mod.y0 >= 0)
    ph,pw = sh
    assert(mod.x0 + mw <= pw)
    assert(mod.y0 + mh <= ph)

    return mod

FLAG_CPU_A   = 1
FLAG_STEPS_A = 2
FLAG_CPU_B   = 4
FLAG_STEPS_B = 8
FLAG_TRIED_C = 0x10
FLAG_CPU_C   = 0x20

def _select_model(chisqs, nparams, galaxy_margin):
    '''
    Returns keepmod
    '''
    keepmod = 'none'

    # This is our "detection threshold": 5-sigma in
    # *parameter-penalized* units; ie, ~5.2-sigma for point sources
    cut = 5.**2
    # Take the best of all models computed
    diff = max([chisqs[name] - nparams[name] for name in chisqs.keys()
                if name != 'none'])

    if diff < cut:
        return keepmod

    # We're going to keep this source!
    if chisqs['ptsrc'] > chisqs['simple']:
        #print('Keeping source; PTSRC is better than SIMPLE')
        keepmod = 'ptsrc'
    else:
        #print('Keeping source; SIMPLE is better than PTSRC')
        keepmod = 'simple'

    if not 'exp' in chisqs:
        return keepmod

    # This is our "upgrade" threshold: how much better a galaxy
    # fit has to be versus ptsrc, and comp versus galaxy.
    cut = galaxy_margin

    # This is the "fractional" upgrade threshold for ptsrc/simple->dev/exp:
    # 2% of ptsrc vs nothing
    fcut = 0.02 * chisqs['ptsrc']
    #print('Cut: max of', cut, 'and', fcut, ' (fraction of chisq_psf=%.1f)'
    # % chisqs['ptsrc'])
    cut = max(cut, fcut)

    expdiff = chisqs['exp'] - chisqs[keepmod]
    devdiff = chisqs['dev'] - chisqs[keepmod]

    #print('EXP vs', keepmod, ':', expdiff)
    #print('DEV vs', keepmod, ':', devdiff)

    if not (expdiff > cut or devdiff > cut):
        #print('Keeping', keepmod)
        return keepmod

    if expdiff > devdiff:
        #print('Upgrading from PTSRC to EXP: diff', expdiff)
        keepmod = 'exp'
    else:
        #print('Upgrading from PTSRC to DEV: diff', expdiff)
        keepmod = 'dev'

    if not 'comp' in chisqs:
        return keepmod

    diff = chisqs['comp'] - chisqs[keepmod]
    #print('Comparing', keepmod, 'to comp.  cut:', cut, 'comp:', diff)
    if diff < cut:
        return keepmod

    #print('Upgrading from dev/exp to composite.')
    keepmod = 'comp'
    return keepmod


def _chisq_improvement(src, chisqs, chisqs_none):
    '''
    chisqs, chisqs_none: dict of band->chisq
    '''
    bright = src.getBrightness()
    bands = chisqs.keys()
    fluxes = dict([(b, bright.getFlux(b)) for b in bands])
    dchisq = 0.
    for b in bands:
        flux = fluxes[b]
        if flux == 0:
            continue
        # this will be positive for an improved model
        d = chisqs_none[b] - chisqs[b]
        if flux > 0:
            dchisq += d
        else:
            dchisq -= np.abs(d)
    return dchisq

def _per_band_chisqs(srctractor, bands):
    chisqs = dict([(b,0) for b in bands])
    for img in srctractor.images:
        chi = srctractor.getChiImage(img=img)
        chisqs[img.band] = chisqs[img.band] + (chi ** 2).sum()
    return chisqs

def _fit_fluxes(cat, tims, bands, use_ceres, alphas):
    # Try fitting fluxes first?
    cat.thawAllRecursive()
    for src in cat:
        src.freezeAllBut('brightness')
    for b in bands:
        for src in cat:
            src.getBrightness().freezeAllBut(b)
        # Images for this band
        btims = [tim for tim in tims if tim.band == b]

        btr = Tractor(btims, cat)
        btr.freezeParam('images')
        done = False
        if use_ceres:
            from tractor.ceres_optimizer import CeresOptimizer
            orig_opt = btr.optimizer
            btr.optimizer = CeresOptimizer(BW=8, BH=8)
            try:
                btr.optimize_forced_photometry(shared_params=False,
                                               wantims=False)
                done = True
            except:
                import traceback
                print('Warning: Optimize_forced_photometry with Ceres failed:')
                traceback.print_exc()
                print('Falling back to LSQR')
            btr.optimizer = orig_opt
        if not done:
            try:
                btr.optimize_forced_photometry(
                    alphas=alphas, shared_params=False, wantims=False)
            except:
                import traceback
                print('Warning: Optimize_forced_photometry failed:')
                traceback.print_exc()
                # carry on

