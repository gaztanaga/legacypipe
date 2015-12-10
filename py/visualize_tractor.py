from matplotlib.patches import Ellipse,Circle,Rectangle
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from astropy.io import fits
import numpy as np
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import argparse

parser = argparse.ArgumentParser(description="test")
parser.add_argument("-imgfn",action="store",help='abs path to image file name')
parser.add_argument("-maskfn",action="store",help='...mask-2')
parser.add_argument("-sefn",action="store",help='...SExtractor catalogue')
parser.add_argument("-psfexfn",action="store",help='...PSFex catalogue')
parser.add_argument("-modelfn",action="store",help='...Tractor model')
parser.add_argument("-tractor_imgfn",action="store",help='...Tractor model')
args = parser.parse_args()

##PSFEX
class KaylanPsfex(object):
    '''
    An object representing a PsfEx PSF model.
    '''
    def __init__(self, psfexfn):
        pf = fits.open(psfexfn, memmap=True)
        ims = pf[1].data['PSF_MASK']
        ims=ims[0]
        hdr = pf[1].header
        print('Got', ims.shape, 'PSF images')
        # PSF distortion bases are polynomials of x,y
        assert(hdr['POLNAME1'].strip() == 'X_IMAGE')
        assert(hdr['POLNAME2'].strip() == 'Y_IMAGE')
        assert(hdr['POLGRP1'] == 1)
        assert(hdr['POLGRP2'] == 1)
        assert(hdr['POLNGRP' ] == 1)
        x0     = hdr.get('POLZERO1')
        xscale = hdr.get('POLSCAL1')
        y0     = hdr.get('POLZERO2')
        yscale = hdr.get('POLSCAL2')
        degree = hdr.get('POLDEG1')
        self.sampling = hdr.get('PSF_SAMP')
        print('PsfEx sampling:', self.sampling)
        # number of terms in polynomial
        ne = (degree + 1) * (degree + 2) / 2
        assert(hdr['PSFAXIS3'] == ne)
        assert(len(ims.shape) == 3)
        assert(ims.shape[0] == ne)
        self.psfbases = ims
        self.xscale, self.yscale = xscale, yscale
        self.degree = degree
        print('PsfEx degree:', self.degree)
        bh,bw = self.psfbases[0].shape
        self.radius = (bh+1)/2.
        self.x0,self.y0 = x0,y0

    def polynomials(self, x, y, powers=False):
        dx = (x - self.x0) / self.xscale
        dy = (y - self.y0) / self.yscale
        nb,h,w = self.psfbases.shape
        terms = np.zeros(nb)

        if powers:
            xpows = np.zeros(nb, int)
            ypows = np.zeros(nb, int)

        for d in range(self.degree + 1):
            # x polynomial degree = j
            # y polynomial degree = k
            for j in range(d+1):
                k = d - j
                amp = dx**j * dy**k
                # PSFEx manual pg. 111 ?
                ii = j + (self.degree+1) * k - (k * (k-1))/ 2
                #print('getPolynomialTerms: j=', j, 'k=', k, 'd=', d, 'ii=', ii)
                # It goes: order 0, order 1, order 2, ...
                # and then j=0, j=1, ...
                terms[ii] = amp
                if powers:
                    xpows[ii] = j
                    ypows[ii] = k
        if powers:
            return (terms, xpows, ypows)
        return terms

    def at(self, x, y, nativeScale=True):
        '''
        Returns an image of the PSF at the given pixel coordinates.
        '''
        psf = np.zeros_like(self.psfbases[0])

        #print('Evaluating PsfEx at', x,y)
        for term,base in zip(self.polynomials(x,y), self.psfbases):
            #print('  polynomial', term, 'x base w/ range', base.min(), base.max())
            psf += term * base

        if nativeScale and self.sampling != 1:
            from scipy.ndimage.interpolation import affine_transform
            ny,nx = psf.shape
            spsf = affine_transform(psf, [1./self.sampling]*2,
                                    offset=nx/2 * (self.sampling - 1.))
            return spsf
            
        return psf


    def plot_bases(self, autoscale=True, stampsize=None):
        '''6 panel plot of the 6 psfex basis vectors'''
        N = len(self.psfbases)
        cols = int(np.ceil(np.sqrt(N)))
        rows = int(np.ceil(N / float(cols)))
        plt.clf()
        plt.subplots_adjust(hspace=0, wspace=0)

        cut = 0
        if stampsize is not None:
            H,W = self.shape
            assert(H == W)
            cut = max(0, (H - stampsize) / 2)

        ima = dict(interpolation='nearest', origin='lower',
                  cmap=plt.get_cmap('gray'))
        if autoscale:
            mx = self.psfbases.max()
            ima.update(vmin=-mx, vmax=mx)
        nil, xpows, ypows = self.polynomials(0., 0., powers=True)
        for i,(xp,yp,b) in enumerate(zip(xpows, ypows, self.psfbases)):
            plt.subplot(rows, cols, i+1)

            if cut > 0:
                b = b[cut:-cut, cut:-cut]
            if autoscale:
                plt.imshow(b, **ima)
            else:
                mx = np.abs(b).max()
                plt.imshow(b, vmin=-mx, vmax=mx, **ima)
            plt.xticks([])
            plt.yticks([])
            plt.title('x^%i y^%i' % (xp,yp))
        plt.suptitle('PsfEx eigen-bases')

##image

########
## copied from skimage/exposure/exposure.py and skimage/util/dtype.py because
## not easy to install scikit-image on NERSC computers
def intensity_range(image, range_values='image', clip_negative=False):
    DTYPE_RANGE = {np.bool_: (False, True),
               np.bool8: (False, True),
               np.uint8: (0, 255),
               np.uint16: (0, 65535),
               np.int8: (-128, 127),
               np.int16: (-32768, 32767),
               np.int64: (-2**63, 2**63 - 1),
               np.uint64: (0, 2**64 - 1),
               np.int32: (-2**31, 2**31 - 1),
               np.uint32: (0, 2**32 - 1),
               np.float32: (-1, 1),
               np.float64: (-1, 1)}
    if range_values == 'dtype':
        range_values = image.dtype.type

    if range_values == 'image':
        i_min = np.min(image)
        i_max = np.max(image)
    elif range_values in DTYPE_RANGE:
        i_min, i_max = DTYPE_RANGE[range_values]
        if clip_negative:
            i_min = 0
    else:
        i_min, i_max = range_values
    return i_min, i_max

def rescale_intensity(image, in_range='image', out_range='dtype'):
    dtype = image.dtype.type

    if in_range is None:
        in_range = 'image'
        msg = "`in_range` should not be set to None. Use {!r} instead."
        warnings.warn(msg.format(in_range))

    if out_range is None:
        out_range = 'dtype'
        msg = "`out_range` should not be set to None. Use {!r} instead."
        warnings.warn(msg.format(out_range))

    imin, imax = intensity_range(image, in_range)
    omin, omax = intensity_range(image, out_range, clip_negative=(imin >= 0))

    image = np.clip(image, imin, imax)

    image = (image - imin) / float(imax - imin)
    return dtype(image * (omax - omin) + omin)
###########
###########

def read_invvar(imgfn,dqfn, clip=False, clipThresh=0.2, **kwargs):
	img= fits.open(imgfn)[0].data
	dq= fits.open(dqfn)[0].data 
	assert(dq.shape == img.shape)
	invvar=np.zeros(img.shape)
	invvar[dq == 0]= np.power(img[dq == 0],-0.5)
	#invvar[dq == 2]= np.power(img[dq == 2],-0.5) #using mask-2, SExtractor ojbect if binary(00010) = 2
	if clip:
	    # Clamp near-zero (incl negative!) invvars to zero.
	    # These arise due to fpack.
	    if clipThresh > 0.:
		med = np.median(invvar[invvar > 0])
		thresh = clipThresh * med
	    else:
		thresh = 0.
	    invvar[invvar < thresh] = 0
	return invvar


def stretch_min_max(img):
    p2, p98 = np.percentile(img, (2, 98))
    return rescale_intensity(img, in_range=(p2, p98))

def stretch_hist_eq(img):
    return exposure.equalize_hist(img)

def plot_images(imgfn,maskfn,xyc):
    fig= plt.figure()
    fig.set_size_inches(20.,20.)
    grid = AxesGrid(fig, 111,  # similar to subplot(143)
                    nrows_ncols=(1, 3),
                    axes_pad=0.2,
                    label_mode="all",
                    share_all=False,
                    cbar_location="bottom",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )
    clim=[]
    for i,fn,title in zip(range(3),[imgfn,maskfn,'junk'],\
                          ['Image (stretched)','Mask (0=good, 1=any flag)','Inverse Variance (stretched)']):
        if i == 2: img= read_invvar(imgfn,maskfn) 
        else: img=fits.open(fn)[0].data
        if i == 1:
            img[img > 0]=1
            print "mask min,max= ",img.min(),img.max()
        else: img= stretch_min_max(img)
        clim.append([img.min(),img.max()])
        im = grid[i].imshow(img,cmap=plt.get_cmap('gray'),origin='lower',\
                            vmin=clim[i][0], vmax=clim[i][1])
        if i == 0:
            for c in xyc:
                bott_left= np.array(c)-100
                rect = Rectangle(bott_left,200,200)
                grid[i].add_artist(rect)
                rect.set_alpha(0.8)
                rect.set_facecolor('none')
                rect.set_edgecolor('b')
        grid[i].set_title(title,fontsize='xx-large')
        if i != 0:
            grid[i].set_xticks([])
            grid[i].set_yticks([])
        grid.cbar_axes[i].colorbar(im)
        grid.cbar_axes[i].set_xticks([clim[i][0], clim[i][1]])
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
    plt.savefig('image_objects.png',dpi=150)    

def tractor_residual(imagefn,modelfn):
    img=fits.open(imagefn)[0].data
    mod=fits.open(modelfn)[0].data
    return img - mod

def plot_tractor(imgfn,modelfn,name=None):
    '''file names are to fits files'''
    fig= plt.figure()
    fig.set_size_inches(20.,20.)
    grid = AxesGrid(fig, 111,  # similar to subplot(143)
                    nrows_ncols=(1, 3),
                    axes_pad=0.2,
                    label_mode="all",
                    share_all=False,
                    cbar_location="bottom",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )
    clim=[]
    for i,fn,title in zip(range(3),[imgfn,modelfn,'junk'],\
                          ['Tractor Image (stretched)','Tractor Model (stretched)','Residual (Image - Model)']):
        if i == 2: img= tractor_residual(imgfn,modelfn) 
        else: img=fits.open(fn)[0].data
        img= stretch_min_max(img)
        clim.append([img.min(),img.max()])
        im = grid[i].imshow(img,cmap=plt.get_cmap('gray'),origin='lower',\
                            vmin=clim[i][0], vmax=clim[i][1])
        grid[i].set_title(title,fontsize='xx-large')
        if i != 0:
            grid[i].set_xticks([])
            grid[i].set_yticks([])
        grid.cbar_axes[i].colorbar(im)
        grid.cbar_axes[i].set_xticks([clim[i][0], clim[i][1]])
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
    if name != None: plt.savefig(name,dpi=150)
    else: plt.savefig('tractor_images.png',dpi=150)    

def plot_tractor_cutouts(imgfn,modelfn,xyc,wh=200):
    #model
    mod= stretch_min_max( fits.open(modelfn)[0].data )
    resid= tractor_residual(imgfn,modelfn)
    resid= stretch_min_max( resid)
    #plot
    fig,axes = plt.subplots(5,6)
    fig.set_size_inches(15.,20.)
    plt.subplots_adjust(hspace=-0.5, wspace=0)
    se_rows,se_cols= 5,3
    #model
    cnt=0
    for r,rpanel in zip(range(se_rows),range(se_rows)[::-1]):  
        yc= xyc[cnt][1]
        for c,cpanel in zip(range(se_cols),np.arange(se_cols)): 
            xc= xyc[cnt][0] 
            xlim,ylim= [xc-wh/2,xc+wh/2],[yc-wh/2,yc+wh/2]
            ax = axes[rpanel,cpanel] #plt.subplot2grid((5, 6), (rpanel, cpanel))
            ax.set_xlim(xlim) #CRITICAL for ellipses to appear
            ax.set_ylim(ylim)
            ax.imshow(mod, cmap=plt.get_cmap('gray'),origin='lower')
            ax.set_xticks([])
            ax.set_yticks([])
            cnt+=1
    #plot PSFex
    cnt=0
    for rpanel in range(se_rows)[::-1]:  
        yc= xyc[cnt][1]
        for cpanel in np.arange(se_cols)+3: 
            xc= xyc[cnt][0] 
            xlim,ylim= [xc-wh/2,xc+wh/2],[yc-wh/2,yc+wh/2]
            ax = axes[rpanel,cpanel] #plt.subplot2grid((5, 6), (rpanel, cpanel))
            ax.set_xlim(xlim) #CRITICAL for ellipses to appear
            ax.set_ylim(ylim)
            ax.imshow(resid, cmap=plt.get_cmap('gray'),origin='lower')
            ax.set_xticks([])
            ax.set_yticks([])
            cnt+=1
    #save plot
    plt.savefig('tractor_cutouts.png',dpi=150)
    plt.close()
 
 
def plot_se(imgfn,maskfn,sefn,psfexfn,wh=200):
    #image
    tmp=fits.open(imgfn)
    img=tmp[0].data
    p2, p98 = np.percentile(img, (2, 98))
    img_rescale = rescale_intensity(img, in_range=(p2, p98))
    #se catalogue
    a=fits.open(sefn)
    a2=a[2].data
    cols= a[2].columns
    #plot
    fig,axes = plt.subplots(5,6)
    fig.set_size_inches(15.,20.)
    plt.subplots_adjust(hspace=-0.5, wspace=0)
#     plt.tight_layout()
    #plot sextractor objects 
    se_rows,se_cols= 5,3
    rowstep,colstep= img.shape[0]/(se_rows+1),img.shape[1]/(se_cols+1)
    xyc=[] #store pixel location of center of each postage stamp
    for r,rpanel in zip(range(se_rows),range(se_rows)[::-1]):  
        yc= rowstep*(r+1)
        for c,cpanel in zip(range(se_cols),np.arange(se_cols)): 
            xc= colstep*(c+1) 
            print "center= ",xc,yc
            xyc.append([xc,yc])#for later use
            xlim,ylim= [xc-wh/2,xc+wh/2],[yc-wh/2,yc+wh/2]
            ax = axes[rpanel,cpanel] #plt.subplot2grid((5, 6), (rpanel, cpanel))
            ax.set_xlim(xlim) #CRITICAL for ellipses to appear
            ax.set_ylim(ylim)
            ax.imshow(img_rescale, cmap=plt.get_cmap('gray'),origin='lower')
            ax.set_xticks([])
            ax.set_yticks([])
#             print "xlim= ",xlim,"ylim= ",ylim
            x_index=np.logical_and(a2['X_IMAGE'] > xlim[0],a2['X_IMAGE'] < xlim[1])
            y_index=np.logical_and(a2['Y_IMAGE'] > ylim[0],a2['Y_IMAGE'] < ylim[1])
            xy=np.logical_and(x_index,y_index)
            for i in range(len(a2['X_IMAGE'][xy])):
    #             ellip= Ellipse(xy= (a2['XWIN_IMAGE'][xy][i],a2['YWIN_IMAGE'][xy][i]),\
    #                                        width=2*a2['AWIN_IMAGE'][xy][i], \
    #                                        height=2*a2['BWIN_IMAGE'][xy][i],\
    #                                        angle=a2['THETAWIN_IMAGE'][xy][i])
                ellip= Ellipse(xy= (a2['X_IMAGE'][xy][i],a2['Y_IMAGE'][xy][i]),\
                                           width=2*a2['A_IMAGE'][xy][i], \
                                           height=2*a2['B_IMAGE'][xy][i],\
                                           angle=a2['THETA_IMAGE'][xy][i])
                ax.add_artist(ellip)
                ellip.set_alpha(0.8)
                ellip.set_facecolor('none')
                ellip.set_edgecolor('b')
    #plot PSFex
    psf= KaylanPsfex(psfexfn)
    cnt=0
    for rpanel in range(se_rows)[::-1]:  
        for cpanel in np.arange(se_cols)+3: 
            center= xyc[cnt]
            img=psf.at(center[0],center[1])
            p2, p98 = np.percentile(img, (2, 98))
            img_rescale = rescale_intensity(img, in_range=(p2, p98))
            ax = axes[rpanel,cpanel] #plt.subplot2grid((5, 6), (rpanel, cpanel))
            ax.imshow(img_rescale, cmap=plt.get_cmap('gray'),origin='lower')
            ax.set_xticks([])
            ax.set_yticks([])
            #
            cnt+=1
    #save plot
    plt.savefig('SExtractor_objects.png',dpi=150)
    plt.close()
    #so can overplot locations on images
    return xyc

   
# plot_img(imgfn)
#xyc= plot_se(args.imgfn,args.maskfn,args.sefn,args.psfexfn,wh=300)
#plot_images(args.imgfn,args.maskfn,xyc)
#plot_tractor(args.tractor_imgfn,args.modelfn,name='tractor2.png')
plot_tractor(args.imgfn,args.modelfn,name='tractor1.png')
#plot_tractor_cutouts(args.tractor_imgfn,args.modelfn,xyc,wh=300)

