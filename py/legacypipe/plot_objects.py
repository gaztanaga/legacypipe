import fitsio
import matplotlib.pyplot as plt
import sys

#def se_ellipse_func(x,y,se):
#    return se['CXXWIN_IMAGE']*(x-XWIN_IMAGE)**2+/
#            se['CYYWIN_IMAGE']*(y-YWIN_IMAGE)**2+/
#            se['CXYWIN_IMAGE']*(x-XWIN_IMAGE)*(y-YWIN_IMAGE) -3**2

def se_objects(se_cat_fn):
    a=fitsio.FITS(se_cat_fn)
    a2=a[2].read()
    keys=a[2].get_colnames()
    for x,y,a,b,theta in zip(a2['XWIN_IMAGE'],a2['YWIN_IMAGE'],a2['AWIN_IMAGE'],\
                            a2['BWIN_IMAGE'],a2['THETAWIN_IMAGE']):
        

if __name__ == '__main__':
    
    sys.exit(main())
  