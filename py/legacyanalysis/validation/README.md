# theValidator
THE software for validating Legacy Survey Data

theValidator is a small collection of python modules, fully integrated into ipython notebooks via NERSCs Jupyter Hub jupyter.nersc.gov, that can read and analyze [Legacy Survey](http://legacysurvey.org/) data products, such as Tractor Catalogues, concantenating them, matching two sets of catalogues, etc. 

Instructions
============

- ssh user@edison.nersc.gov (or Cori)
- git clone https://github.com/legacysurvey/legacypipe.git (do this in home directory)
- open browser, go to jupyter.nersc.gov, login
- navigate to $HOME/legacypipe/py/legacyanalysis/validation
- click on README.ipynb and follow instructions

Features
========

- Anyone with NERSC account can run at jupyter.nersc.gov 
- Simple dependencies -- matplotlib, numpy, scipy, astropy (no astrometry.net, tractor, legacypipe)
- Tested -- concatenating Tractor Catalogues with astropy identical results to Dustins fits_tables

Future
========

- healpix maps -- depth, seeing, airmass, etc (Ashely R. did some of this already) 
- psql -- example code for using the PostgreSQL DB at NERSC
- decals, mzls, bass -- telescope specific validation
- compact galaxies -- identify with dchisq + PSF 
- cosmic rays 
