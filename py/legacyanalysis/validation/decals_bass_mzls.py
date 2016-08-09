#!/usr/bin/env python

"""do DECaLS to BASS/MzLS comparison
"""

from __future__ import division, print_function

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import argparse
import pickle

from legacyanalysis.validation.pathnames import get_outdir
from legacyanalysis.validation.combine_cats import get_matched_dataset
import legacyanalysis.validation.common_plots as plots

parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='DECaLS simulations.')
parser.add_argument('--decals_list', type=str, help='list of tractor cats',required=True)
parser.add_argument('--bassmos_list', type=str, help='ditto',required=True) 
args = parser.parse_args()

# Get data for Matched sources
d= get_matched_dataset(args.decals_list, args.bassmos_list, \
                       comparison='bmd',debug=False)

# All matched 
d.apply_cut(['all'])
plots.nobs(d.ref.data['tractor'], outname=os.path.join(d.outdir,'nobs_decals.png'))
plots.nobs(d.test.data['tractor'], outname=os.path.join(d.outdir,'nobs_bass.png'))
plots.radec(d.ref.data['tractor'], outname=os.path.join(d.outdir,'radec_decals.png'))
plots.radec(d.test.data['tractor'], outname=os.path.join(d.outdir,'radec_bass.png'))

# Clean and psf 
d.apply_cut(['clean'])
plots.confusion_matrix(d.ref.data['tractor'],d.test.data['tractor'], outname=os.path.join(d.outdir,'conf.png'),\
                       ref_name='DECaLS',test_name='BASS_MzLS')
plots.stacked_confusion_matrix(d.ref.data['tractor'],d.test.data['tractor'],\
                               d.ref.data['extra'],d.test.data['extra'],\
                               ref_name='DECaLS',test_name='BASS_MzLS',\
                               outname=os.path.join(d.outdir,'conf_stacked.png'))

plots.delta_mag_vs_mag(d.ref.data['extra'],d.test.data['extra'], ylim=[-0.5,0.5],\
                       outname=os.path.join(d.outdir,'delta_mag.png'))
plots.chi_v_gaussian(d.ref.data['tractor'],d.test.data['tractor'],\
                     d.ref.data['extra'],d.test.data['extra'],\
                     low=-8.,hi=8., outname=os.path.join(d.outdir,'chi.png'))

#plots.nobs(d.test_matched, name='BASS_MzLS')
#plots.radec(d.test_missed, name='Missed_Test')
#plots.matched_dist(d.ref_matched,d.meta['d_matched'], name='')
#plots.hist_types(d.ref_matched, name='Ref')
#plots.hist_types(d.test_matched, name='Test')
#plots.n_per_deg2(d.ref_matched, deg2=d.meta['deg2']['ref'], name='Ref_MatchedDefault')
#plots.n_per_deg2(d.test_matched, deg2=d.meta['deg2']['test'], name='Test_MatchedDefault')
#
## Change mask to union of default masks
#d.ref_matched.add_to_current_mask(d.test_matched.masks['default'])
#d.test_matched.add_to_current_mask(d.ref_matched.masks['default'])
#plots.radec(d.ref_matched, name='Matched_DefaultUnion')
#plots.hist_types(d.ref_matched, name='Ref_DefaultUnion')
#plots.hist_types(d.test_matched, name='Test_DefaultUnion')
#
## union of default + PSF
#d.ref_matched.apply_mask_by_names(['current','psf'])
#d.test_matched.apply_mask_by_names(['current','psf'])
#
## Change mask to default + PSF
#d.ref_matched.apply_mask_by_names(['default','psf'])
#d.test_matched.apply_mask_by_names(['default','psf'])
#plots.sn_vs_mag(d.ref_matched, name="Ref_DefaultPsf")
#plots.sn_vs_mag(d.test_matched, name="Test_DefaultPsf")

print('finished DECaLS to BASS/MzLS comparison')
