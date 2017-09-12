#!/usr/bin/env python

import lxdata


mjd_min = 55000
mjd_max = 57200
#mjd_step = 0.1
mjd_step = 0.5
nbuf = int((mjd_max-mjd_min)/mjd_step)

lxdat_gsc = lxdata.Lxdata()
lxdat_gsc.init_buf(mjd_min, mjd_max, nbuf)

distance = 10.0
convfact = 2.0
f2lx = convfact*distance**2*0.12 ### e37 erg
snr_llim = -1.

#lxdat_gsc.read_fitsfile("glcbin24.0h_regbg_hv0.fits", "LCDAT_PIBAND3", "MJD", "RATE")
lxdat_gsc.read_gsclcfile("glcbin24.0h_regbg_hv0.fits", "LCDAT_PIBAND4", snr_llim, f2lx)

lxdat_gsc.set_integlxpow(0.0, 0.0, 6./7.)

for idx in range(nbuf) :
    print lxdat_gsc.get_mjd(idx), lxdat_gsc.get_lx(idx), lxdat_gsc.get_integlxpow(idx)
    

