#!/usr/bin/env python

import math
import array
import ROOT
import lxdata


mjd_min = 55000
mjd_max = 57200
#mjd_step = 0.1
mjd_step = 0.5
nbuf = int((mjd_max-mjd_min)/mjd_step)

lxdat_gsc = lxdata.Lxdata()
lxdat_gsc.init_buf(mjd_min, mjd_max, nbuf)

#lxdat_gsc.read_fitsfile("glcbin24.0h_regbg_hv0.fits", "LCDAT_PIBAND3", "MJD", "RATE")
lxdat_gsc.read_fitsfile("glcbin24.0h_regbg_hv0.fits", "LCDAT_PIBAND3", "MJD", "RATE", "RERR", 1.0)
lxdat_gsc.set_integlxpow(0.0, 0.0, 6./7.)


vmjd = array.array('d', [0.0]*nbuf)
vlx  = array.array('d', [0.0]*nbuf)
vinteglx = array.array('d', [0.0]*nbuf)
vzero = array.array('d', [0.0]*nbuf)

for idx in range(nbuf) :
    vmjd[idx] = lxdat_gsc.get_mjd(idx)
    vlx[idx]  = lxdat_gsc.get_lx(idx)
    vinteglx[idx] = -lxdat_gsc.get_integlxpow(idx)
    

grph_lx = ROOT.TGraph(nbuf, vmjd, vlx)
grph_inlx = ROOT.TGraph(nbuf, vmjd, vinteglx)

lxmin = min(vlx)
lxmax = max(vlx)
inlxmin = min(vinteglx)
inlxmax = max(vinteglx)


ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1", "", 600, 600)
c1.Draw()

padlist = []
hlist = []
xmin = math.floor(mjd_min)
xmax = math.ceil(mjd_max)
nx = int(xmax-xmin)

npad = 2
pxmin = 0.0
pxmax = 1.0


### pad0
c1.cd()
ipad = 0
padname = "pad%d"%(ipad)
pymin = float(npad-ipad-1)/npad
pymax = float(npad-ipad)/npad
padlist.append(ROOT.TPad(padname, "", pxmin, pymin, pxmax, pymax))
padlist[-1].Draw()
padlist[-1].cd()
hname = "hdum%d"%(ipad)
hlist.append(ROOT.TH1F(hname, "", nx, xmin, xmax))
ymin = lxmin - (lxmax-lxmin)*0.05
ymax = lxmax + (lxmax-lxmin)*0.05
hlist[-1].SetMinimum(ymin)
hlist[-1].SetMaximum(ymax)
hlist[-1].Draw()
grph_lx.Draw("Lsame")

### pad1
c1.cd()
ipad = 1
padname = "pad%d"%(ipad)
pymin = float(npad-ipad-1)/npad
pymax = float(npad-ipad)/npad
padlist.append(ROOT.TPad(padname, "", pxmin, pymin, pxmax, pymax))
padlist[-1].Draw()
padlist[-1].cd()
hname = "hdum%d"%(ipad)
hlist.append(ROOT.TH1F(hname, "", nx, xmin, xmax))
ymin = inlxmin - (inlxmax-inlxmin)*0.05
ymax = inlxmax + (inlxmax-inlxmin)*0.05
hlist[-1].SetMinimum(ymin)
hlist[-1].SetMaximum(ymax)
hlist[-1].Draw()
grph_inlx.Draw("Lsame")

