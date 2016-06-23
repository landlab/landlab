# -*- coding: utf-8 -*-
"""
A script of VV's potentiality flow routing method.

Created on Fri Feb 20 13:45:52 2015

@author: danhobley
"""
from __future__ import print_function

from six.moves import range

#from landlab import RasterModelGrid
#from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel
from matplotlib.ticker import MaxNLocator

sqrt = np.sqrt

n = 50
#mg = RasterModelGrid(n, n, 1.)
nt = 20000
width = 1.
p_thresh = 0.000001
core = (slice(1,-1),slice(1,-1))
timp=np.zeros((nt, 2), dtype=int)
dtwidth = 0.2
hR = np.zeros((n+2,n+2), dtype=float)
pR=np.zeros_like(hR)
p=pR.view().ravel()
qsourceR=np.zeros_like(hR)
qsource=qsourceR[core].view().ravel()
qspR=np.zeros_like(hR)
qsp = qspR[core].view().ravel()
qsourceR[core][0,-1]=.9*sqrt(2.)*.33
qsourceR[core][0,0]=.9*sqrt(2.)
qsourceR[core][n-1,n//2-1]=1
qspR[core][0,-1]=sqrt(2)
qspR[core][0,0]=sqrt(2)
qspR[core][n-1,n//2-1]=1
slope=0.1

hgradEx = np.zeros_like(hR)
hgradWx = np.zeros_like(hR)
hgradNx = np.zeros_like(hR)
hgradSx = np.zeros_like(hR)
pgradEx = np.zeros_like(hR)
pgradWx = np.zeros_like(hR)
pgradNx = np.zeros_like(hR)
pgradSx = np.zeros_like(hR)
hgradEy = np.zeros_like(hR)
hgradWy = np.zeros_like(hR)
hgradNy = np.zeros_like(hR)
hgradSy = np.zeros_like(hR)
pgradEy = np.zeros_like(hR)
pgradWy = np.zeros_like(hR)
pgradNy = np.zeros_like(hR)
pgradSy = np.zeros_like(hR)
CslopeE = np.zeros_like(hR)
CslopeW = np.zeros_like(hR)
CslopeN = np.zeros_like(hR)
CslopeS = np.zeros_like(hR)
thetaE = np.zeros_like(hR)
thetaW = np.zeros_like(hR)
thetaN = np.zeros_like(hR)
thetaS = np.zeros_like(hR)
vE = np.zeros_like(hR)
vW = np.zeros_like(hR)
vN = np.zeros_like(hR)
vS = np.zeros_like(hR)

#set up slice offsets:
Es = (slice(1,-1),slice(2,n+2))
NEs = (slice(2,n+2),slice(2,n+2))
Ns = (slice(2,n+2),slice(1,-1))
NWs = (slice(2,n+2),slice(0,-2))
Ws = (slice(1,-1),slice(0,-2))
SWs = (slice(0,-2),slice(0,-2))
Ss = (slice(0,-2),slice(1,-1))
SEs = (slice(0,-2),slice(2,n+2))

for i in range(nt):
    if i%100==0:
        print(i)
    qE = np.zeros_like(hR)
    qW = np.zeros_like(hR)
    qN = np.zeros_like(hR)
    qS = np.zeros_like(hR)

    #update the dummy edges of our variables:
    hR[0,1:-1] = hR[1,1:-1]
    hR[-1,1:-1] = hR[-2,1:-1]
    hR[1:-1,0] = hR[1:-1,1]
    hR[1:-1,-1] = hR[1:-1,-2]
    pR[0,1:-1] = pR[1,1:-1]
    pR[-1,1:-1] = pR[-2,1:-1]
    pR[1:-1,0] = pR[1:-1,1]
    pR[1:-1,-1] = pR[1:-1,-2]

    hgradEx[core] = (hR[core]-hR[Es])#/width
    hgradEy[core] = hR[SEs]-hR[NEs]+hR[Ss]-hR[Ns]
    hgradEy[core] *= 0.25
    CslopeE[core] = sqrt(np.square(hgradEx[core])+np.square(hgradEy[core]))
    thetaE[core] = np.arctan(np.fabs(hgradEy[core])/(np.fabs(hgradEx[core])+1.e-10))
    pgradEx[core] = (pR[core]-pR[Es])#/width
    pgradEy[core] = pR[SEs]-pR[NEs]+pR[Ss]-pR[Ns]
    pgradEy[core] *= 0.25
    vE[core] = sqrt(np.square(pgradEx[core])+np.square(pgradEy[core]))
    qE[core] = np.sign(hgradEx[core])*vE[core]*(CslopeE[core]-slope).clip(0.)*np.cos(thetaE[core])
###the clip should deal with the eastern edge, but return here to check if probs

    hgradWx[core] = (hR[Ws]-hR[core])#/width
    hgradWy[core] = hR[SWs]-hR[NWs]+hR[Ss]-hR[Ns]
    hgradWy[core] *= 0.25
    CslopeW[core] = sqrt(np.square(hgradWx[core])+np.square(hgradWy[core]))
    thetaW[core] = np.arctan(np.fabs(hgradWy[core])/(np.fabs(hgradWx[core])+1.e-10))
    pgradWx[core] = (pR[Ws]-pR[core])#/width
    pgradWy[core] = pR[SWs]-pR[NWs]+pR[Ss]-pR[Ns]
    pgradWy[core] *= 0.25
    vW[core] = sqrt(np.square(pgradWx[core])+np.square(pgradWy[core]))
    qW[core] = np.sign(hgradWx[core])*vW[core]*(CslopeW[core]-slope).clip(0.)*np.cos(thetaW[core])

    hgradNx[core] = hR[NWs]-hR[NEs]+hR[Ws]-hR[Es]
    hgradNx[core] *= 0.25
    hgradNy[core] = (hR[core]-hR[Ns])#/width
    CslopeN[core] = sqrt(np.square(hgradNx[core])+np.square(hgradNy[core]))
    thetaN[core] = np.arctan(np.fabs(hgradNy[core])/(np.fabs(hgradNx[core])+1.e-10))
    pgradNx[core] = pR[NWs]-pR[NEs]+pR[Ws]-pR[Es]
    pgradNx[core] *= 0.25
    pgradNy[core] = (pR[core]-pR[Ns])#/width
    vN[core] = sqrt(np.square(pgradNx[core])+np.square(pgradNy[core]))
    qN[core] = np.sign(hgradNy[core])*vN[core]*(CslopeN[core]-slope).clip(0.)*np.sin(thetaN[core])

    hgradSx[core] = hR[SWs]-hR[SEs]+hR[Ws]-hR[Es]
    hgradSx[core] *= 0.25
    hgradSy[core] = (hR[Ss]-hR[core])#/width
    CslopeS[core] = sqrt(np.square(hgradSx[core])+np.square(hgradSy[core]))
    thetaS[core] = np.arctan(np.fabs(hgradSy[core])/(np.fabs(hgradSx[core])+1.e-10))
    pgradSx[core] = pR[SWs]-pR[SEs]+pR[Ws]-pR[Es]
    pgradSx[core] *= 0.25
    pgradSy[core] = (pR[Ss]-pR[core])#/width
    vS[core] = sqrt(np.square(pgradSx[core])+np.square(pgradSy[core]))
    qS[core] = np.sign(hgradSy[core])*vS[core]*(CslopeS[core]-slope).clip(0.)*np.sin(thetaS[core])

    hR[core] += dtwidth*(qS[core]+qW[core]-qN[core]-qE[core]+qsourceR[core])

    #while 1:
    #mask for which core nodes get updated:
    mask = (hR[core]<p_thresh)
    for j in range(100):
        pR[core] = pR[Ns]+pR[Ss]+pR[Es]+pR[Ws]+qspR[core]
        pR[core] *= 0.25

        pR[core][mask] = 0.

X,Y = np.meshgrid(np.arange(n),np.arange(n))
#imshow_node_grid(mg, h)
figure(1)
f1 = imshow(hR[core])
figure(2)
f2 = contour(X,Y,hR[core], locator=MaxNLocator(nbins=100))
clabel(f2)
figure(3)
f3 = contour(X,Y,pR[core], locator=MaxNLocator(nbins=100))
clabel(f3)
figure(4)
contour(X,Y,hR[core], locator=MaxNLocator(nbins=100))
contour(X,Y,pR[core], locator=MaxNLocator(nbins=100))
