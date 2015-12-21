# -*- coding: utf-8 -*-
"""test_script_fr2
A script of VV's potentiality flow routing method.
This version attempts to weight the potentials by slopes, equivalent to varying
the depth of the box into which the "injection molding" occurs.

Created on Fri Feb 20 13:45:52 2015

@author: danhobley
"""
from __future__ import print_function

from six.moves import range

#from landlab import RasterModelGrid
#from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel, quiver
from matplotlib.ticker import MaxNLocator

sqrt = np.sqrt

n = 50
#mg = RasterModelGrid(n, n, 1.)
nt = 10000
width = 1.
p_thresh = 0.000001
diffusivity_offsetter = 100000000.
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
fixdis = 0.4 #this is the p contour at which flow is "pinned"
variable_K = False

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
theta_vE = np.zeros_like(hR)
theta_vW = np.zeros_like(hR)
theta_vN = np.zeros_like(hR)
theta_vS = np.zeros_like(hR)
vmagE = np.zeros_like(hR)
vmagW = np.zeros_like(hR)
vmagN = np.zeros_like(hR)
vmagS = np.zeros_like(hR)
uE = np.zeros_like(hR)
uW = np.zeros_like(hR)
uN = np.zeros_like(hR)
uS = np.zeros_like(hR)

Wchanged = np.zeros_like(hR, dtype=bool)
Echanged = np.zeros_like(Wchanged)
Nchanged = np.zeros_like(Wchanged)
Schanged = np.zeros_like(Wchanged)

#uval = np.zeros_like(hR)
#vval = np.zeros_like(hR)


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
    pgradEx[core] = uE[core] #pgrad is VV's vv, a velocity
    pgradEy[core] = uN[core]+uS[core]+uN[Es]+uS[Es]
    pgradEy[core] *= 0.25
    vmagE[core] = sqrt(np.square(pgradEx[core])+np.square(pgradEy[core]))
    #now resolve the effective flow magnitudes to downhill
    theta_vE[core] = np.arctan(np.fabs(pgradEy[core])/(np.fabs(pgradEx[core])+1.e-10))
    vmagE[core] *= np.cos(np.fabs(thetaE[core]-theta_vE[core]))
    qE[core] = np.sign(hgradEx[core])*vmagE[core]*(CslopeE[core]-slope).clip(0.)*np.cos(thetaE[core])
    #the clip should deal with the eastern edge, but return here to check if probs


    hgradWx[core] = (hR[Ws]-hR[core])#/width
    hgradWy[core] = hR[SWs]-hR[NWs]+hR[Ss]-hR[Ns]
    hgradWy[core] *= 0.25
    CslopeW[core] = sqrt(np.square(hgradWx[core])+np.square(hgradWy[core]))
    thetaW[core] = np.arctan(np.fabs(hgradWy[core])/(np.fabs(hgradWx[core])+1.e-10))
    pgradWx[core] = uW[core]#/width
    pgradWy[core] = uN[core]+uS[core]+uN[Ws]+uS[Ws]
    pgradWy[core] *= 0.25
    vmagW[core] = sqrt(np.square(pgradWx[core])+np.square(pgradWy[core]))
    theta_vW[core] = np.arctan(np.fabs(pgradWy[core])/(np.fabs(pgradWx[core])+1.e-10))
    vmagW[core] *= np.cos(np.fabs(thetaW[core]-theta_vW[core]))
    qW[core] = np.sign(hgradWx[core])*vmagW[core]*(CslopeW[core]-slope).clip(0.)*np.cos(thetaW[core])

    hgradNx[core] = hR[NWs]-hR[NEs]+hR[Ws]-hR[Es]
    hgradNx[core] *= 0.25
    hgradNy[core] = (hR[core]-hR[Ns])#/width
    CslopeN[core] = sqrt(np.square(hgradNx[core])+np.square(hgradNy[core]))
    thetaN[core] = np.arctan(np.fabs(hgradNy[core])/(np.fabs(hgradNx[core])+1.e-10))
    pgradNx[core] = uE[core]+uW[core]+uE[Ns]+uW[Ns]
    pgradNx[core] *= 0.25
    pgradNy[core] = uN[core]#/width
    vmagN[core] = sqrt(np.square(pgradNx[core])+np.square(pgradNy[core]))
    theta_vN[core] = np.arctan(np.fabs(pgradNy[core])/(np.fabs(pgradNx[core])+1.e-10))
    vmagN[core] *= np.cos(np.fabs(thetaN[core]-theta_vN[core]))
    qN[core] = np.sign(hgradNy[core])*vmagN[core]*(CslopeN[core]-slope).clip(0.)*np.sin(thetaN[core])

    hgradSx[core] = hR[SWs]-hR[SEs]+hR[Ws]-hR[Es]
    hgradSx[core] *= 0.25
    hgradSy[core] = (hR[Ss]-hR[core])#/width
    CslopeS[core] = sqrt(np.square(hgradSx[core])+np.square(hgradSy[core]))
    thetaS[core] = np.arctan(np.fabs(hgradSy[core])/(np.fabs(hgradSx[core])+1.e-10))
    pgradSx[core] = uE[core]+uW[core]+uE[Ss]+uW[Ss]
    pgradSx[core] *= 0.25
    pgradSy[core] = uS[core]#/width
    vmagS[core] = sqrt(np.square(pgradSx[core])+np.square(pgradSy[core]))
    theta_vS[core] = np.arctan(np.fabs(pgradSy[core])/(np.fabs(pgradSx[core])+1.e-10))
    vmagS[core] *= np.cos(np.fabs(thetaS[core]-theta_vS[core]))
    qS[core] = np.sign(hgradSy[core])*vmagS[core]*(CslopeS[core]-slope).clip(0.)*np.sin(thetaS[core])

    hR[core] += dtwidth*(qS[core]+qW[core]-qN[core]-qE[core]+qsourceR[core])


    ###P SOLVER

    #mask for which core nodes get updated:
    mask = (hR[core]<p_thresh)
    #not_mask = np.logical_not(mask)
    if variable_K:
        kE = 1.+np.fabs(hR[Es]-hR[core])#+diffusivity_offsetter)#/diffusivity_offsetter
        kW = 1.+np.fabs(hR[Ws]-hR[core])#+diffusivity_offsetter)#/diffusivity_offsetter
        kS = 1.+np.fabs(hR[Ss]-hR[core])#+diffusivity_offsetter)#/diffusivity_offsetter
        kN = 1.+np.fabs(hR[Ns]-hR[core])#+diffusivity_offsetter)#/diffusivity_offsetter
    else:
        kE = 1.
        kW = 1.
        kN = 1.
        kS = 1.

    Wchanged[core] = np.less(hR[Ws]+hR[core],fixdis)
    Echanged[core] = np.less(hR[Es]+hR[core],fixdis)
    Nchanged[core] = np.less(hR[Ns]+hR[core],fixdis)
    Schanged[core] = np.less(hR[Ss]+hR[core],fixdis)

    for j in range(10):

        uW[Wchanged] = kW*(pR[Ws]-pR[core])[Wchanged[core]]
        uE[Echanged] = kE*(pR[Es]-pR[core])[Echanged[core]]
        uN[Nchanged] = kN*(pR[Ns]-pR[core])[Nchanged[core]]
        uS[Schanged] = kS*(pR[Ss]-pR[core])[Schanged[core]]

        pR[core] += uW[core]
        pR[core] += uS[core]
        pR[core] -= uE[core]
        pR[core] -= uN[core]
        pR[core] += qspR[core]
        pR[core] /= kN+kS+kE+kW

        pR[core][mask] = 0.

X,Y = np.meshgrid(np.arange(n),np.arange(n))
uval = uW[core]+uE[core]
vval = uN[core]+uS[core]
#imshow_node_grid(mg, h)
figure(1)
f1 = imshow(hR[core])
figure(2)
f2 = contour(X,Y,hR[core], locator=MaxNLocator(nbins=100))
clabel(f2)
quiver(X,Y,uval,vval)
figure(3)
f3 = contour(X,Y,pR[core], locator=MaxNLocator(nbins=100))
clabel(f3)
quiver(X,Y,uval,vval)
figure(4)
contour(X,Y,hR[core], locator=MaxNLocator(nbins=100))
contour(X,Y,pR[core], locator=MaxNLocator(nbins=100))
