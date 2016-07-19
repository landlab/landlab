# -*- coding: utf-8 -*-
"""test_script_fr3
A script of our potentiality "ghost field" flow routing method.

Created on Fri Feb 20 13:45:52 2015

@author: danhobley
"""
from __future__ import print_function

from six.moves import range

#from landlab import RasterModelGrid
#from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel, quiver
from landlab.plot import imshow_grid_at_node
from landlab import RasterModelGrid
from matplotlib.ticker import MaxNLocator

sqrt = np.sqrt

nrows = 50
ncols = 50
#mg = RasterModelGrid(n, n, 1.)
#nt = 13000
nt = 15500
width = 1.
slope=0.1
core = (slice(1,-1),slice(1,-1))
dtwidth = 0.2
hR = np.zeros((nrows+2,ncols+2), dtype=float)
qwater_inR=np.zeros_like(hR) #WATER
qsed_inR=np.zeros_like(hR) #SED
#qwater_inR[core][0,-1]=np.pi/2.
qwater_inR[core][0,0]=1.
qwater_inR[core][0,-1]=0.5
qwater_inR[core][-1,24]=1.
#qsourceR[core][0,0]=.9*sqrt(2.)
#qsourceR[core][-1,n//2-1]=1
#qspR[core][0,-1]=np.pi/2.*(1.-slope)
#qsed_inR[core][0,0]=np.pi/2.*(1.-slope)
qsed_inR[core][0,0]=1.
qsed_inR[core][0,-1]=0.5
qsed_inR[core][-1,24]=1.
#qspR[core][0,0]=sqrt(2)
#qspR[core][-1,n//2-1]=1
flat_threshold = 0.00001

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
#coeffs for K solver:
aPP = np.zeros_like(hR)
aWW = np.zeros_like(hR)
aWP = np.zeros_like(hR)
aEE = np.zeros_like(hR)
aEP = np.zeros_like(hR)
aNN = np.zeros_like(hR)
aNP = np.zeros_like(hR)
aSS = np.zeros_like(hR)
aSP = np.zeros_like(hR)

qsedE = np.zeros_like(hR)
qsedW = np.zeros_like(hR)
qsedN = np.zeros_like(hR)
qsedS = np.zeros_like(hR)

K = np.zeros_like(hR)
not_flat = np.zeros((nrows,ncols), dtype=bool)

#Wchanged = np.zeros_like(hR, dtype=bool)
#Echanged = np.zeros_like(Wchanged)
#Nchanged = np.zeros_like(Wchanged)
#Schanged = np.zeros_like(Wchanged)

#uxval = np.zeros_like(hR)
#uyval = np.zeros_like(hR)


#set up slice offsets:
Es = (slice(1,-1),slice(2,ncols+2))
NEs = (slice(2,nrows+2),slice(2,ncols+2))
Ns = (slice(2,nrows+2),slice(1,-1))
NWs = (slice(2,nrows+2),slice(0,-2))
Ws = (slice(1,-1),slice(0,-2))
SWs = (slice(0,-2),slice(0,-2))
Ss = (slice(0,-2),slice(1,-1))
SEs = (slice(0,-2),slice(2,ncols+2))

for i in range(nt):
    if i%100==0:
        print(i)
    qsedE.fill(0.)
    qsedW.fill(0.)
    qsedN.fill(0.)
    qsedS.fill(0.)

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
    qsedE[core] = np.sign(hgradEx[core])*vmagE[core]*(CslopeE[core]-slope).clip(0.)*np.cos(thetaE[core])
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
    qsedW[core] = np.sign(hgradWx[core])*vmagW[core]*(CslopeW[core]-slope).clip(0.)*np.cos(thetaW[core])

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
    qsedN[core] = np.sign(hgradNy[core])*vmagN[core]*(CslopeN[core]-slope).clip(0.)*np.sin(thetaN[core])

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
    qsedS[core] = np.sign(hgradSy[core])*vmagS[core]*(CslopeS[core]-slope).clip(0.)*np.sin(thetaS[core])

    hR[core] += dtwidth*(qsedS[core]+qsedW[core]-qsedN[core]-qsedE[core]+qsed_inR[core])

    #update the dummy edges of our variables:
    hR[0,1:-1] = hR[1,1:-1]
    hR[-1,1:-1] = hR[-2,1:-1]
    hR[1:-1,0] = hR[1:-1,1]
    hR[1:-1,-1] = hR[1:-1,-2]
    hR[(0,-1,0,-1),(0,-1,-1,0)] = hR[(1,-2,1,-2),(1,-2,-2,1)]

    ###P SOLVER

    #not_flat = np.greater(hR[core],flat_threshold)
    #not_mask = np.logical_not(mask)
    aNN[core] = (-hR[core]+hR[Ns]).clip(0.)
    aNP[core] =  (hR[core]-hR[Ns]).clip(0.)
    aSS[core] = (-hR[core]+hR[Ss]).clip(0.)
    aSP[core] =  (hR[core]-hR[Ss]).clip(0.)
    aEE[core] = (-hR[core]+hR[Es]).clip(0.)
    aEP[core] =  (hR[core]-hR[Es]).clip(0.)
    aWW[core] = (-hR[core]+hR[Ws]).clip(0.)
    aWP[core] =  (hR[core]-hR[Ws]).clip(0.)
    aPP[core] = aWP[core]+aEP[core]+aSP[core]+aNP[core]+1.e-6

    for j in range(15):

        #assert np.all(np.greater(aPP[core][not_flat],0.)) #this is here to eliminate a divby0
        #K[core][not_flat] = ((aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
        #                    +qwater_inR[core])[not_flat])/aPP[core][not_flat]
        K[core] = (aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
                            +qwater_inR[core])/aPP[core]

        for BC in (K,):
            BC[0,1:-1] = BC[1,1:-1]
            BC[-1,1:-1] = BC[-2,1:-1]
            BC[1:-1,0] = BC[1:-1,1]
            BC[1:-1,-1] = BC[1:-1,-2]
            BC[(0,-1,0,-1),(0,-1,-1,0)] = BC[(1,-2,1,-2),(1,-2,-2,1)]

        uW[core] = aWW[core]*K[Ws]-aWP[core]*K[core]
        uE[core] = -aEE[core]*K[Es]+aEP[core]*K[core]
        uN[core] = -aNN[core]*K[Ns]+aNP[core]*K[core]
        uS[core] = aSS[core]*K[Ss]-aSP[core]*K[core]
        #update the u BCs
        for BC in (uW,uE,uN,uS):
            BC[0,1:-1] = BC[1,1:-1]
            BC[-1,1:-1] = BC[-2,1:-1]
            BC[1:-1,0] = BC[1:-1,1]
            BC[1:-1,-1] = BC[1:-1,-2]
            BC[(0,-1,0,-1),(0,-1,-1,0)] = BC[(1,-2,1,-2),(1,-2,-2,1)]

X,Y = np.meshgrid(np.arange(ncols),np.arange(nrows))
uval = uW[core]+uE[core]
vval = uN[core]+uS[core]
#velmag = sqrt(uval**2 + vval**2)
#uval /= velmag
#vval /= velmag
#imshow_node_grid(mg, h)
figure(1)
mg = RasterModelGrid((nrows, ncols))
f1 = imshow_grid_at_node(mg, hR[core].flatten(), grid_units=('m', 'm'))
figure(2)
f2 = contour(X,Y,hR[core], locator=MaxNLocator(nbins=100))
# f2 = contour(X, Y, np.sqrt(uval**2+vval**2), locator=MaxNLocator(nbins=10))
clabel(f2)
quiver(X,Y,uval,vval)
