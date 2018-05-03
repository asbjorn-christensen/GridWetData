#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package netcdf_aux
#        Auxillary functions for convenient netcdf I/O
#        TODO: migrate to GSW package for wider usage of seawater physics
#######################################################################################################################
import os
import exceptions
from numpy import *
import numbers   # python 2.6+
from scipy.interpolate import * 

_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info




## Find the highest gradient of a scalar and the depth of the maximum gradient for a set of vertical columns
#  @param z      shape (npt,nz) : npt vertical columns of depths (ascending order), each with nz points corresponding to y
#  @param y      shape (npt,nz) : npt vertical columns of scalar, each with nz points
#  @param blay   shape (npt,)   : highest wet point index in column; blay < 0 indicates a dry column
#  @param zstep                 : (optional) z step for spline curvature refit
#  @param zdry                  : (optional) depth value to apply in output for dry points
#  @param graddry               : (optional) max-gradient value to apply in output for dry points
#
#  fit cubic spline to (z, y). Refit cubic spline to second derivative and detect max gradients.
#  of no interior max gradients, select end point with max gradient; linear limit: select mid point
def find_highest_gradient(z, y, blay, zstep=1.0, zdry=0.0, graddry=0.0):
    gradmax = []
    zmax    = []
    # ---- loop over all vertical columns ----
    for zv,yv,ib in zip(z, y, blay):
        if ib < 0:    # --- dry point
            gradmax.append(graddry)
            zmax.append(   zdry   )
        elif ib == 0: # --- single layer situation (no gradient assigned)
            gradmax.append( 0.0  ) # only sensible value
            zmax.append(   zv[0] ) # assign gradient to layer position
        elif ib == 1: # --- double layer situation (set mid point)
            gradmax.append((yv[1]-yv[0])/(zv[1]-zv[0]) ) # linear interpolation
            zmax.append(    0.5*(zv[0]+zv[1]) )          # midpoint
        elif ib == 2: # --- parabolic situation - max grad at upper/lower bound
            det = (-zv[1] + zv[2])*zv[0]**2 + (zv[0] - zv[2])*zv[1]**2 + (-zv[0] + zv[1])*zv[2]**2
            if abs(det)>1e-20: # non singular situation
                a   = (zv[2]*(yv[0] - yv[1]) + zv[0]*(yv[1] - yv[2]) + zv[1]*(-yv[0] + yv[2]))/det    
                b   = ((-yv[1] + yv[2])*zv[0]**2 + (yv[0] - yv[2])*zv[1]**2 + (-yv[0] + yv[1])*zv[2]**2) / det   
                c   = ((zv[2]*yv[1] - zv[1]*yv[2])*zv[0]**2 + (-zv[2]*yv[0] + zv[0]*yv[2])*zv[1]**2 +
                       (zv[1]*yv[0] - zv[0]*yv[1])*zv[2]**2) / det   
                dydz0 = 2*a*zv[0] + b
                dydz2 = 2*a*zv[2] + b
                if dydz0>dydz2:
                    gradmax.append( dydz0 )
                    zmax.append(    zv[0] )
                else:
                    gradmax.append( dydz2 )
                    zmax.append(    zv[2] )
            else: # singular situation -> apply linear limit
                gradmax.append((yv[2]-yv[0])/(zv[2]-zv[0]) ) # linear interpolation
                zmax.append(    0.5*(zv[0]+zv[2]) )          # midpoint
        else: # --- apply cubic spline for more than 3 data points (most cases should end here)
            zv = zv[:(ib+1)] # waste dry parts of water column
            yv = yv[:(ib+1)] # waste dry parts of water column
            spli       = UnivariateSpline(zv, yv, s=0)  # s=0: no smoothing
            dspli_dz   = spli.derivative(1)      
            d2spli_dz2 = spli.derivative(2)
            zsampl     = arange(zv[0], zv[-1], zstep)
            if len(zsampl)<4:
                zsampl = linspace(zv[0], zv[-1], 4) # sampling distance < zstep
            d2refit = UnivariateSpline(zsampl, d2spli_dz2(zsampl), s=0) # required to apply root()
            roots = d2refit.roots()   # max/min gradients
            if len(roots)>0: # interior max/min gradients
                dydz = dspli_dz(roots)
                imax = argmax(dydz)
                gradmax.append( dydz[imax] )
                zmax.append(    roots[imax] )
            else:  # max/min gradients at end points
                zendpt = [zv[0], zv[-1]]
                dydz   =  dspli_dz( zendpt )
                imax = argmax(dydz)
                gradmax.append( dydz[imax] )
                zmax.append(    zendpt[imax] )
    # for yv,zv,ib
    return array(zmax, float), array(gradmax, float)


## Find the lowest gradient of a scalar and the depth of the minimum gradient for a set of vertical columns
#  @param z      shape (npt,nz) : npt vertical columns of depths (ascending order), each with nz points corresponding to y
#  @param y      shape (npt,nz) : npt vertical columns of scalar, each with nz points
#  @param blay   shape (npt,)   : highest wet point index in column; blay < 0 indicates a dry column
#  @param zstep                 : (optional) z step for spline curvature refit
#  @param zdry                  : (optional) depth value to apply in output for dry points
#  @param graddry               : (optional) max-gradient value to apply in output for dry points
#
#  fit cubic spline to (z, y). Refit cubic spline to second derivative and detect max gradients.
#  of no interior max gradients, select end point with max gradient; linear limit: select mid point
def find_lowest_gradient(z, y, blay, zstep=1.0, zdry=0.0, graddry=0.0):
    gradmin = []
    zmin    = []
    # ---- loop over all vertical columns ----
    for zv,yv,ib in zip(z, y, blay):
        if ib < 0:    # --- dry point
            gradmin.append(graddry)
            zmin.append(   zdry   )
        elif ib == 0: # --- single layer situation (no gradient assigned)
            gradmin.append( 0.0  ) # only sensible value
            zmin.append(   zv[0] ) # assign gradient to layer position
        elif ib == 1: # --- double layer situation (set mid point)
            gradmin.append((yv[1]-yv[0])/(zv[1]-zv[0]) ) # linear interpolation
            zmin.append(    0.5*(zv[0]+zv[1]) )          # midpoint
        elif ib == 2: # --- parabolic situation - min grad at upper/lower bound
            det = (-zv[1] + zv[2])*zv[0]**2 + (zv[0] - zv[2])*zv[1]**2 + (-zv[0] + zv[1])*zv[2]**2
            if abs(det)>1e-20: # non singular situation
                a   = (zv[2]*(yv[0] - yv[1]) + zv[0]*(yv[1] - yv[2]) + zv[1]*(-yv[0] + yv[2]))/det    
                b   = ((-yv[1] + yv[2])*zv[0]**2 + (yv[0] - yv[2])*zv[1]**2 + (-yv[0] + yv[1])*zv[2]**2) / det   
                c   = ((zv[2]*yv[1] - zv[1]*yv[2])*zv[0]**2 + (-zv[2]*yv[0] + zv[0]*yv[2])*zv[1]**2 +
                       (zv[1]*yv[0] - zv[0]*yv[1])*zv[2]**2) / det   
                dydz0 = 2*a*zv[0] + b
                dydz2 = 2*a*zv[2] + b
                if dydz0<dydz2:  # NB: opposite find_highest_gradient
                    gradmin.append( dydz0 )
                    zmin.append(    zv[0] )
                else:
                    gradmin.append( dydz2 )
                    zmin.append(    zv[2] )
            else: # singular situation -> apply linear limit
                gradmin.append((yv[2]-yv[0])/(zv[2]-zv[0]) ) # linear interpolation
                zmin.append(    0.5*(zv[0]+zv[2]) )          # midpoint
        else: # --- apply cubic spline for more than 3 data points (most cases should end here)
            zv = zv[:(ib+1)] # waste dry parts of water column
            yv = yv[:(ib+1)] # waste dry parts of water column
            spli       = UnivariateSpline(zv, yv, s=0)  # s=0: no smoothing
            dspli_dz   = spli.derivative(1)      
            d2spli_dz2 = spli.derivative(2)
            zsampl     = arange(zv[0], zv[-1], zstep)
            if len(zsampl)<4:
                zsampl = linspace(zv[0], zv[-1], 4) # sampling distance < zstep
            d2refit = UnivariateSpline(zsampl, d2spli_dz2(zsampl), s=0) # required to apply root()
            roots = d2refit.roots()   # max/min gradients
            if len(roots)>0: # interior max/min gradients
                dydz = dspli_dz(roots)
                imin = argmin(dydz)  # NB: opposite find_highest_gradient
                gradmin.append( dydz[imin] )
                zmin.append(    roots[imin] )
            else:  # max/min gradients at end points
                zendpt = [zv[0], zv[-1]]
                dydz   =  dspli_dz( zendpt )
                imin = argmin(dydz)  # NB: opposite find_highest_gradient
                gradmin.append( dydz[imin] )
                zmin.append(    zendpt[imin] )
    # for yv,zv,ib
    return array(zmin, float), array(gradmin, float)



## Evaluate water density for a set of npt vertical columns top-down
#  @param z       shape (npt,nz)  : npt vertical columns of depths [m] (ascending order), each with nz points corresponding to s,t
#  @param s       shape (npt,nz)  : npt vertical columns of salinity [PSU], each with nz points
#  @param t       shape (npt,nz)  : npt vertical columns of temperature [Celcius], each with nz points
#  @param patm    (optional) patm : atmospheric pressure [bar]
#  @param maxiter (optional)      : iteration cap on pressure loop
#  @param dpmax   (optional)      : max change in pressure [bar] in pressure loop
#  @return rhow   shape (npt,nz)  : npt vertical columns of water density [kg/m3] corresponding to input at positions z
#
#  vectorize over columns - it can be ignored whether cells are wet/dry, because evaluation is top-down
#  iterate column pressure to self consistency
#  Water density and pressure are determined by two coupled, non-linear relations: 
#          1) dp   = g rhow(p) dz
#          2) rhow = rho_UNESCO(s,T,p)
#  (rhow, p=pressure) are cell-centered and considered cell-wise constants and therefore
#  eqs. 1,2 reads (iz = 0 .. nz-1)
#          3)  p[iz]  = p[iz-1] + 0.5*g * dz(iz-1) * (rhow(iz)-rhow(iz-1))
#          4)  rhow   = rho_UNESCO(s,T,p)
#  3,4 are iterated to selfconsistency for each layer, after which it is solved 
#  for next layer below, iz+1. The iterative loop usually converges to numerical 
#  precision within 3-4 iteratative cycles of eqs 3,4
def evaluate_water_density(z, s, t, patm = 1.01325, maxiter=10, dpmax=1e-20):
      gfac   = 9.81/1.e5  # g * conversion from Pa to bar
      npt,nz = z.shape
      
      dz    = z[:,1:] - z[:,:-1]  # spacing between cell centers in each column, counting from 1st-2nd layer
      dz    = where(dz>0, dz, 0)  # in case z is not an increasing function for dry layers (it should not be decreasing)
      s     = where(s>0,       s,   1e-10)  # avoid numerical exceptions at dry points at vector evaluation
      t     = where(t>-273.16, t, -273.16)  # avoid numerical exceptions at dry points at vector evaluation
      p     = patm*ones(z.shape, float) # pressure at points z
      p_old = 1.0*p                     # convergence monitoring
      for it in range(maxiter):
            rhow   = rho_UNESCO(s,t,p)
            p[:,0] = patm + gfac*z[:,0]*rhow[:,0] # at center point of first layer
            for iz in range(1,nz):
                  p[:,iz]  = p[:,iz-1] + 0.5 * gfac * dz[:,iz-1] * (rhow[:,iz]-rhow[:,iz-1])
            # check convergence       
            dp = amax(abs(p-p_old))
            #print it, dp
            if dp < dpmax:
                  break
            else:
                  p_old = 1.0*p # force copy
      return rhow
            


##  calculates the density of sea water at a given point using the UNESCO equation of state (from Per Berg, per@dmi.dk)
#   transcripted from IBMlib version in fortran to python 
#   @param  ss :  salinity    [PSU]
#   @param  tt : temperature [Celcius]
#   @param  pp : pressure    [bar] 1 bar = 100 000 pascals (Pa) 
#   @return   local density of water at (ss,tt,pp) [kg/m3]
def rho_UNESCO(ss,tt,pp):
      # ----- parameters -----
      a0 = 999.842594                            
      a1 =   6.793952e-2                            
      a2 =  -9.095290e-3                            
      a3 =   1.001685e-4                            
      a4 =  -1.120083e-6                            
      a5 =   6.536332e-9                            
      b0 =  +8.24493e-1                             
      b1 =  -4.0899e-3                              
      b2 =  +7.6438e-5                              
      b3 =  -8.2467e-7                              
      b4 =  +5.3875e-9                              
      c0 =  -5.72466e-3                             
      c1 =  +1.0227e-4                              
      c2 =  -1.6546e-6                              
      d0 =  +4.8314e-4                              
      e0 = 19652.21                              
      e1 =   148.4206                            
      e2 =    -2.327105                          
      e3 =     1.360477e-2                          
      e4 =    -5.155288e-5
      f0 =    54.6746                            
      f1 =    -0.603459                          
      f2 =     1.09987e-2                           
      f3 =    -6.1670e-5                            
      g0 =     7.944e-2                             
      g1 =     1.6483e-2                            
      g2 =    -5.3009e-4                            
      h0 =  3.239908                                 
      h1 =  1.43713e-3                              
      h2 =  1.16092e-4                              
      h3 = -5.77905e-7                              
      i0 =  2.2838e-3                               
      i1 = -1.0981e-5                               
      i2 = -1.6078e-6                               
      j0 =  1.91075e-4                              
      k0 =  8.50935e-5                              
      k1 = -6.12293e-6                              
      k2 =  5.2787e-8                               
      l0 = -9.9348e-7                               
      l1 =  2.0816e-8                               
      l2 =  9.1697e-10
      # ----- variables -----
      s1 = ss
      s2 = ss*ss
      s3 = s2*ss
      s15= sqrt(s3)
      t1 = tt
      t2 = tt*tt
      t3 = t2*tt
      t4 = t2*t2
      t5 = t4*tt
 
      rhonul = a0+a1*t1+a2*t2+a3*t3+a4*t4+a5*t5 + \
        (b0+b1*t1+b2*t2+b3*t3+b4*t4)*s1      + \
        (c0+c1*t1+c2*t2)*s15                 + \
        d0*s2

      sbm    =   e0+e1*t1+e2*t2+e3*t3+e4*t4 + \
        (f0+f1*t1+f2*t2+f3*t3)*s1 + \
        (g0+g1*t1+g2*t2)*s15 + \
        (  h0+h1*t1+h2*t2+h3*t3 + \
            (i0+i1*t1+i2*t2)*s1 + \
        j0*s15                  ) * pp + \
        (  k0+k1*t1+k2*t2 + \
        (l0+l1*t1+l2*t2)*s1     ) * pp * pp

      return rhonul/(1.e0-pp/sbm)

    

# --- test 
z = linspace(0, 40, 100)
s = linspace(30, 35, 100)
t = linspace(15, 5, 100)
#y = 1.0/(1+exp(-0.3*(z-17)))
#print find_max_vertical_gradient([z], [y], [len(z)-1])
#for zz,yy in zip(z,y):
#    print zz,yy 
#print find_lowest_gradient([z], [y], [len(z)-1])
#
#rhow = evaluate_water_density(array([z]), array([s]), array([t]))
#for (zz,dens) in zip(z,rhow[0]):
#    print zz,dens
