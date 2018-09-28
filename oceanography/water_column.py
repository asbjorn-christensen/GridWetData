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

## fit 3 pts to parabola by naive 3x3 fit as y(z) = a*z**2 + b*z + c
#  @param z      : z0,z1,z2
#  @param y      : y0,y1,y2
#  @return a,b,c,problem_solvable
def fit_3pt_to_parabola(z, y):
    assert len(z) == 3 
    assert len(y) == 3 
    det = (-z[1] + z[2])*z[0]**2 + (z[0] - z[2])*z[1]**2 + (-z[0] + z[1])*z[2]**2
    problem_solvable = abs(det)>1e-20
    if problem_solvable: # non singular situation
        a   = (z[2]*(y[0] - y[1]) + z[0]*(y[1] - y[2]) + z[1]*(-y[0] + y[2]))/det    
        b   = ((-y[1] + y[2])*z[0]**2 + (y[0] - y[2])*z[1]**2 + (-y[0] + y[1])*z[2]**2) / det   
        c   = ((z[2]*y[1] - z[1]*y[2])*z[0]**2 + (-z[2]*y[0] + z[0]*y[2])*z[1]**2 +
               (z[1]*y[0] - z[0]*y[1])*z[2]**2) / det   
    else: # singular situation
        a = 0.
        b = 0.
        c = 0.
    return a,b,c,problem_solvable

def point_with_largest_value(z0,v0,z1,v1):
    if v1>v0:
        return z1,v1
    else:
        return z0,v0
    
## Find the maximum value of a scalar and the depth of the maximum for a set of vertical columns
#  @param zv      shape (nz,) : vertical columns of depths (ascending order) with nz points
#  @param yv      shape (nz,) : vertical columns of scalar with nz values (at points zv)
#  @param ib                  : highest wet point index in column; ib < 0 indicates a dry column
#  @param zstep               : (optional) z step for spline curvature refit
#  @param zdry                : (optional) depth value to apply in output for dry points
#  @param dryval              : (optional) max-gradient value to apply in output for dry points
#  @return zmax, maxval       : depth at maximum yv and value of maximum yv
#  
#  fit cubic spline to (z, y). If no interior maxima, select end point with max value
def find_maximum_value(zv, yv, ib, zstep=1.0, zdry=0.0, dryval=0.0):
    #
    # --- branch out depending on number of points ib in water column
    #
    if ib < 0:    # --- dry point
        
        return zdry, dryval
    
    elif ib == 0: # --- single layer situation (assign that value + depth)
        
        return zv[0], yv[0] # assign that depth and value
        
    elif ib == 1: # --- double layer situation (linear case, select max end point)
        
        return point_with_largest_value(zv[0],yv[0],zv[1],yv[1])
        
    elif ib == 2: # --- parabolic situation - max at bound or top
        
        a,b,c,problem_solvable = fit_3pt_to_parabola(zv[:3], yv[:3])
        # non singular situation + concave + interior max (rely on left-to-right evaluation)
        if problem_solvable and a < -1e-20 and (zv[0] < -b/2/a < zv[2]):
            ztop = -b/2/a
            return ztop, a*ztop**2 + b*ztop + c
        else: # all other cases, max at end point
            return point_with_largest_value(zv[0],yv[0],zv[2],yv[2])
            
    else: # --- apply cubic spline for more than 3 data points (most cases should end here)
        
        zv = zv[:(ib+1)] # waste dry parts of water column
        yv = yv[:(ib+1)] # waste dry parts of water column
        spli       = UnivariateSpline(zv, yv, s=0)  # s=0: no smoothing
        dspli_dz   = spli.derivative(1)      
        d2spli_dz2 = spli.derivative(2)
        # currently finding roots unsupported for non-cubic splines, so refit dspli_dz
        zsampl     = arange(zv[0], zv[-1], zstep)
        if len(zsampl)<4:
            zsampl = linspace(zv[0], zv[-1], 4) # sampling distance < zstep
        drefit = UnivariateSpline(zsampl, dspli_dz(zsampl), s=0) # required to apply root()
        roots      = drefit.roots()   # max/min
        if len(roots)>0: # interior max/min values
            best_maxval = -1e20
            z_at_maxval = None
            for z in roots:
                value = spli(z)
                if d2spli_dz2(z)<0 and value > best_maxval:
                    best_maxval = value
                    z_at_maxval = z
            if z_at_maxval is None: # no maxima in roots of dspli_dz
                return point_with_largest_value(zv[0],yv[0],zv[-1],yv[-1])
            else:
                return z_at_maxval, best_maxval # z for best max and best max
        else:  # max/min at end points
            return  point_with_largest_value(zv[0],yv[0],zv[-1],yv[-1])
            


## Find the highest gradient of a scalar and the depth of the maximum gradient for a set of vertical columns
#  @param zv      shape (nz,) : vertical columns of depths (ascending order) with nz points
#  @param yv      shape (nz,) : vertical columns of scalar with nz values (at points zv)
#  @param ib                  : highest wet point index in column; ib < 0 indicates a dry column
#  @param zstep                 : (optional) z step for spline curvature refit
#  @param zdry                  : (optional) depth value to apply in output for dry points
#  @param graddry               : (optional) max-gradient value to apply in output for dry points
#  @return zmax, gradmax        : depth of the maximum gradient and value of maximum gradient
#
#  fit cubic spline to (z, y). Refit cubic spline to second derivative and detect max gradients.
#  of no interior max gradients, select end point with max gradient; linear limit: select mid point
def find_highest_gradient(zv, yv, ib, zstep=1.0, zdry=0.0, graddry=0.0):
    #
    # --- branch out depending on number of points ib in water column
    #
    if ib < 0:    # --- dry point

        return zdry, graddry
        
    elif ib == 0: # --- single layer situation (no gradient assigned)

        return zv[0], 0.0 # assign gradient to layer position and only sensible value
       
    elif ib == 1: # --- double layer situation (set mid point)

        return 0.5*(zv[0]+zv[1]), (yv[1]-yv[0])/(zv[1]-zv[0]) # midpoint and linear interpolation
        
    elif ib == 2: # --- parabolic situation - max grad at upper/lower bound
        a,b,c,problem_solvable = fit_3pt_to_parabola(zv[:3], yv[:3])
        if problem_solvable: # non singular situation
            dydz0 = 2*a*zv[0] + b
            dydz2 = 2*a*zv[2] + b
            if dydz0>dydz2:
                return zv[0], dydz0
            else:
                return zv[2], dydz2
        else: # singular situation -> apply linear limit
            return 0.5*(zv[0]+zv[2]), (yv[2]-yv[0])/(zv[2]-zv[0]) # midpoint and linear interpolation
            
    else: # --- apply cubic spline for more than 3 data points (most cases should end here)
        
        zv = zv[:(ib+1)] # waste dry parts of water column
        yv = yv[:(ib+1)] # waste dry parts of water column
        spli       = UnivariateSpline(zv, yv, s=0)  # s=0: no smoothing
        dspli_dz   = spli.derivative(1)      
        d2spli_dz2 = spli.derivative(2)
        zsampl     = arange(zv[0], zv[-1], zstep)
        if len(zsampl)<4:
            zsampl = linspace(zv[0], zv[-1], 4) # sampling distance < zstep
        # currently finding roots unsupported for non-cubic splines, so refit dspli_dz
        d2refit = UnivariateSpline(zsampl, d2spli_dz2(zsampl), s=0) # required to apply root()
        roots = d2refit.roots()   # max/min gradients
        if len(roots)>0: # interior max/min gradients
            dydz = dspli_dz(roots)
            imax = argmax(dydz)
            return roots[imax], dydz[imax]
        else:  # max/min gradients at end points
            zendpt = [zv[0], zv[-1]]
            dydz   =  dspli_dz( zendpt )
            imax = argmax(dydz)
            return zendpt[imax], dydz[imax]



## Find the lowest gradient of a scalar and the depth of the minimum gradient for a set of vertical columns
#  @param zv      shape (nz,) : vertical columns of depths (ascending order) with nz points
#  @param yv      shape (nz,) : vertical columns of scalar with nz values (at points zv)
#  @param ib                  : highest wet point index in column; ib < 0 indicates a dry column
#  @param zstep               : (optional) z step for spline curvature refit
#  @param zdry                : (optional) depth value to apply in output for dry points
#  @param graddry             : (optional) max-gradient value to apply in output for dry points
#  @return zmin, gradmin      : depth at minimum gradient and value of minimum gradient of yv
#
#  fit cubic spline to (z, y). Refit cubic spline to second derivative and detect max gradients.
#  of no interior max gradients, select end point with max gradient; linear limit: select mid point
def find_lowest_gradient(zv, yv, ib, zstep=1.0, zdry=0.0, graddry=0.0):
    #
    # --- branch out depending on number of points ib in water column
    #
    if ib < 0:    # --- dry point

        return zdry, graddry
        
    elif ib == 0: # --- single layer situation (no gradient assigned)

        return zv[0], 0.0  # assign gradient to layer position and only sensible value
        
    elif ib == 1: # --- double layer situation (set mid point)

        return 0.5*(zv[0]+zv[1]), (yv[1]-yv[0])/(zv[1]-zv[0]) # midpoint and linear interpolation
       
    elif ib == 2: # --- parabolic situation - min grad at upper/lower bound
        
        a,b,c,problem_solvable = fit_3pt_to_parabola(zv[:3], yv[:3])
        if problem_solvable: # non singular situation
            dydz0 = 2*a*zv[0] + b
            dydz2 = 2*a*zv[2] + b
            if dydz0<dydz2:  # NB: opposite find_highest_gradient
                return zv[0], dydz0
            else:
                return zv[2], dydz2
        else: # singular situation -> apply linear limit
            return 0.5*(zv[0]+zv[2]), (yv[2]-yv[0])/(zv[2]-zv[0]) # midpoint and linear interpolation
            
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
            return roots[imin], dydz[imin] 
        else:  # max/min gradients at end points
            zendpt = [zv[0], zv[-1]]
            dydz   =  dspli_dz( zendpt )
            imin = argmin(dydz)  # NB: opposite find_highest_gradient
            return zendpt[imin], dydz[imin]



## Find vertical integral for a vertical column
#  @param zv      shape (nz,) : vertical column of depths (ascending order) with nz points
#  @param yv      shape (nz,) : vertical column of scalar with nz values (at points zv)
#  @param cwv                 : cell widths of vertical column (from surface and down)
#  @param ib                  : highest wet point index in column; ib < 0 indicates a dry column
#
#  trapez summation, with constant extrapolation at end points (half of upper cell + lowest wet cell)
def find_vertical_integral(zv, yv, cwv, ib):
    if ib < 0:    # --- dry point
        return 0.0
    elif ib == 0: # --- single layer situation: value * cell width
        return cwv[0]*yv[0] # layer width * value
    else:
        trapez = 0.5*( yv[1:(ib+1)] + yv[0:ib] )*( zv[1:(ib+1)] - zv[0:ib] ) # trapez sum parts + extrapolation
        intg   = 0.5*cwv[0]*yv[0] + sum(trapez) + 0.5*cwv[ib]*yv[ib] # trapez sum + tails
        return intg
    

## Evaluate water density for a set of npt vertical columns top-down (single column version)
#  @param z       shape (mz,)     : layer center positions [m] of vertical column 
#  @param s       shape (mz,)     : salinity [PSU] at z positions
#  @param t       shape (mz,)     : temperature [Celcius]  at z positions
#  @param patm    (optional) patm : atmospheric pressure [bar]
#  @param maxiter (optional)      : iteration cap on pressure loop
#  @param dpmax   (optional)      : max change in pressure [bar] in pressure loop
#  @return rhow   shape (mz,)     : water density [kg/m3] corresponding to input at positions z
#
#  iterate column pressure to self consistency
#  Water density and pressure are determined by two coupled, non-linear relations: 
#          1) dp   = g rhow(p) dz
#          2) rhow = rho_UNESCO(s,T,p)
#  (rhow, p=pressure) are cell-centered and considered cell-wise constants and therefore
#  eqs. 1,2 reads (iz = 0 .. mz-1)
#          3)  p[iz]  = p[iz-1] + 0.5*g * dz(iz-1) * (rhow(iz)-rhow(iz-1))
#          4)  rhow   = rho_UNESCO(s,T,p)
#  3,4 are iterated to selfconsistency for each layer, after which it is solved 
#  for next layer below, iz+1. The iterative loop usually converges to numerical 
#  precision within 3-4 iteratative cycles of eqs 3,4
def evaluate_water_density(z, s, t, patm = 1.01325, maxiter=10, dpmax=1e-20):
      gfac   = 9.81/1.e5  # g * conversion from Pa to bar
      mz    = len(z)
      
      dz    = z[1:] - z[:-1]  # spacing between cell centers in each column, counting from 1st-2nd layer
      dz    = where(dz>0, dz, 0)  # in case z is not an increasing function for dry layers (it should not be decreasing)
      s     = where(s>0,       s,   1e-10)  # avoid numerical exceptions at dry points at vector evaluation
      t     = where(t>-273.16, t, -273.16)  # avoid numerical exceptions at dry points at vector evaluation
      p     = patm*ones(len(z), float)  # pressure at points z
      p_old = 1.0*p                     # convergence monitoring
      for it in range(maxiter):
            rhow   = rho_UNESCO(s,t,p)
            p[0] = patm + gfac*z[0]*rhow[0] # at center point of first layer
            for iz in range(1,mz):
                  p[iz]  = p[iz-1] + 0.5 * gfac * dz[iz-1] * (rhow[iz]-rhow[iz-1])
            # check convergence       
            dp = amax(abs(p-p_old))
            #print it, dp
            if dp < dpmax:
                  break
            else:
                  p_old = 1.0*p # force copy
      return rhow
            

  
## Evaluate water density for a set of npt vertical columns top-down (vectorized over columns)
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
def evaluate_water_density_vectorized(z, s, t, patm = 1.01325, maxiter=10, dpmax=1e-20):
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
if __name__ == "__main__":
    z = linspace(0, 40, 100)
    #s = linspace(30, 35, 100)
    #t = linspace(15, 5, 100)
    #y = 1.0/(1+exp(-0.3*(z-17)))
    #print find_max_vertical_gradient([z], [y], [len(z)-1])
    #for zz,yy in zip(z,y):
    #    print zz,yy 
    #print find_lowest_gradient([z], [y], [len(z)-1])
    #
    #rhow = evaluate_water_density(array([z]), array([s]), array([t]))
    #for (zz,dens) in zip(z,rhow[0]):
    #    print zz,dens
    #y = z*z*exp(-z)
    acc = arange(100)**1.4
    z  = 0.5*(acc[1:]+acc[:-1])
    cw = acc[1:]-acc[:-1]
    y  = 3 + 4*z
    #print find_maximum_value([z], [y], [len(z)-1])
    print find_vertical_integrals(z, y, cw, len(z)-1)
    zm = acc[-1]
    print 3*zm + 2*zm**2
