MATLAB codes for analyzing and drawing magnetic fields  

Poloidal vector spherical harmonics used in the codes below are definitions from  
Bullard and Gellman (1954) dx.doi.org/10.1098/rsta.1954.0018 Eq. 13  
Spherical harmonic definitions use Schmidt semi-normalized associated Legendre functions  
from MATLAB's built in legendre() multiplied by (-1)^m to remove Condon-Shortley phase.  

=======================================
Specific to UMD 3m Geodynamo Experiment
=======================================  

gcoeff3m.m (formerly getcoeffaxl.m) : gets gauss coefficients from 3m experiment radial hall probe array time series, requires positions of 31 probes  
probepos.m : outputs probe positions for gcoeff3m.m  
probecalcbias.mat : radial applied fields (in Gauss per ampere of magnet current) calculated at probe positions for debiasing based on magnet current  
gaussorder.txt : order of gauss coefficients in output time series of gcoeff3m.m  
eqw_coilcoords.mat : xc,yc,zc cartesian coordinates of the white equatorial coil for early 3m experiments with applied field; used in coilcalc.m  
  
===============================================
General Vector Spherical Harmonic Drawing Codes
===============================================
  
bline_spag.m : fast vectorized code for calculating external field lines from Gauss coefficients and line start positions  
gm3cart_gauss.m : required by bline_spag, sums gm3cart.m calculations over all vector spherical harmonics weighted by gauss coefficients  
gm3cart.m : calculates Cartesian magnetic field components of a single vector spherical harmonic  
gmode3m.m : calculates just the spherical radial component, used in surface plots  
gm3comp.m : calculates three spherical components of a single harmonic directly from spherical coordinates  
bs3m.m : plots a Mollweide projection at a chosen radius from a vector of gauss coefficients  
bs3msphere.m : plots a spherical surface plot from gauss coefficient vector  
Mollweidecoords.mat : x,y Mollweide coordinates for bs3m.m  

================================
Additional Useful Magnetic Codes  
================================

bcart2bsph.m : takes cartesian locations x,y,z and cartesian vector field components Bx, By, Bz and outputs spherical coords & components  
bcart2bsph_linemasks.m : same as bcart2bsph.m but also adds masks for spherical radial field >0 and <0 used to color lines by polarity  
coilcalc.m : Biot-Savart finite element code, field at point x,y,z given a coil represented as a parametric curve xc,yc,zc and a current  
selrules.m : calculates Bullard and Gellman 1954 selection rules and outputs TeX/MathJAX formatted lists of allowed couplings.


