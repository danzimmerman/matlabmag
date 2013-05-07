%gm3cart_gauss sums vector spherical harmonic fields from gm3cart.m in cartesian components for field line drawings

function [Sx Sy Sz] = gm3cart_gauss(g,x,y,z);
if ~(nargin==4)
    error('Usage: [Sx Sy Sz] = gm3cart_gauss(g,x,y,z)')
end
[Sx1 Sy1 Sz1] = gm3cart(1,x,y,z); %initialize matrices with l=1 m=0 
Sx = g(1)*Sx1;
Sy = g(1)*Sy1;
Sz = g(1)*Sz1;
NCOEFF = length(g);
for gn = 2:NCOEFF %do the rest
	[Sxn Syn Szn] = gm3cart(gn,x,y,z);
	Sx = Sx + g(gn)*Sxn;
	Sy = Sy + g(gn)*Syn;
	Sz = Sz + g(gn)*Szn;
end
