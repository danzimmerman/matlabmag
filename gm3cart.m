%gm3cart generates vector spherical harmonic fields using cartesian components for field line drawings
%change history
%============
% modified 050213 4:57pm by Axl to use different recursion formula to simplify indexing
%http://upload.wikimedia.org/math/0/f/e/0fe861fc364a17a41be2ae89c14da71d.png
%
%$(x^2-1)dP_l^m/dx = -(l+m)(l-m+1)\sqrt{1-x^2}P_l^{m-1}(x)-mxP_l^m(x)$
%
%also seems that Schmidt normalization or 'norm' normalization are UNACCEPTABLE for recursion relations
%this version seems okay but then needs UN-NORMALIZED gauss coefficients :(
%add back in the Schmidt coefficients by hand to match with getcoeffaxl, etc
%===========
%modified 050313 9:00am by Axl to fix sqrt(2) issue with "by hand" Schmidt normalization.
%if m=0 then Plm does NOT get multiplied by sqrt(2*(l-m)!/(l+m)!) in MATLAB
%now this agrees with the results from, say, bs3m.m
%============
%modified 050313 11:46am by Axl to shape vectors or matrices so that legendre() fits shape of other things
%this is to help vectorize the overall code & allow 1D,2D or 3D input arrays of points
%============


function [Sx Sy Sz] = gm3cart(n,x,y,z);
if ~(nargin==4)
    error('Usage: [Sx Sy Sz] = gm3cart(n,x,y,z)')
end

if (isvector(x) & (~(isrow(x) & isrow(y) & isrow(z))))
	x = x';
	y = y';
	z = z';
end


r = sqrt(x.^2+y.^2+z.^2);
phi = atan2(y,x);
theta = acos(z./r);
rs = size(r);
phs = size(phi);
ths = size(theta);

switch n %all in order

case 1 
l=1; m=0; %g10
phifun = ones(size(phi)); %cos(m*phi) for m=0
dphifun = zeros(size(phi)); 

case 2
l=1; m=1; %g11c
phifun = cos(m*phi);
dphifun = -m*sin(m*phi); %d_\phi cos(m*phi) for phi component 

case 3
l=1; m=1; %g11s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 4
l=2; m=0; %g20
phifun = 1;
dphifun = 0;

case 5
l=2; m=1; %g21c
phifun = cos(m*phi);
dphifun = -m*sin(m*phi);

case 6
l=2; m=1; %g21s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 7
l=2; m=2; %g22c
phifun = cos(m*phi);
dphifun = -m*sin(m*phi);

case 8
l=2; m=2; %g22s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 9
l=3; m=0; %g30
phifun = 1;
dphifun = 0;

case 10
l=3; m=1; %g31c
phifun = cos(m*phi);
dphifun = -m*sin(m*phi);

case 11
l=3; m=1; %g31s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 12
l=3; m=2; %g32c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 13
l=3; m=2; %g32s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 14
l=3; m=3; %g33c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 15
l=3; m=3; %g33s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 16
l=4; m=0; %g40
phifun = 1;
dphifun = 0;

case 17
l=4; m=1; %g41c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 18
l=4; m=1; %g41s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 19
l=4; m=2; %g42c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 20
l=4; m=2; %g42s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 21
l=4; m=3; %g43c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 22
l=4; m=3; %g43s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 23
l=4; m=4; %g44c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 24
l=4; m=4; %g44s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 25
l=5; m=0; %g50s
phifun = 1;
dphifun = 0;

case 26
l=5; m=1; %g51c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 27
l=5; m=1; %g51s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 28
l=5; m=2; %g52c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 29
l=5; m=2; %g52s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 30
l=5; m=3; %g53c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 31
l=5; m=3; %g53s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 32
l=5; m=4; %g54c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 33
l=5; m=4; %g54s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 34
l=5; m=5; %g55c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 35
l=5; m=5; %g55s
phifun = sin(m*phi)

case 36
l=6; m=0; %g60s
phifun = 1;
dphifun = 0;

case 37
l=6; m=1; %g61c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 38
l=6; m=1; %g61s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 39
l=6; m=2; %g62c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 40
l=6; m=2; %g62s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 41
l=6; m=3; %g63c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 42
l=6; m=3; %g63s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 43
l=6; m=4; %g64c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 44
l=6; m=4; %g64s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 45
l=6; m=5; %g65c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 46
l=6; m=5; %g65s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);

case 47
l=6; m=6; %g66c
phifun = cos(m*phi);
dphifun =-m*sin(m*phi);

case 48
l=6; m=6; %g66s
phifun = sin(m*phi);
dphifun = m*cos(m*phi);



end %of switch/case

legmatrix = legendre(l,cos(theta)); %this is all we need now

Plm = squeeze(legmatrix(m+1,:,:,:)); %m+1th row desired

%for theta derivative we now need P_l^{m-1} which for m=0 requires us to get P_l^{-1}
%using http://en.citizendium.org/images/math/b/f/9/bf9b81c08165403abea9eda69fbcfd3d.png
%this is related to P_l^{-m} =  (-1)^|m| (l-|m|)! /(l+|m|)! P_l^m or  -(l-1)!/(l+1)! * P_l^1
if (m==0)
	Plmm1 = -factorial(l-1)/factorial(l+1) * squeeze(legmatrix(2,:,:,:)); % second row is |m| = 1 
else
	Plmm1 = squeeze(legmatrix(m,:,:,:)); %the mth row is then P_l^{m-1} which is okay for all m>0
end


if ~(m==0)
	Plm = Plm*sqrt(2*factorial(l-m)/factorial(l+m)); %make them Schmidt semi-normalized by HAND
	Plmm1 = Plmm1*sqrt(2*factorial(l-m)/factorial(l+m)); %then this should be okay!
end


%rsize = size(r)
%Plmsize = size(Plm)
%phisize = size(phifun)
Sr = (l*(l+1)./(r.^(l+2))).*Plm.*phifun; % radial component
Stheta = (-l./(r.^(l+2))).*(((l+m)*(l-m+1)*sin(theta).*Plmm1)+(m*cos(theta).*Plm))./(-sin(theta)).*phifun; %theta component
Sphi = (-l./(sin(theta).*r.^(l+2))).*Plm.*dphifun; %phi component

Sx = Sr.*sin(theta).*cos(phi)+Stheta.*cos(theta).*cos(phi)-Sphi.*sin(phi);
Sy = Sr.*sin(theta).*sin(phi)+Stheta.*cos(theta).*sin(phi)+Sphi.*cos(phi);
Sz = Sr.*cos(theta)-Stheta.*sin(theta);

end %of function
