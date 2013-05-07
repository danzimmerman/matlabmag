%coilcalc.m
%Biot Savart routine to calculate magnetic field from 3m coil
%2 August 2012 by Axl
%Output Bx By Bz full 3D to take into account deformed magnet if necessary
%and as a check
%X,Y,Z are the coordinates on which to calculate the magnetic field 
%(XC,YC,ZC) is a set of points that describe the magnet geometry
%like a parametric curve
%coil current is I 
%output is in GAUSS if distance is in meters and current I is in amperes
%Updated 14 August: wrong sign of y-component of cross product - didn't show up in initial verification in XZ plane
%check for verification in /data/tech_info/magcalcs

function [Bx By Bz] = coilcalc(X,Y,Z,I,XC,YC,ZC);

%circshift wants column vectors: 
if isrow(XC)
XC = XC';
end
if isrow(YC)
YC = YC';
end
if isrow(ZC)
ZC = ZC';
end

mu0 = 0.0125663706; %permeability of free space in GAUSS*meter/amp

%below all variables are vectores where each element corresponds to
%a point on the parametrically defined magnet curve

%calculate wire segment midpoints: average j+1th and jth coord
xc_pr = (circshift(XC,-1)+XC)/2;
yc_pr = (circshift(YC,-1)+YC)/2;
zc_pr = (circshift(ZC,-1)+ZC)/2;

%calculate displacement vector components
rx = X-xc_pr;
ry = Y-yc_pr;
rz = Z-zc_pr;

%calculate denominator

magr3 = (rx.^2+ry.^2+rz.^2).^(3/2);

%calculate wire end displacement vector dL

dlx = (circshift(XC,-1)-XC);
dly = (circshift(YC,-1)-YC);
dlz = (circshift(ZC,-1)-ZC);

%calculate integrand
%three components of integrand: dl cross r over r3

dBx = (dly.*rz-dlz.*ry)./magr3;
dBy = (dlz.*rx-dlx.*rz)./magr3; %FIXED 8/14/12 had wrong sign of dBy! 
dBz = (dlx.*ry-dly.*rx)./magr3; 

%now sum up the whole curve:

Bx = mu0/(4*pi)*I*sum(dBx);
By = mu0/(4*pi)*I*sum(dBy); 
Bz = mu0/(4*pi)*I*sum(dBz);


