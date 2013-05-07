function b=bs3m(g,r)
% b=bs3m(g);
%
% Given a vector of 24 Gauss coefficients g, bs calculates the spherical
% radial component of the magnetic field and produces a Mollweide
% projection of the field. Specify the (r,theta,phi) locations for
% calculations in the code below.  

% Original bs written 9 April 2008 by Doug Kelley. This version Axl 041912


scalefac = 1; % Should match value used in gc.m!
%r = 20;
%r=1460/1536;  %1460/1536 "core mantle boundary" when r=1 from getcoeff_HCN.m
theta=[0:pi/200:199/200*pi]; % length must match size expected by expmodeproj2
phi=[0:2*pi/200:199/200*2*pi]; % length must match size expected by expmodeproj2
phase=0; % rotates the plot; see expmodeproj2. 

if nargin<1
    error('Usage: b=bs3m(g)');
end
 if length(g)~=24
     error('Argument must be a vector of length 24.');
 end

s=warning('off','MATLAB:divideByZero');
[ttheta,pphi]=meshgrid(theta,phi);
load('Mollweidecoords.mat','xx');
load('Mollweidecoords.mat','yy');
b=zeros(size(ttheta)); % b(theta,phi)
for ii=1:length(g)
    b=b+g(ii)*gmode3m(ii,r,ttheta,pphi);
end
b=b/scalefac;
figure;
pcolor(xx,yy,b);shading flat;
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'Box','off',...
    'Visible','off','Xlim',[-3 3],'YLim',[-1.5 1.5],'YDir','normal');
warning(s.state,'MATLAB:divideByZero');



end 
