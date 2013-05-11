function b=bs3msphere(g,npoints,varargin)
% b=bs3msphere(g);
%
% Given a vector of 24 Gauss coefficients g, bs calculates the spherical
% radial component of the magnetic field and produces a 3D sphere plot
%of the field. Specify the (r,theta,phi) locations for
% calculations in the code below.  

% Original bs written 9 April 2008 by Doug Kelley. This version Axl 041912


scalefac = 1; % Should match value used in gc.m!
%r = 20;
%r=1460/1536;  %1460/1536 "core mantle boundary" when r=1 from getcoeff_HCN.m
r=1
if nargin>2
    r = varargin{1};
end

[xs ys zs] = sphere(npoints);
xs = r*xs; ys = r*ys; zs = r*zs;
phi = atan2(ys,xs);
theta = acos(zs./sqrt(xs.^2+ys.^2+zs.^2));

phase=0; % rotates the plot; see expmodeproj2. 

if nargin<1
    error('Usage: b=bs3m(g)');
end
 if length(g)~=24
     error('Argument must be a vector of length 24.');
 end

s=warning('off','MATLAB:divideByZero');
b=zeros(size(theta)); % b(theta,phi)
for ii=1:length(g)
    b=b+g(ii)*gmode3m(ii,r,theta,phi);
end
b=b/scalefac;
figure;
surf(xs,ys,zs,b);shading flat;
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'Box','off',...
    'Visible','on','Xlim',[-1.5 1.5],'YLim',[-1.5 1.5],'YDir','normal');
axis equal;
axis off;
colormap redblue;
warning(s.state,'MATLAB:divideByZero');



end 
