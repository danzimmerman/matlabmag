function J = redblue(m,whiteloc)
% Usage: J=redblue([m],[whiteloc])
% redblue: the red-white-blue colormap that Nick likes! Optionally specify
% the colormap size in "m" and the normalized location of the white point 
% in "whiteloc", with 0 being the minimum and 1 being the maximum. See also
% redblackblue.m. 

% Written 17 July 2009 by Doug Kelley.
% Updated 21 April 2010 to allow shifting the white point.
% Fixed to allow whiteloc=0 and whiteloc=1 24 March 2011. 

whitelocdefault=0.5;

if ~exist('m','var') || isempty(m)
   m = size(get(gcf,'colormap'),1);
end
if ~exist('whiteloc','var') || isempty(whiteloc)
    whiteloc=whitelocdefault;
end

if floor(whiteloc*m)==0
    J=[ ones(m,1) linspace(1,0,m)' linspace(1,0,m)' ];
elseif floor(whiteloc*m)==m
    J=[ linspace(0,1,m)' linspace(0,1,m)' ones(m,1)];
else
    J=[ [linspace(0,1,floor(whiteloc*m))'; ones(m-floor(whiteloc*m),1)]...
        [linspace(0,1,floor(whiteloc*m))'; linspace(1,0,m-floor(whiteloc*m))']...
        [ones(floor(whiteloc*m),1); linspace(1,0,m-floor(whiteloc*m))'] ];
end

