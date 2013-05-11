function J = glatzmap(m,blackloc)
% Usage: J=glatzmap([m],[blackloc])
% glatzmap, colormap similar to Glatzmaier et al dynamo simulations. Optionally specify the colormap 
% size in "m" and the normalized location of the black point in "blackloc", 
% with 0 being the minimum and 1 being the maximum. See also redblue.m. 

% Modified by Axl from the "redblue" map by Doug Kelley


blacklocdefault=0.5;

if ~exist('m','var') || isempty(m)
   m = size(get(gcf,'colormap'),1);
end
if ~exist('blackloc','var') || isempty(blackloc)
    blackloc=blacklocdefault;
end

if floor(blackloc*m)==0
    J=[ linspace(0,1,m)' zeros(m,1) zeros(m,1) ];
elseif floor(blackloc*m)==m
    J=[ zeros(m,1) zeros(m,1) linspace(1,0,m)'];
else
    J=[ 1*[zeros(floor(blackloc*m),1);linspace(0,1,m-floor(blackloc*m))']... % red channel
        0.82*[zeros(floor(blackloc*m),1);linspace(0,1,m-floor(blackloc*m))']+1*[linspace(1,0,floor(blackloc*m))';zeros(m-floor(blackloc*m),1)]  ... % green channel
        1*[linspace(1,0,floor(blackloc*m))';zeros(m-floor(blackloc*m),1)] ]; % blue channel
end

