%plots gauss coefficients as l*(l+1)*g_lm / B_0
%param,gauss,and current should be all the same length
%pass desired gauss coefficients as gaussindex for automatic labeling
function fighandle = gaussnormplot(param,gauss,gaussindex,current,LFLAG,NEWPLOT)

if ~(nargin==6)
	error('Usage : fighandle = gaussnormplot(param,gauss,gaussindex,current,LFLAG,NEWPLOTFLAG)')
end
lstruct = load('lvector.mat');
l = lstruct.l;
B0 = 0.529*current;
gauselect=gauss(:,gaussindex)/0.0315;
vec = 1:24;
clist = 0.35 * (1+[sin(5/3*pi*vec/24); sin(5/3*pi*vec/24+2*pi/3); sin(5/3*pi*vec/24+4*pi/3)]);
%choose symbols in a sensible way for structure of data
symbols = char(1,24);
labels = cell(1,24);

symbols([1 4 9 16]) = 'o'; %'o' for axisymmetric
symbols([2 5 10 17]) = 'v'; %'v' for m=1 cosine
symbols([3 6 11 18]) = '^'; %'^' for m=1 sine
symbols([7 12 19]) = '>'; %'>' for m=2 cosine
symbols([8 13 20]) = '<'; %'<' for m=2 sine
symbols([14 21]) = 's'; %square for m=3 cosine
symbols([15 22]) = 'd'; %diamond for m=3 sine
symbols(23) = 'x'; %'x' for m=4 cosine
symbols(24) = '+'; %'+' for m=4 sine
symsize = 8*ones(1,24);
symsize(23:24) = 10;

syminuse = symbols(gaussindex);
if NEWPLOT
	fighandle = figure; subplot(2,1,1); hold on
end
load sphlabels.mat
for j = 1:length(gaussindex);
	if LFLAG
		gaumult = l(j)*(l(j)+1)./B0; %l(l+1)
	else
		gaumult = 1;
	end
	
	plot(param,gaumult.*gauselect(:,j),'.-','linewidth',2,'marker',symbols(gaussindex(j)),'markersize',symsize(gaussindex(j)),'color',clist(:,gaussindex(j)),'markerfacecolor',clist(:,gaussindex(j)));
end
if LFLAG
	legend(blabeltex{gaussindex})
else
	legend(glabeltex{gaussindex})
end
h = findobj('tag','legend')
set(h,'location','northwest')

end %of function 
