%notes: realized that initial point cloud is drifting at \omega when it should really drift at \omega/m
%here's a new top view glatzmaier colors 
set(0,'defaultfigurevisible','off');
load /data/3m/082212/movies/fulltqswitchmovdata.mat
KT = find(tq>53650 & tq<53700); %times between 53650-53700
gf = gb(KT,:); %gb is bandpass filtered 0.1-6Hz
tq = tq(KT,:);
SAVEDATFILES = 1; %set to 1 for saving, 0 for not
DSFAC = 10;
NSPHERE = 6;
gd = downsample(gf,DSFAC);
td = downsample(tq,DSFAC);
DIRNAME = '/data/3m/082212/movies/spag2';
AXLIMTOP=5;
AXLIMSIDE = 3;
STEPSIZE = 0.01; %field line integration step size
NEGINTFAC = 1; %factor for inward integration
NSTART = 1; %if you need to pause and restart
if exist(DIRNAME)
	cd(DIRNAME)
	system(['cp /home/axl/matlabgit/mag/spagmovscript.m ' DIRNAME '/thisinstance_spagmovscript.m'])

else
	mkdir(DIRNAME)
	cd(DIRNAME)
	system(['cp /home/axl/matlabgit/mag/spagmovscript.m ' DIRNAME '/thisinstance_spagmovscript.m'])
end

if ~(exist('top'))
	mkdir('top');
end
if ~(exist('side'))
	mkdir('side');
end
if ~(exist('dave'))
    mkdir('dave');
end

[x0 y0 z0] = sphere(NSPHERE);
WAVEFREQ = 0.04; %rough guess from movie
m = 1;
DRIFTDIR = -1; %(CCW from top bs3msphere)
NUMTIMESTEP = length(td);
td0 = td-td(1);
thetat = repmat(DRIFTDIR*2*pi*WAVEFREQ*td0/m,[1 size(x0)]);

x0t = ones([length(td) size(x0)]);
y0t = ones([length(td) size(x0)]);
z0t = ones([length(td) size(x0)]);

for j = 1:length(td)
x0t(j,:,:) = x0.*squeeze(cos(thetat(j,:,:)))-y0.*squeeze(sin(thetat(j,:,:)));
y0t(j,:,:) = x0.*squeeze(sin(thetat(j,:,:)))+y0.*squeeze(cos(thetat(j,:,:)));
z0t(j,:,:) = z0;
end
ttotal = 0

for j = NSTART:length(td)
	tic
	[b figj] = bs3msphere_glatzcol(gd(j,:),200); %"glatzmeier" style color map
	hold on;
	x0j = squeeze(x0t(j,:,:));
	y0j = squeeze(y0t(j,:,:));
	z0j = squeeze(z0t(j,:,:));
	[xt yt zt Bxt Byt Bzt xll yll zll] = bline_spag(gd(j,:),3*x0j(:),3*y0j(:),3*z0j(:),STEPSIZE*NEGINTFAC,5000,-1,1,5);
	rll = sqrt(xll.^2+yll.^2+zll.^2);
	rlmin = min(rll(~isnan(rll)))
	%rll = rll/rlmin;
	xll = 1.01*xll(1<rll<1.05)./rll;
	yll = 1.01*yll(1<rll<1.05)./rll;
	zll = 1.01*zll(1<rll<1.05)./rll;
	[xt yt zt Bxt Byt Bzt xlast ylast zlast] = bline_spag(gd(j,:),xll(:),yll(:),zll(:),STEPSIZE,5000,1,1,5);
	[r theta phi Br Btheta Bphi Bmag maskp maskn] = bcart2bsph_linemasks(xt,yt,zt,Bxt,Byt,Bzt);
	hold on
	plot3((maskp.*xt),(maskp.*squeeze(yt)),(maskp.*squeeze(zt)),'color',[0.8 0.66 0],'linewidth',2);
	plot3((maskn.*squeeze(xt)),(maskn.*squeeze(yt)),(maskn.*squeeze(zt)),'color',[0 0.7 0.7],'linewidth',2);
	plot3([0 0],[0 0],[-1.5 1.5],'w','linewidth',3);
	plot3(0,0,1.5,'^w','markersize',6,'markerfacecolor','w');
	
	xlim([-AXLIMSIDE AXLIMSIDE]); ylim([-AXLIMSIDE AXLIMSIDE]); zlim([-AXLIMSIDE AXLIMSIDE]);
	set(figj,'InvertHardCopy','off','color','k');
	fn = sprintf([DIRNAME '/side/%05d.png'],j);
	saveas(figj,fn);
	set(gca,'view',[0 90]); %top view!
	xlim([-AXLIMTOP AXLIMTOP]); ylim([-AXLIMTOP AXLIMTOP]); zlim([-AXLIMTOP AXLIMTOP]);
	fn = sprintf([DIRNAME '/top/%05d.png'],j);
	saveas(figj,fn);
	if SAVEDATFILES
		fndat = sprintf([DIRNAME '/dave/%05d.dat'],j);
		unwind(xt,yt,zt,fndat);
		close
	    figure; pcolor(b); 
	    shading flat
	    colormap glatzmap
	    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
	    fnsurf = sprintf([DIRNAME '/dave/%05d.png'],j);
	    saveas(gcf,fnsurf);
	    close   
    end
    
	tj = toc
	ttotal = ttotal +tj
	meant = ttotal/j
end
set(0,'defaultfigurevisible','on');
