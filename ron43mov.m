
load moviedata.mat

[x0 y0 z0] = sphere(10);
WAVEFREQ = 0.64666667;
DRIFTDIR = 1; %(CCW from top bs3msphere)
NUMTIMESTEP = length(td);
td0 = td-td(1);
thetat = repmat(2*pi*0.646666667*td0,[1 size(x0)]);

x0t = ones([length(td) size(x0)]);
y0t = ones([length(td) size(x0)]);
z0t = ones([length(td) size(x0)]);

for j = 1:length(td)
x0t(j,:,:) = x0.*squeeze(cos(thetat(j,:,:)))-y0.*squeeze(sin(thetat(j,:,:)));
y0t(j,:,:) = x0.*squeeze(sin(thetat(j,:,:)))+y0.*squeeze(cos(thetat(j,:,:)));
z0t(j,:,:) = z0;
end

for j = 1:length(td)
	bs3msphere(gd(j,:),200);
	colormap redblue; hold on;
	x0j = squeeze(x0t(j,:,:));
	y0j = squeeze(y0t(j,:,:));
	z0j = squeeze(z0t(j,:,:));
	[xt yt zt Bxt Byt Bzt xll yll zll] = bline_spag(gd(j,:),4*x0j(:),4*y0j(:),4*z0j(:),0.2,5000,-1,1,5);
	figure; plot3(xt,yt,zt);
	rll = sqrt(xll.^2+yll.^2+zll.^2);
	xll = xll(0.99<rll<1.01);
	yll = yll(0.99<rll<1.01);
	zll = zll(0.99<rll<1.01);
	zls = size(zll)
	[xt yt zt Bxt Byt Bzt xlast ylast zlast] = bline_spag(gd(j,:),squeeze(xll),squeeze(yll),squeeze(zll),0.05,5000,1,1,5);
	[r theta phi Br Btheta Bphi Bmag maskp maskn] = bcart2bsph_linemasks(squeeze(xt),squeeze(yt),squeeze(zt),squeeze(Bxt),squeeze(Byt),squeeze(Bzt));
	hold on
	xts = size(xt)
	plot3((maskp.*squeeze(xt)),(maskp.*squeeze(yt)),(maskp.*squeeze(zt)),'color',[0.8 0 0],'linewidth',2);
	plot3((maskn.*squeeze(xt)),(maskn.*squeeze(yt)),(maskn.*squeeze(zt)),'color',[0 0 0.8],'linewidth',2);
	xlim([-3 3]); ylim([-3 3]); zlim([-3 3]);
	fn = sprintf('spag2/%05d.png',j);
	saveas(gcf,fn);
	close
end

