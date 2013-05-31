%spherediag.m is a script for drawing a 3D representation of the 3m system with equatorial magnet, probes, torque sensor, ports
%requires probepos.m for probe positions


%approximate dimensions for drawing the sphere with top lid (not exact references!)
%top lid/outer flange radius 0.854m
%inside neck flange radius 0.76m
%lid thickness 7.94cm
%neck flange height 14cm
%sphere radius 1.46m
%top lid surface height from center
%total outer sphere inner surface to top of flange 24.44cm

OUTOUTR = 1.46+0.025; %sphere outer surface radius 1.46m + 1 inch
OUTINR = 1.46;
THTH = 0.5274; %top hole theta 0.5274 radians
MAGOUTR = 2; %magnet outer radius 2m.
MAGINR = 1.8; %magnet inner radius 1.8m
ZMAGB = -0.066;
ZMAGT = 0.066; %+/-6.6cm from equator.

SPHERECOL = -0.7;
SPHEREICOL = -0.9;
PROBECOL = 0.7;
ISCOL = 0;
WIRECOL = -0.4;
SHOWINNER = 1; %SHOW INNER SPHERE IN CUTOUT
NECKBOTTOM = (OUTOUTR*cos(THTH)-0.0254);
NECKTOP =  NECKBOTTOM+0.2444;
FLANGEOUTR = 0.854;
FLANGEINR = 0.76;
FLANGEHEIGHT = 0.14;
FLANGETOP = NECKTOP; %flange top from center
FLANGEBOTTOM = FLANGETOP-0.14;

RSHAFT = 0.084; %shaft radius 8.4cm
LIDTOP = FLANGETOP+.0794;


[xs ys zs] = sphere(300); %basic sphere






Jg = [linspace(0,1,11); linspace(0,1,11); linspace(0,1,11)]';
Jr = [linspace(0,1,11); 0*linspace(0,1,11); 0*linspace(0,1,11)]';
Jb = [0*linspace(0,1,11); 0*linspace(0,1,11); linspace(0,1,11)]';
J = [Jg; Jr; Jb]

%plot outer sphere and inner surface outer sphere


xo = OUTOUTR*xs;
yo = OUTOUTR*ys;
zo = OUTOUTR*zs; %OUTOUTRm  sphere OUTER diameter
%cut a hole in the top of outer sphere
TOPHOLE = 0;
if TOPHOLE
	hole = find(zo>0 & (xo.^2+yo.^2)<0.76.^2);
	xo(hole) = NaN;
	yo(hole) = NaN;
	zo(hole) = NaN;
end

if SHOWINNER
	windowk = find(xo<0 & yo<0 & zo>0);
	phik = atan2(yo(windowk),xo(windowk));
	thetak = acos(zo(windowk)./sqrt(xo(windowk).^2+yo(windowk).^2+zo(windowk).^2));
	thetamax = max(max(thetak));
	thetamin = min(min(thetak));
	thetas = linspace(thetamin-pi/100,thetamax+pi/100,50);
	thetam = [thetas; thetas];
	
	phimax = max(max(phik));
	phimin = min(min(phik));
	phis = linspace(phimin-pi/100,phimax+pi/100,50);
	phim = [phis; phis];
	rwall = [OUTINR*ones(size(phis)); OUTOUTR*ones(size(phis))];
	xw = rwall.*cos(phim);
	yw = rwall.*sin(phim);
	zw = zeros(size(xw));
	zv = rwall.*cos(thetam);
	xvmax = rwall.*sin(thetam).*cos(phimax+pi/100);
	yvmax = rwall.*sin(thetam).*sin(phimax+pi/100);
	xvmin = rwall.*sin(thetam).*cos(phimin-pi/100);
	yvmin = rwall.*sin(thetam).*sin(phimin-pi/100);
	
	WALLCOL = SPHERECOL;
	
	xo(windowk) = NaN;

end

figure; surf(xo,yo,zo,SPHERECOL*ones(size(zo)));
colormap(J); axis equal; axis off;shading flat; 
hold on

shading flat
surf(xo*OUTINR/OUTOUTR,yo*OUTINR/OUTOUTR,zo*OUTINR/OUTOUTR,SPHEREICOL*ones(size(zo)));
if SHOWINNER
	surf(xw,yw,zw,WALLCOL*ones(size(zw)));
	surf(xvmax,yvmax,zv,WALLCOL*ones(size(zv)));
	surf(xvmin,yvmin,zv,WALLCOL*ones(size(zv)));

end
surf(xh,yh,zh,sign(zh)); %to make easy colormapping plot a little sphere that goes from -1 to 1. 
%plot neck
rhedge = 0.76;
thetahole = [0:pi/100:2*pi];
xhedge = 0.76*cos(thetahole);
yhedge = 0.76*sin(thetahole);
ztopedge = FLANGEBOTTOM*ones(size(xhedge));
zbottomedge = (OUTOUTR*cos(THTH)-0.03)*ones(size(xhedge));
Mhx = [xhedge; xhedge];
Mhy = [yhedge; yhedge];
Mhz = [ztopedge; zbottomedge];
surf(Mhx,Mhy,Mhz,SPHERECOL*ones(size(Mhz)));

%plot hall sensors as tiny spheres

ppos = probepos;
Rp = OUTOUTR*ppos(:,1);
thetap = ppos(:,2);
phip = ppos(:,3);
xp = Rp.*sin(thetap).*cos(phip);
yp = Rp.*sin(thetap).*sin(phip);
zp = Rp.*cos(thetap);
[xh yh zh] = sphere(30);
xh = 0.05*xh;
yh = 0.05*yh;
zh = 0.05*zh;
for j = 1:31
	xpr{j} = xp(j) + xh;
	ypr{j} = yp(j) + yh;
	zpr{j} = zp(j) + zh;
end

for j = 1:31
	surf(xpr{j},ypr{j},zpr{j},PROBECOL*ones(size(zpr{j}))); shading flat
end

%inner sphere and shafts
xi = 0.35*OUTINR*xs;
yi = 0.35*OUTINR*ys;
zi = 0.35*OUTINR*zs;

[xsh ysh zsh] = cylinder(RSHAFT,50);

surf(xi,yi,zi,ISCOL*ones(size(zi))); shading flat
surf(xsh,ysh,zsh-1,ISCOL*ones(size(zsh))); shading flat
surf(xsh,ysh,1.4*zsh+0.5,ISCOL*ones(size(zsh))); shading flat

%make a rectangular tube for the magnet
rmag = [1.8*ones(1,10) linspace(MAGINR,MAGOUTR,10) MAGOUTR*ones(1,10) linspace(MAGOUTR,MAGINR,10)]; %set of points on rectangle ... magnet bdry xsec
zmag = [linspace(ZMAGB,ZMAGT,10) ZMAGT*ones(1,10) linspace(ZMAGT,ZMAGB,10) ZMAGB*ones(1,10)];
theta = [0:pi/100:2*pi];
Mxm = zeros(length(theta),length(rmag));
Mym = zeros(length(theta),length(rmag));
Mzm = zeros(length(theta),length(rmag));
for j = 1:length(theta)
	Mxm(j,:) = rmag*cos(theta(j));
	Mym(j,:) = rmag*sin(theta(j));
	Mzm(j,:) = zmag;
end
%make a rectangular tube for the flange
rf = [FLANGEINR*ones(1,10) linspace(FLANGEINR,FLANGEOUTR,10) FLANGEOUTR*ones(1,10) linspace(FLANGEOUTR,FLANGEINR,10)]; %set of points on rectangle ... magnet bdry xsec
zf = [linspace(FLANGEBOTTOM,FLANGETOP,10) FLANGETOP*ones(1,10) linspace(FLANGETOP,FLANGEBOTTOM,10) FLANGEBOTTOM*ones(1,10)];
theta = [0:pi/100:2*pi];
Mxf = zeros(length(theta),length(rf));
Myf = zeros(length(theta),length(rf));
Mzf = zeros(length(theta),length(rf));
for j = 1:length(theta)
	Mxf(j,:) = rf*cos(theta(j));
	Myf(j,:) = rf*sin(theta(j));
	Mzf(j,:) = zf;
end

%make a lid
rlid = [RSHAFT*ones(1,10) linspace(RSHAFT,FLANGEOUTR,30) FLANGEOUTR*ones(1,10) linspace(FLANGEOUTR,RSHAFT,30)]; %set of points on rectangle ... magnet bdry xsec
zlid = [linspace(FLANGETOP,LIDTOP,10) LIDTOP*ones(1,30) linspace(FLANGETOP,FLANGEBOTTOM,10) FLANGEBOTTOM*ones(1,30)];
theta = [0:pi/100:2*pi];
Mxl = zeros(length(theta),length(rlid));
Myl = zeros(length(theta),length(rlid));
Mzl = zeros(length(theta),length(rlid));
for j = 1:length(theta)
	Mxl(j,:) = rlid*cos(theta(j));
	Myl(j,:) = rlid*sin(theta(j));
	Mzl(j,:) = zlid;
end

%make the bell housing
BHBOT=1.8;
BHTOP = 2.1;
BHR = 0.16;
BHHOLE = 0.01;
BHCOL=-0.9
rbh = [BHHOLE*ones(1,10) linspace(BHHOLE,BHR,10) BHR*ones(1,10) linspace(BHR,BHHOLE,10)]; 
zbh = [linspace(BHBOT,BHTOP,10) BHTOP*ones(1,10) linspace(BHBOT,BHTOP,10) BHBOT*ones(1,10)];
theta = [0:pi/100:2*pi];
Mxb = zeros(length(theta),length(rbh));
Myb = zeros(length(theta),length(rbh));
Mzb = zeros(length(theta),length(rbh));

for j = 1:length(theta)
	Mxb(j,:) = rbh*cos(theta(j));
	Myb(j,:) = rbh*sin(theta(j));
	Mzb(j,:) = zbh;
end


%make the torque sensor
TSBOT=BHTOP;
TSTOP = BHTOP+0.15;
TSR = 0.075;
TSHOLE = 0.03;
TSCOL=PROBECOL+0.3;
rts = [TSHOLE*ones(1,10) linspace(TSHOLE,TSR,10) TSR*ones(1,10) linspace(TSR,TSHOLE,10)]; 
zts = [linspace(TSBOT,TSTOP,10) TSTOP*ones(1,10) linspace(TSBOT,TSTOP,10) TSBOT*ones(1,10)];
theta = [0:pi/100:2*pi];
Mxt = zeros(length(theta),length(rts));
Myt = zeros(length(theta),length(rts));
Mzt = zeros(length(theta),length(rts));

for j = 1:length(theta)
	Mxt(j,:) = rts*cos(theta(j));
	Myt(j,:) = rts*sin(theta(j));
	Mzt(j,:) = zts;
end

discr = linspace(0.001,1,10);
diskx = zeros(length(theta),length(discr));
disky = zeros(length(theta),length(discr));
for j = 1:length(theta);
	diskx(j,:) = discr*cos(theta(j));
	disky(j,:) = discr*sin(theta(j));
end

PORTR = 0.16/2;
portx = PORTR*diskx;
porty = PORTR*disky;
portposx = [0 0.6 0 -0.6];
portposy = [0.6 0 -0.6 0];
portposz = [LIDTOP LIDTOP LIDTOP LIDTOP]+0.01;
thetabc = [0:pi/24:2*pi-pi/24];
NBOLTS = length(thetabc)
bx = 0.0254*diskx/2;
by = 0.0254*disky/2;
BOLTR = 0.8255;
bpx = BOLTR*cos(thetabc);
bpy = BOLTR*sin(thetabc);
bpz = ones(size(thetabc))*(LIDTOP+0.01);




%plot magnet, flange, lid, torque coupler, ports, bolt holes
surf(Mxm,Mym,Mzm,WIRECOL*ones(size(Mzm))); %magnet
surf(Mxf,Myf,Mzf,SPHERECOL*ones(size(Mzf))); %flange
surf(Mxl,Myl,Mzl,SPHERECOL*ones(size(Mzl))); %lid
surf(Mxb,Myb,Mzb,BHCOL*ones(size(Mzb))); %bell housing
surf(Mxt,Myt,Mzt,TSCOL*ones(size(Mzt))); %torque sensor

for j = 1:length(portposx);
	surf(portx+portposx(j),porty+portposy(j),portposz(j)*ones(size(portx)),PROBECOL*portposz(j)*ones(size(portx)));
end
shading flat

for j = 1:length(bpx);
	surf(bx+bpx(j),by+bpy(j),bpz(j)*ones(size(bx)),zeros(size(bx)));
end
%set angle, viewport, and lighting
set(gca,'view',[-20 30]);
camlight
xlim([-4 4]);
zlim([-4 4]);
ylim([-4 4]);
zoom(3);

