%this version integrates B if the initial point has B_r > 0 and -B if the initial point has B_r<0
%this sums a g vector
%change history 050313 started history written by Axl
% this version is vectorized over field lines for speed

%x0,y0,z0 should be all one size. for now they can be 1D vector lists or 2D like meshgrid
%g is a vector of 24 gauss coefficients

function [xt yt zt Bxt Byt Bzt xlast ylast zlast] = bline_spag_reverse_gauss(g,x0,y0,z0,stepsize,maxnumberofsteps,fieldintdir,MINRADIUS,MAXRADIUS);
if ~(nargin==9)
	error('Usage [xt yt zt Bxt Byt Bzt xlast ylast zlast]=bline_spag_reverse_gauss(g,x0,y0,z0,stepsize,maxnumberofsteps,fieldintdir,MINRADIUS,MAXRADIUS)');
end

DELT = stepsize;
SIZELINEBLOCK = size(x0);
FIELDINTDIR = fieldintdir; %1 to reverse Br<0 -1 to reverse Br>0)
NFIELDLINES = length(x0(:));
	%if (isrow(x0) | isrow(y0) | isrow(z0)) %this also tests for "is a vector" 
		%x0 = x0';
		%y0 = y0';
		%z0 = z0';
	%end
MAXRADVEC = MAXRADIUS*ones([1 SIZELINEBLOCK]);
MINRADVEC = MINRADIUS*ones([1 SIZELINEBLOCK]);
lastradius = zeros([1 SIZELINEBLOCK]);
NT = maxnumberofsteps;
xt = zeros([NT SIZELINEBLOCK]); %preallocate  
yt = zeros([NT SIZELINEBLOCK]);
zt = zeros([NT SIZELINEBLOCK]);

Bxt = zeros([NT SIZELINEBLOCK]);
Byt = zeros([NT SIZELINEBLOCK]);
Bzt = zeros([NT SIZELINEBLOCK]);
% variables used at each timestep

Bmag = zeros([1 SIZELINEBLOCK]);
vx = zeros([1 SIZELINEBLOCK]);
vy = zeros([1 SIZELINEBLOCK]);
vz = zeros([1 SIZELINEBLOCK]);
xtemp = zeros([1 SIZELINEBLOCK]);
ytemp = zeros([1 SIZELINEBLOCK]);
ztemp = zeros([1 SIZELINEBLOCK]);
active_lines = logical(ones(size(MAXRADVEC)));
mrvs0 = size(MAXRADVEC);
acls0 = size(active_lines);
%to hold the termination points of lines at the point they become inactive
xlast = zeros([1 SIZELINEBLOCK]);
ylast = zeros([1 SIZELINEBLOCK]);
zlast = zeros([1 SIZELINEBLOCK]);
%initialize first timestep
xt(1,:,:,:) = x0; %if x0 is a vector or 2D, the extra ":" in (1,:,:,:) don't matter but will allow 1,2 or 3D arrays to be used
yt(1,:,:,:) = y0;
zt(1,:,:,:) = z0;

	[Bxt(1,active_lines) Byt(1,active_lines) Bzt(1,active_lines)] = gm3cart_gauss(g,xt(1,active_lines),yt(1,active_lines),zt(1,active_lines)); %initialize (Bx,By,Bz) at (x0,y0,z0); 


%code to reverse "velocity field" if point starts in a B_r<0 region
r0 = sqrt(xt(1,:,:,:).^2+yt(1,:,:,:).^2+zt(1,:,:,:).^2); 
theta0 = acos(zt(1,:,:,:)./r0);
phi0 = atan2(yt(1,:,:,:),xt(1,:,:,:));
Br0 = Bxt(1,:,:,:).*sin(theta0).*cos(phi0)+Byt(1,:,:,:).*sin(theta0).*sin(phi0)+Bzt(1,:,:,:).*cos(theta0);
idr = zeros([1 SIZELINEBLOCK]);%make direction matrices the size of the position blocks
idr(Br0>=0) = 1*FIELDINTDIR; 
idr(Br0<0) = -1*FIELDINTDIR;
Bmag(1,active_lines) = sqrt(Bxt(1,active_lines).^2+Byt(1,active_lines).^2+Bzt(1, active_lines).^2);


for tk = 1:(NT-1)
	if (mod(tk,floor(NT/100)) == 0) %display percentage done
		disp(num2str(tk/NT*100));
	end

%calculate the field


lastradius(1,active_lines) = sqrt(xt(tk,active_lines).^2+yt(tk,active_lines).^2+zt(tk,active_lines).^2);
active_lines =  (~isnan(lastradius)) & ((lastradius>=MINRADVEC) & (lastradius<=MAXRADVEC)); %a logical array with 1 if it's inside the bounding spheres, 0 otherwise


%make three velocity components from current field value (and idr initial direction matrix)
Bmag(1,active_lines) = sqrt(Bxt(tk,active_lines).^2+Byt(tk,active_lines).^2+Bzt(tk, active_lines).^2);
vx(1,active_lines) = idr(1,active_lines).*Bxt(tk,active_lines)./Bmag(1,active_lines);
vx(1,~(active_lines)) = NaN;
vy(1,active_lines) = idr(1,active_lines).*Byt(tk,active_lines)./Bmag(1,active_lines);
vy(1,~(active_lines)) = NaN;
vz(1,active_lines) = idr(1,active_lines).*Bzt(tk,active_lines)./Bmag(1,active_lines);
vz(1,~(active_lines)) = NaN;

%calculate temporary new positions
xtemp(1,:) = xt(tk,:)+vx(1,:)*DELT;
ytemp(1,:) = yt(tk,:)+vy(1,:)*DELT;
ztemp(1,:) = zt(tk,:)+vz(1,:)*DELT;

%break before calculating the field
inactive_number = sum(sum(sum(~(active_lines)))); %the total number of not active lines
if (inactive_number==NFIELDLINES) %if NO lines are active then break loop
	['no active lines left, breaking loop at ' num2str(tk) ' timesteps']
	lasttimestep = tk;
	xt(tk,~(active_lines)) = NaN;
	yt(tk,~(active_lines)) = NaN;
	zt(tk,~(active_lines)) = NaN;
	Bxt(tk,~(active_lines)) = NaN;
	Byt(tk,~(active_lines)) = NaN;
	Bzt(tk,~(active_lines)) = NaN;
	break
end

%calc field at temporary new positions
[Bxt(tk+1,active_lines) Byt(tk+1,active_lines) Bzt(tk+1,active_lines)] = gm3cart_gauss(g,xtemp(1,active_lines),ytemp(1,active_lines),ztemp(1,active_lines));

%remake three velocity components from current field value (and idr initial direction matrix)
Bmag(1,active_lines) = sqrt(Bxt(tk+1,active_lines).^2 + Byt(tk+1,active_lines).^2+Bzt(tk+1, active_lines).^2);
vx(1,active_lines) = idr(1,active_lines).*Bxt(tk+1,active_lines)./Bmag(1,active_lines);
vx(1,~(active_lines)) = NaN;
vy(1,active_lines) = idr(1,active_lines).*Byt(tk+1,active_lines)./Bmag(1,active_lines);
vy(1,~(active_lines)) = NaN;
vz(1,active_lines) = idr(1,active_lines).*Bzt(tk+1,active_lines)./Bmag(1,active_lines);
vz(1,~(active_lines)) = NaN;

xt(tk+1,active_lines) = xt(tk,active_lines)+vx(1,active_lines)*DELT;
yt(tk+1,active_lines) = yt(tk,active_lines)+vy(1,active_lines)*DELT;
zt(tk+1,active_lines) = zt(tk,active_lines)+vz(1,active_lines)*DELT;
xt(tk,~(active_lines)) = NaN;
yt(tk,~(active_lines)) = NaN;
zt(tk,~(active_lines)) = NaN;

xlast(1,active_lines) = xt(tk+1,active_lines);
ylast(1,active_lines) = yt(tk+1,active_lines);
zlast(1,active_lines) = zt(tk+1,active_lines);


end

%finish with a whole line of NaNs but truncate matrices
xt = xt(1:lasttimestep,:,:,:);
yt = yt(1:lasttimestep,:,:,:);
zt = zt(1:lasttimestep,:,:,:);
Bxt = Bxt(1:lasttimestep,:,:,:);
Byt = Byt(1:lasttimestep,:,:,:);
Bzt = Bzt(1:lasttimestep,:,:,:);
xlast = squeeze(xlast);
ylast = squeeze(ylast);
zlast = squeeze(zlast);




end %of function
