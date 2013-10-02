%here is an example analysis script 
RELOAD = 0;
DELAYFROMSTART = 60; %delay beyond tstarts for settling at beginning of ramp step
ro = 1.46; %outer sphere inner radius 1.46m
ri = 0.51; %inner sphere outer radius 0.51m 
LGAP = ro-ri; %gap width in m
ETA = 0.079; %magnetic diffusivity in m^2/s
MU0 =  1.25663706e-6; %permeability of free space in Tesla meters per amp
GPERA = 0.529; %B0 at reference location in Gauss per amp (0.529 for experiment center, single eq. coil)
GPERAS = 0.068; %Bs0 at stick location in Gauss per amp (0.068 for single eq. coil);
GPERAZ = 0.287;	%Bz0 at stick in Gauss per amp (0.287 for single eq. coil) 
HALLV = 0.0315; %nominal 31.5mV/Gauss for SS94A1F probes with 10V supply, could be a bit different?
DIRNAME = '/data/3m/092713/'; %full path to working directory WITH trailing slash
INDEXNAME = 'MAT/ramptimes_Ro6mr.mat'; %tstarts & tends index file with leading MAT/
DAQFILENAME = 'mr_i350_o050.daq'; %daq file name
BIASFILENAME = 'MAT/bias1.mat' %hall standing bias file with leading MAT/
TORQUECUTOFF = 0.5; %frequency cutoff in Hz for torque filter
SEALTORQUE = 16.5; %seal torque in N*m
PARAMCALC = 1; %Calculate all dim'less parameters.
PLOTS = 1; %plot some stuff automatically

if (~exist('vmr') | RELOAD)
	load([DIRNAME INDEXNAME])
end

%=== Calculate mean and std. sodium temperature and density and viscosity ===
if (~(isfield(vmr,'mtemp') & isfield(vmr,'stemp') & isfield(vmr,'dens') & isfield(vmr,'visc')) | RELOAD)
	'calculating mean temp, density, viscosity'
	wt = load([DIRNAME 'wtemp.log']);
	[vmr.mtemp vmr.stemp] = parmstd(wt(:,2),wt(:,1),vmr,0);
	[vmr.dens vmr.visc] = naprop(vmr.mtemp);
end

%=== Calculate mean and std. of magnet current ===
if (~(isfield(vmr,'mI') & isfield(vmr,'sI')) | RELOAD)
	'calculating mean and standard deviation of magnet current'
	ml = load([DIRNAME 'magnet.log']);
	[vmr.mI vmr.sI] = parmstd(ml(:,3),ml(:,1),vmr,DELAYFROMSTART); 
end

%=== Calculate mean and std. of torque (fluctuations below TORQUECUTOFF frequency!) ===
if (~(isfield(vmr,'mT') & isfield(vmr,'sT')) | RELOAD)
	['calculating mean and standard deviation of torque with ' num2str(TORQUECUTOFF) ' Hz cutoff']
	[bt at] = butter(2,TORQUECUTOFF/16,'low');
	load([DIRNAME 'torque.dat']);
	tqf = filtfilt(bt,at,torque(:,2));
	tqt = torque(:,1);
	tqfac = -1130/1801990; %1130 Newton Meters (maximum strain gauge load) at 1801190 bits (max ADC value is 2^21) 
	[vmr.mT vmr.sT] = parmstd(tqfac*tqf,tqt,vmr,DELAYFROMSTART);
end

%=== Calculate mean and std. of inner and outer frequencies as reported by the drives ===
if (~(isfield(vmr,'fo') & isfield(vmr,'fi')) | RELOAD)
	'calculating mean and standard deviation of drive-reported sphere frequencies'
	dl = load('/data/3m/092713/control.log');
	[vmr.fo vmr.sfo] = parmstd(48/400*dl(:,19),dl(:,1),vmr,DELAYFROMSTART); %column 19 in drivelog is MOTOR freq, mult by 48/400 for os freq
	[vmr.fi vmr.sfi] = parmstd(dl(:,13),dl(:,1),vmr,DELAYFROMSTART); %inner sphere direct drive
end

%=== Calculate magnetic quantities and dimensionless parameters and quantities for all ramp steps ===
if PARAMCALC
	'calculating various quantities'
	%Rossby number
	vmr.Ro = (vmr.fi-vmr.fo)./vmr.fo;
	%applied S-component at stick in Gauss calculated from Biot-Savart
	vmr.Bs0 = GPERAS*vmr.mI;
	%applied reference field at experiment center
	vmr.B0 = GPERA*vmr.mI;
	%applied reference field in Tesla
	vmr.B0T = GPERA*vmr.mI/10000;
	%Lundquist number $S = \frac{B_0 (r_o-r_i)}{\sqrt{\rho \mu_0}\eta}$
	vmr.S = LGAP*vmr.B0T./(sqrt(vmr.dens*MU0)*ETA); 
	%Alfven speed 
	vmr.vA = vmr.B0T./sqrt(vmr.dens*MU0);   
	%Lehnert number $\lambda = \frac{V_A}{\Omega (r_o-r_i)}$
	vmr.lambda = vmr.vA./(2*pi*abs(vmr.fo)*LGAP); 
	%Elsasser number $\Lambda = \frac{B_0^2}{\rho\mu_0\eta\Omega}$
	vmr.Lambda = (vmr.B0T.^2)./(vmr.dens*MU0*ETA*2*pi.*abs(vmr.fo));
	%Hartmann number $Ha = \frac{B_0 (r_o-r_i)}{\sqrt{\rho \mu_0 \eta \nu}}$
	vmr.Ha = (vmr.B0T*LGAP)./sqrt(vmr.dens*MU0*ETA.*vmr.visc);
	%Magnetic Reynolds number $Rm = \frac{\Delta\Omega (r_o-r_i)^2}{\eta}$
	vmr.Rm = (2*pi*abs(vmr.fi-vmr.fo)*LGAP^2)/ETA;
	%Magnetic Reynolds number $Re = \frac{\Delta\Omega (r_o-r_i)^2}{\nu}$
	vmr.Re = (2*pi*abs(vmr.fi-vmr.fo)*LGAP^2)./vmr.visc;
	%Ekman number $E = \frac{\nu}{\Omega (r_o-r_i)^2}$
	vmr.E = vmr.visc./(2*pi*vmr.fo*LGAP^2);
	%Dimensionless torque normalization factor:  G = Torque/gfac where gfac = \rho\nu^2 r_i
	vmr.gfac = (ri*vmr.dens.*vmr.visc.^2);
end

%=== Calculate mean and std. of Gauss coeffs, internal fields, and individual array probes. Gauss coeff & Hall array are debiased ===
if (~(isfield(vmr,'mGV') & isfield(vmr,'sGV') & isfield(vmr,'mstk') & isfield(vmr,'sstk') & isfield(vmr,'mhall') & isfield(vmr,'shall')) | RELOAD)

	if (~(exist('time') & exist('data') & exist('abstime')) | RELOAD)
		['loading DAQ file' DIRNAME DAQFILENAME]
		[data time abstime] = daqread([DIRNAME DAQFILENAME]);
	end
	
	if (~(exist('ddb') & exist('stickdb')) | RELOAD)
		'removing zero-field  array and internal stick biases'
		load([DIRNAME BIASFILENAME])
		preproc;
	end
	
	if (~(exist('dout')) | RELOAD)
		['removing applied contribution from magnet ramp using debias_magramp (090612 mean solid body ramp data)']
		dout = debias_magramp(ddb,time,vmr);
	end

	if (~(exist('g')) | RELOAD)
		'calculating Gauss coefficients from fully debiased array data'
		g = gcoeff3m(dout,probepos);
	end

	'calculating mean & std. of Gauss coeffs, internal fields, array probes'
	[vmr.mGV vmr.sGV] = parmstd(g,time,vmr,DELAYFROMSTART); %gauss coeff, APPLIED subtracted
	[vmr.mstk vmr.sstk] = parmstd(stickdb,time,vmr,DELAYFROMSTART); %bs and bphi APPLIED NOT SUBTRACTED
	[vmr.mhall vmr.shall] = parmstd(dout,time,vmr,DELAYFROMSTART); %mean and std hall probes with APPLIED subtracted


end



	
	%[sp f t p] = spectrogram(dout(:,16),100*256,80*256,100*256,256); %MLS1 spectrogram
	%figure; pcolor(t,f(1:1000)/0.5,log10(p(1:1000,:))); shading flat
	%[sp f t pg] = spectrogram(gn(:,12),100*256,80*256,100*256,256); %g32 spectrogram
	%figure; pcolor(t,f(1:1000)/0.5,log10(pg(1:1000,:))); shading flat
	%figure; plot(vmr.S,vmr.mstk(:,1)./(0.0315*vmr.Bs0),'v-'); %divide SS94A1F voltage by 0.0315 for B in gauss
	%hold on
	%plot(vmr.S,vmr.mstk(:,2)./(0.0315*vmr.Bs0),'o-');
	%legend('B_s','B_\phi');
	%bigaxrot('S','B_{int}/B_0');
	%gaussnormplot(vmr.S,vmr.sGV,[1 5 6 12 13 19 20 21 22],vmr.mI,1,1);
	%gn = oscomb(g,0.5,10,0.05,256);
	%[vmr.mGVn vmr.sGVn] = parmstd(gn,time,vmr,DELAYFROMSTART);
	
	
	

