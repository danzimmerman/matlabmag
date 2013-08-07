%this function takes a time series of hall sensors (which should ALREADY have the zero-field bias removed with getbias, etc!)
%and removes the mean DC field from the magnet measured with the 090612 0.1Hz rotating magnetic ramp
%this could also be modified to use the PREDICTED DC offsets in /data/tech_info/pbias.mat but the measured
%DC levels from applied field should be more accurate, especially near the equator
%ddb is the hall time series with zero-field biases subtracted, ivec should contain:

%ivec.mI : mean current at each step
%ivec.tstarts : ramp step start times
%ivec.tends : ramp step end times


function dout = debias_magramp(ddb,time,ivec);
if ~(nargin==3)
	error('Usage: dout = debias_magramp(ddb,time,ivec)')
end

temp = load('/data/3m/090612/MAT/ramp_index.mat');
mpatt = temp.v.BVpatt; %mean voltage pattern, volts of probe response per amp of magnet current for 31 probes
mps = size(mpatt)
dout = zeros(size(ddb));
for j = 1:length(ivec.mI);
k = find(time>ivec.tstarts(j) & time<ivec.tends(j));
dout(k,:) = ddb(k,:)-ivec.mI(j)*ones(length(k),1)*mpatt;
end

end %of function
