function [] = unwind(xt,yt,zt, fname)
%
% David Meichle May 30 2013
%
% xt,yt,zt are input lines data
% fname specifies a file to store the lines data for python to read. 
% if fname ends in '.txt' it will write an ASCII file with dlmwrite, otherwise 
%	will print formatted binary 
%read binary file in python with x = np.fromfile(fname) 

if( fname( (end-3):end)  == '.txt') 
ASCII = true;
else
ASCII = false; 
end

%make single column
%sxt = xt(:); 
%yt = yt(:); 
%zt = zt(:); 

M = zeros([size(xt),3]); 
M(:,:,1) = xt; 
M(:,:,2) = yt;
M(:,:,3) = zt;

NLINES = size(xt,2);
MAXNT = size(xt,1);

x = zeros([MAXNT*NLINES,3]); %declare variable to hold 1D representation. oversize now to maximum possible size, will crop later. 


s = 1; 
for c = 1:3 %choose between X,Y,Z
s = 1; %reassign s = 1 for loop on c over X,Y,Z coordinates

	for l = 1:NLINES
		row = M(:,l,c); 
		inds = find(~isnan(row));
		if length(inds) > 0
			xvals = row(inds); 
			count = length(xvals); 
			x(s,:) = repmat(count, [1,3]); 
			rng = (s+1):(s+count); 
			x(rng,c) = xvals; 
			s = rng(end) + 1; 
		end

	end
end
 
 x = x( 1:s, :); 
 
 %make strictly 1D for easy saving 
xx = [x(:,1)', x(:,2)', x(:,3)']';


%%%% read in python with:
	%x = np.fromfile(fname)
	%x = np.transpose(np.reshape(x, [3, np.shape(x)[0]/3]) )

%write to file
if(ASCII)
	dlmwrite(fname, x);
	disp(['Done writing to ASCII file: ', fname]);
else
	fid = fopen(fname,'w');
	fwrite(fid, xx, 'double');
	fclose(fid);
	disp(['Done writing to binary double file: ', fname]);
end

return 

