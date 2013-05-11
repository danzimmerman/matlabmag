function [] = unwind(xt,yt,zt, fname)

%if fname specified ends in '.txt' write as ASCII with dlmwrite
%	otherwise write as 8byte floats in binary 
%read binary file in python with x = np.fromfile(fname)
%by David Meichle

if( fname( (end-3):end)  == '.txt') 
ASCII = true;
else
ASCII = false; 
end

%make single column
xt = xt(:); 
yt = yt(:); 
zt = zt(:); 

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
		xvals = row(inds); 
		count = length(xvals); 
		x(s,:) = repmat(count, [1,3]); 
		rng = (s+1):(s+count); 
		x(rng,c) = xvals; 
		s = rng(end) + 1; 
	end
end

%crop to actual data size. size(x,1) will always be <= MAXNT*NLINES
x = x(1:(s-3), :);

%write to file
if(ASCII)
	dlmwrite(fname, x);
	disp(['Done writing to ASCII file: ', fname]);
else
	fid = fopen(fname,'w');
	fwrite(fid, x, 'double');
	fclose(fid);
	disp(['Done writing to binary double file: ', fname]);
end

return 

