%type can be :
% 1 == rule set one ==
%sSS, sTT, tST
% 2 == rule set two ==
%sTS, sST, tSS, tTT


function [selm llist]= selrules(lmax,type,PLOTFLAG,lvobs,mvobs,lBobs,mBobs);

l = zeros(lmax*(lmax+3)/2,1);
m = zeros(size(l));
for ll = 1:lmax
	ls = (ll-1)*(ll+3-1)/2+1;
	lt = (ll)*(ll+3)/2;
	l(ls:lt) = ll;
	m(ls:lt) = 0:ll;
end

%varargin should be four variables, 
if nargin==3
	lvobs=l;
	lBobs=l;
	mvobs=m;
	mBobs=m;
end

if ~((nargin==3) | (nargin==7))
	error('Usage: [selectionmatrix connectionlist] = selrules(lmax,typestring,PLOTFLAG,{lvobs,mvobs,lBobs,mBobs}');
end

if ( strcmp(type,'sSS') | strcmp(type,'sTT') | strcmp(type,'tST'))
	
	tn = type;
	type =1;
	['Bullard and Gellman Case 1: Vel ' tn(1) ' Coupling Mag ' tn(2) ' to Mag ' tn(3)]
elseif (strcmp(type,'sTS') | strcmp(type,'sST') | strcmp(type,'tSS') | strcmp(type,'tTT'))
	
	tn = type;
	type =2;
	['Bullard and Gellman Case 2: Vel ' tn(1) ' Coupling Mag ' tn(2) ' to Mag ' tn(3)]
elseif strcmp(type,'tTS')
	error('tTS is always zero!')
else
	error('Invalid Spherical Harmonic Type Specification')
end
	


for j = 1:length(lvobs)
	axlab1{j} = [tn(1) '_' num2str(lvobs(j)) '^' num2str(mvobs(j))];
end
for j = 1:length(lBobs)
	axlab2{j} = [tn(2) '_' num2str(lBobs(j)) '^' num2str(mBobs(j))];
end
for j = 1:length(l)
	axlab3{j} = [tn(3) '_' num2str(l(j)) '^' num2str(m(j))];
end

[L1 L2 L3] = ndgrid(lvobs,lBobs,l);
[M1 M2 M3] = ndgrid(mvobs,mBobs,m);

Lsize = size(L1)
Msize = size(M1)


switch type
	case 1 %Bullard and Gellman (1): sSS, sTT, tST

	'case 1 switchcase' 
	selm = (~mod(L1+L2+L3,2))  & ((L1+L2>=L3) & (L2+L3>=L1) & (L1+L3>=L2)) & (((M1+M2+M3)==0)|((M1-M2+M3)==0)|((M1+M2-M3)==0)|((M1-M2-M3==0)));
	Selsize= size(selm)
	llist=[];
	for j1 = 1:length(lvobs)
		for j2 = 1:length(lBobs)
			for j3 = 1:length(l)
				if selm(j1,j2,j3)
					llist = [llist; ['$$' tn(1) '_' num2str(lvobs(j1)) '^' num2str(mvobs(j1)) tn(2) '_' num2str(lBobs(j2)) '^' num2str(mBobs(j2)) tn(3) '_' num2str(l(j3)) '^' num2str(m(j3)) '$$']];
				end
			end
		end
	end

	case 2 %Bullard and Gellman (2): sTS, sST, tSS, tTT
	'case 2 switchcase'
	selm = (mod(L1+L2+L3,2))  & ((L1+L2>=L3) & (L2+L3>=L1) & (L1+L3>=L2)) & (((M1+M2+M3)==0)|((M1-M2+M3)==0)|((M1+M2-M3)==0)|((M1-M2-M3==0)));
	%in this section, identical harmonics are not allowed
	selm = selm & ~(((L1==L2) & (M1==M2) & (upper(tn(1))==upper(tn(2))) ) | ((L2==L3) & (M2==M3) & (upper(tn(2))==upper(tn(3)))) | ((L1==L3) & (M1==M3) & (upper(tn(1))==upper(tn(3)))));
	llist=[];
	for j1 = 1:length(lvobs)
		for j2 = 1:length(lBobs)
			for j3 = 1:length(l)
				if selm(j1,j2,j3)
					llist = [llist; ['$$' tn(1) '_' num2str(lvobs(j1)) '^' num2str(mvobs(j1)) tn(2) '_' num2str(lBobs(j2)) '^' num2str(mBobs(j2)) tn(3) '_' num2str(l(j3)) '^' num2str(m(j3)) '$$']];
					
				end
			end
		end
	end

end %of switch/case



if PLOTFLAG
[mm nn oo] = size(L1);
[x y z] = meshgrid(1:mm,1:nn,1:oo);
maxr = sqrt(3)*lmax; 
nanm = (selm+0);
nanm(~selm)=NaN;
figure; scatter3(nanm(:).*x(:),y(:),z(:),40*selm(:)+1,sqrt(L1(:).^2+L2(:).^2+L3(:).^2)./lmax-0.5,'filled')
set(gca,'xtick',1:length(l),'ytick',1:length(l),'ztick',1:length(l),'xticklabel',axlab1,'yticklabel',axlab2,'zticklabel',axlab3)
axis equal
xlabel(tn(1),'fontsize',18); ylabel(tn(2),'fontsize',18); zlabel(tn(3),'fontsize',18,'rotation',0);
end
