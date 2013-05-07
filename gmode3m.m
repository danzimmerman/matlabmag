function ff = gmode3m(n,r,theta,phi);
if ~(nargin==4)
    error('Usage: ff = gmode3m(n,r,theta,phi)')
end

switch n

case 1 
l=1; m=0; %g10
phifun = 1; %cos(m*phi) for m=0

case 2
l=1; m=1; %g11c
phifun = cos(m*phi);

case 3
l=1; m=1; %g11s
phifun = sin(m*phi);

case 4
l=2; m=0; %g20
phifun = 1;

case 5
l=2; m=1; %g21c
phifun = cos(m*phi);

case 6
l=2; m=1; %g21s
phifun = sin(m*phi);

case 7
l=2; m=2; %g22c
phifun = cos(m*phi);

case 8
l=2; m=2; %g22s
phifun = sin(m*phi);

case 9
l=3; m=0; %g30
phifun = 1;

case 10
l=3; m=1; %g31c
phifun = cos(m*phi);

case 11
l=3; m=1; %g31s
phifun = sin(m*phi);

case 12
l=3; m=2; %g32c
phifun = cos(m*phi);

case 13
l=3; m=2; %g32s
phifun = sin(m*phi);

case 14
l=3; m=3; %g33c
phifun = cos(m*phi);

case 15
l=3; m=3; %g33s
phifun = sin(m*phi);

case 16
l=4; m=0; %g40
phifun = 1;

case 17
l=4; m=1; %g41c
phifun = cos(m*phi);

case 18
l=4; m=1; %g41s
phifun = sin(m*phi);

case 19
l=4; m=2; %g42c
phifun = cos(m*phi);

case 20
l=4; m=2; %g42s
phifun = sin(m*phi);

case 21
l=4; m=3; %g43c
phifun = cos(m*phi);

case 22
l=4; m=3; %g43s
phifun = sin(m*phi);

case 23
l=4; m=4; %g44c
phifun = cos(m*phi);

case 24
l=4; m=4; %g44s
phifun = sin(m*phi);

case 25
l=5; m=0; %g50s
phifun = 1;

case 26
l=5; m=1; %g51c
phifun = cos(m*phi);

case 27
l=5; m=1; %g51s
phifun = sin(m*phi);

case 28
l=5; m=2; %g52c
phifun = cos(m*phi);

case 29
l=5; m=2; %g52s
phifun = sin(m*phi);

case 30
l=5; m=3; %g53c
phifun = cos(m*phi);

case 31
l=5; m=3; %g53s
phifun = sin(m*phi);

case 32
l=5; m=4; %g54c
phifun = cos(m*phi);

case 33
l=5; m=4; %g54s
phifun = sin(m*phi);

case 34
l=5; m=5; %g55c
phifun = cos(m*phi);

case 35
l=5; m=5; %g55s
phifun = sin(m*phi)

case 36
l=6; m=0; %g60s
phifun = 1;

case 37
l=6; m=1; %g61c
phifun = cos(m*phi);

case 38
l=6; m=1; %g61s
phifun = sin(m*phi);

case 39
l=6; m=2; %g62c
phifun = cos(m*phi);

case 40
l=6; m=2; %g62s
phifun = sin(m*phi);

case 41
l=6; m=3; %g63c
phifun = cos(m*phi);

case 42
l=6; m=3; %g63s
phifun = sin(m*phi);

case 43
l=6; m=4; %g64c
phifun = cos(m*phi);

case 44
l=6; m=4; %g64s
phifun = sin(m*phi);

case 45
l=6; m=5; %g65c
phifun = cos(m*phi);

case 46
l=6; m=5; %g65s
phifun = sin(m*phi);

case 47
l=6; m=6; %g66c
phifun = cos(m*phi);

case 48
l=6; m=6; %g66s
phifun = sin(m*phi);




end %of switch/case

legmatrix = legendre(l,cos(theta),'sch');
Plm = squeeze(legmatrix(m+1,:,:)); %m+1th row is the desired one

ff = l*(l+1)*(r.^-(l+2)).*Plm.*phifun*(-1)^m;

end %of function

