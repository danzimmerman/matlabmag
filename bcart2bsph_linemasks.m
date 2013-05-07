%this converts x,y,z field line trajectories and their field components into spherical
%also outputs masks based on radial field for drawing colored lines. 1 if Br is "p"ositive or "n"egative, NaN otherwise

function [r theta phi Br Btheta Bphi Bmag maskp maskn] = bcart2bsph(x,y,z,Bx,By,Bz);

r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi = atan2(y,x);

%Br = B dot rhat = Bx*rhat_x + By*rhat_y + Bz*rhat_z etc

Br = Bx.*sin(theta).*cos(phi)+By.*sin(theta).*sin(phi)+Bz.*cos(theta);
Btheta = Bx.*cos(theta).*cos(phi)+By.*cos(theta).*sin(phi)-Bz.*sin(theta);
Bphi = -Bx.*sin(phi)+By.*cos(phi);
Bmag = sqrt(Bx.^2+By.^2+Bz.^2);

maskp = zeros(size(x));
maskn = zeros(size(x));
maskp(Br>=0) = 1;
maskn(Br<=0) = 1;
maskp(Br<0) = NaN;
maskn(Br>0) = NaN;

end



