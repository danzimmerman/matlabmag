function [r theta phi Br Btheta Bphi] = bcart2bsph(x,y,z,Bx,By,Bz);

r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi = atan2(y,x);

%Br = B dot rhat = Bx*rhat_x + By*rhat_y + Bz*rhat_z etc

Br = Bx.*sin(theta).*cos(phi)+By.*sin(theta).*sin(phi)+Bz.*cos(theta);
Btheta = Bx.*cos(theta).*cos(phi)+By.*cos(theta).*sin(phi)-Bz.*sin(theta);
Bphi = -Bx.*sin(phi)+By.*cos(phi); 

end



