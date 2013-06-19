%naprop outputs density, kinematic viscosity, and electrical conductivity as a
%function of temperature in degrees celcius
%expressions taken from Fink and Leibowitz 1995
%"Thermodynamic and transport properties of sodium liquid and vapor"
%Argonne Report ANL/RE--95/2
%dynamic viscosity Section 2.2 Eq. 1
%density Eq. 18 page 20
%written by Axl 6/19/13

function [rho nu] = naprop(TNa);
if ~(nargin==1)
	error('Usage: [rho nu] = naprop(TNa) ([density viscosity] = naprop(sodium_temperature_in_Celcius))');
end

TK = TNa+273.16;
Tc = 2503.7;
rhoc = 219;
densf = 275.32;
densg = 511.58;
densh = 0.5;

%calculate density
rho = rhoc + densf*(1-TK/Tc)+densg*(1-TK/Tc).^densh; %in kg/m^3

%calculate dynamic viscosity (remember log is natural log)

mu = exp(-6.4406-0.3958*log(TK)+556.835./TK); %in Pa*s

%calculate kinematic viscosity = dynamic viscosity/density

nu = mu./rho; %in m^2/s

end %of function


