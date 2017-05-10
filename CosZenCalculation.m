function [coszen,Declin] = CosZenCalculation(JulDay,Lat,Lon)
%{
Return the cosine of the solar zenith angle. Assumes 365.0 days/year and 
JulDay equals fractional days of the year from 1.xx to 365.xx.

Equation for declination taken from technical notes for LSM1 by Bonan (1996).
%}

Lat = pi*Lat/180;
Lon = pi*Lon/180;

Rad = 2*pi*JulDay/365;

Declin = 0.006918 - 0.399912*cos(Rad) + 0.070257*sin(Rad) - 0.006758*cos(2*Rad)...
    + 0.000907*sin(2*Rad) - 0.002697*cos(3*Rad) + 0.00148*sin(3*Rad);

SolarHour = JulDay*2*pi + Lon;    % if JulDay with respect to Greenwich Time
SolarHour = JulDay*2*pi;          % if JulDay already local time

coszen = sin(Lat)*sin(Declin) - cos(Lat)*cos(Declin)*cos(SolarHour);

end