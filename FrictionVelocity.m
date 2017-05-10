function [ustar,temp1,temp2,temp12m,temp22m] = FrictionVelocity(displa,z0m,z0h,obu,um,forc_hgt_u_pft,forc_hgt_t_pft)
%{
DESCRIPTION:
Calculation of the friction velocity, relation for potential temperature
and humidity profiles of surface boundary layer.
%}

vkc = 0.4;      % von Karman constant
zetam = 1.574;  % transition point of flux-gradient relation (wind profile)
zetat = 0.465;  % transition point of flux-gradient relation (temperature profile)

% wind profile
zldis = forc_hgt_u_pft - displa;
zeta = zldis/obu;
if zeta < -zetam
    ustar = vkc*um/(log(-zetam*obu/z0m) - StabFunc1(-zetam) + StabFunc1(z0m/obu)...
        + 1.14*((-zeta)^(1/3) - zetam^(1/3)));
elseif zeta < 0
    ustar = vkc*um/(log(zldis/z0m) - StabFunc1(zeta) + StabFunc1(z0m/obu));
elseif zeta <= 1
    ustar = vkc*um/(log(zldis/z0m) + 5*zeta - 5*z0m/obu);
else
    ustar = vkc*um/(log(obu/z0m) + 5 - 5*z0m/obu + (5*log(zeta) + zeta - 1));
end
if zeta < 0
    vds_temp = 2*10^(-3) * ustar * (1 + (300/(-obu))^(2/3));
else
    vds_temp = 2*10^(-3) * ustar;
end
vds = vds_temp;

% calculation of 10m wind added for CLM4.5 -> not necessary here

% temperature profile
zldis = forc_hgt_t_pft - displa;
zeta = zldis/obu;
if zeta < -zetat
    temp1 = vkc/(log(-zetat*obu/z0h) - StabFunc2(-zetat) + StabFunc2(z0h/obu)...
        + 0.8*(zetat^(-1/3) - (-zeta)^(-1/3)));
elseif zeta < 0
    temp1 = vkc/(log(zldis/z0h) - StabFunc2(zeta) + StabFunc2(z0h/obu));
elseif zeta <= 1
    temp1 = vkc/(log(zldis/z0h) + 5*zeta - 5*z0h/obu);
else
    temp1 = vkc/(log(obu/z0h) + 5 - 5*z0h/obu + (5*log(zeta) + zeta - 1));
end


% humidity profile !!! assuming equal measurement heights and roughness lengths for temperature and humidity !!!
temp2 = temp1;


% temperature profile applied at 2m
zldis = 2 + z0h;
zeta = zldis/obu;
if zeta < -zetat
    temp12m = vkc/(log(-zetat*obu/z0h) - StabFunc2(-zetat) + StabFunc2(z0h/obu)...
        + 0.8*(zetat^(-1/3) - (-zeta)^(-1/3)));
elseif zeta < 0
    temp12m = vkc/(log(zldis/z0h) - StabFunc2(zeta) + StabFunc2(z0h/obu));
elseif zeta <= 1
    temp12m = vkc/(log(zldis/z0h) + 5*zeta - 5*z0h/obu);
else
    temp12m = vkc/(log(obu/z0h) + 5 - 5*z0h/obu + (5*log(zeta) + zeta - 1));
end


% humidity profile applied at 2m, s.a.
temp22m = temp12m;

end