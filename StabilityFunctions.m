function [psi_m,psi_h] = StabilityFunctions(xi)

% STABILITY FUNCTIONS FOR MOMENTUM
if xi <= 0
    % unstable case from Paulsen et al 1970
    x = (1 - 19*xi)^0.25; % 19 from H�gstr�m, 1996
    psi_m = log((1+x^2)*(1+x)^2/8) - 2*atan(x) + pi/2;
else
    % stable case from Holstlag and Bruin, following Beljaars & Holtslag 1991
    a = 1;
	b = 2/3;
    c = 5;
    d = 0.35;
	psi_m = -(a*xi + b*(xi - c/d)*exp(-d*xi) + b*c/d);
end

% STABILITY FUNCTIONS FOR HEAT
if xi <= 0
    % unstable case. Integrated by Paulsen, 1970 from phi-functions found by Businger et al, 1971
    x = (1 - 11.6*xi)^0.25; % 11.6 from H�gstr�m, 1996
    psi_h = 2*log((1+x^2)/2);
else
    % stable case, func=1, equation from Holtslag and De Bruin following Beljaars & Holstlag, 1991
    a = 1;
	b = 2/3;
    c = 5;
    d = 0.35;
	psi_h = -((1 + a*xi*2/3)^1.5) + b*(xi - c/d)*exp(-d*xi) + b*c/d - 1;
end





end