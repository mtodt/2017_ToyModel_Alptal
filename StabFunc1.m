function bla = StabFunc1(zeta)
%{
DESCRIPTION:
Stability function for rib < 0.
%}
chik2 = sqrt(1 - 16*zeta);
chik = sqrt(chik2);

bla = 2*log((1+chik)/2) + log((1+chik2)/2) - 2*atan(chik) + pi/2;

end