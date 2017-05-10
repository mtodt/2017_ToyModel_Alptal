function bla = StabFunc2(zeta)
%{
DESCRIPTION:
Stability function for rib < 0.
%}
chik2 = sqrt(1 - 16*zeta);

bla = 2*log((1+chik2)/2);

end