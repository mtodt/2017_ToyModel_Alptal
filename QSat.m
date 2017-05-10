function [es,esdT,qs,qsdT] = QSat(T,p)
%{
DESCRIPTION:
Computes saturation mixing ratio and the change in saturation mixing ratio
with respect to temperature.
Reference:  Polynomial approximations from:
            Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
            vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
%}

% for water vapor (temperature range 0C-100C)
a0 =  6.11213476;
a1 =  0.444007856;
a2 =  0.143064234*10^(-1);
a3 =  0.264461437*10^(-3);
a4 =  0.305903558*10^(-5);
a5 =  0.196237241*10^(-7);
a6 =  0.892344772*10^(-10);
a7 = -0.373208410*10^(-12);
a8 =  0.209339997*10^(-15);

% for derivative:water vapor
b0 =  0.444017302;
b1 =  0.286064092*10^(-1);
b2 =  0.794683137*10^(-3);
b3 =  0.121211669*10^(-4);
b4 =  0.103354611*10^(-6);
b5 =  0.404125005*10^(-9);
b6 = -0.788037859*10^(-12);
b7 = -0.114596802*10^(-13);
b8 =  0.381294516*10^(-16);

% for ice (temperature range -75C-0C)
c0 =  6.11123516;
c1 =  0.503109514;
c2 =  0.188369801*10^(-1);
c3 =  0.420547422*10^(-3);
c4 =  0.614396778*10^(-5);
c5 =  0.602780717*10^(-7);
c6 =  0.387940929*10^(-9);
c7 =  0.149436277*10^(-11);
c8 =  0.262655803*10^(-14);

% for derivative:ice
d0 =  0.503277922;
d1 =  0.377289173*10^(-1);
d2 =  0.126801703*10^(-2);
d3 =  0.249468427*10^(-4);
d4 =  0.313703411*10^(-6);
d5 =  0.257180651*10^(-8);
d6 =  0.133268878*10^(-10);
d7 =  0.394116744*10^(-13);
d8 =  0.498070196*10^(-16);


T_limit = T - 273.15;
if T_limit > 100
    T_limit=100;
end
if T_limit < -75
    T_limit = -75;
end

td = T_limit;
if td >= 0
    es = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 + td*(a5 + td*(a6 + td*(a7 + td*a8)))))));
    esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 + td*(b5 + td*(b6 + td*(b7 + td*b8)))))));
else
    es = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 + td*(c5 + td*(c6 + td*(c7 + td*c8)))))));
    esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 + td*(d5 + td*(d6 + td*(d7 + td*d8)))))));
end

es = es*100;        % [Pa]
esdT = esdT*100;    % [Pa K^{-1}]

vp = 1/(p - 0.378*es);
vp1 = 0.622*vp;
vp2 = vp1*vp;

qs = es*vp1;        % [kg kg^{-1}]
qsdT = esdT*vp2*p;  % [k^{-1}]
end