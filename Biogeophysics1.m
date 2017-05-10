function [qg,soilbeta,btran,btran2] = ...
    Biogeophysics1(ProxySoilWater,om_frac,sand,clay,frac_sno,t_grnd,forc_pbot,forc_q)
%{
Variable wx (partial volume of ice and water of surface layer) is unitless
[m/m] and could be approximated by the height of ground water level, which
is available for Alptal. However, it still will probably not be a sufficient
solution since other observational sites might not have measurements of soil
water, etc.
Variables
- watsat (volumetric soil water at saturation),
- sucsat (minimum soil suction [mm]),
- bsw (Clapp and Hornberger "b"),
- smpmin (restriction for min of soil potential [mm]) and
- hksat (hydraulic conductivity at saturation [mm H2O s^{-1}])
are calculated from soil composition data.
Calculations partially taken from iniTimeConst.f90.
%}
grav = 9.80616; % acceleration of gravity [m s^{-2}]

watsat = (1-om_frac)*(0.489-0.00126*sand*100) + 0.88*om_frac;    % 0.88 - om_watsat -> [0.83,0.93]
sucsat = (1-om_frac)*(10*10^(1.88-0.0131*sand*100)) + 10.2*om_frac;    % 10.2 - om_sucsat -> [10.1,10.3]
bsw = (1-om_frac)*(2.91+0.159*clay*100) + om_frac*5.35;    % 5.35 - om_b -> [2.7,12]
xksat = 0.0070556*10^(-0.884+0.0153*sand*100);
if om_frac>0.5
    perc_frac = (1-0.5)^(-0.139) * (om_frac-0.5)^0.139;
else
    perc_frac = 0;
end
if om_frac<1
    uncon_hksat = ((1-om_frac)+(1-perc_frac)*om_frac)...
        /((1-om_frac)/xksat + ((1-perc_frac)*om_frac)/0.1);
else
    uncon_hksat = 0;
end
hksat = ((1-om_frac)+(1-perc_frac)*om_frac)*uncon_hksat + (perc_frac*om_frac)*0.14;    % 0.14 - om_hksat -> [0.0001,0.28]
watfc = watsat*(0.1/(hksat*86400))^(1/(2*bsw+3));
smpmin = -1*10^8;
wx = ProxySoilWater*watsat; % Proxy for Alptal just unitless fraction, not partial volume
fac = max(0.01,min(1,wx/watsat)); % wx relative to watsat <-'
psit = max(smpmin,-sucsat*fac^(-bsw));
roverg = (6.02214*10^26*1.38065*10^(-23)/18.016)/grav * 1000; % [mm K^{-1}]
hr = exp(psit/roverg/t_grnd);
qred = (1-frac_sno)*hr + frac_sno;  % frac_h2osfc disregarded
if wx < watfc
    fac_fc = max(0.01,min(1,wx/watfc));
    soilbeta = (1-frac_sno)*0.25*(1-cos(pi*fac_fc))^2 + frac_sno;  % frac_h2osfc disregarded
else
    soilbeta = 1;
end
soilalpha = qred;
[eg,degdT,qsatg,qsatgdT] = QSat(t_grnd,forc_pbot);
qg = qred*qsatg;
dqgdT = qred*qsatgdT;
if qsatg > forc_q && forc_q > qred*qsatg   % only works for frac_sno either 1 or 0...
    qg = forc_q;                           % ...otherwise split into snow and soil cases
    dqgdT = 0;
end


%%%%%%%%%%%%   moved from CanopyFluxes.m to Biogeophysics1.m   %%%%%%%%%%%%
% transpiration wetness factor (0 to 1)
btran  = 0;
btran2  = 0;
%{
"Effective porosity of soil, partial volume of ice and liquid (needed for
btran) and root resistance factors"
...are calculated here to get btran, which is necessary for the calculation
of transpiration, etc. later on. A value of btran=0 would limit the
calculations and thus be disadvantageous for the assessment of CLM4.5.
Therefore, btran is calculated out of proxy values describing the soil
water content.
%}
smpsc = -255000;    % soil water potential at full stomatal closure [mm], value for NBTs, -224000 for DBTs
smpso = -66000;     % soil water potential at full stomatal opening [mm], value for NBTs, -35000 for DBTs
tfrz = 273.15;
% conversion of water content proxy to fractions of solid and liquid water
h2osoi_vol = watsat*ProxySoilWater;
%{
if t_grnd <= tfrz   % temperatures and water content per layer in CLM -> how to include here?
    vol_ice_proxy = h2osoi_vol;
    vol_liq_proxy = 0;
else
    vol_ice_proxy = 0; 
    vol_liq_proxy = h2osoi_vol;
end
%}
% temperatures and water content per layer in CLM -> here assumed surface
% temperature has to be quite low for soil to freeze completely
if t_grnd <= tfrz-15
    vol_ice_proxy = h2osoi_vol;
    vol_liq_proxy = 0;
    htvp = 2.501*10^6 + 3.337*10^5; % arbitrarily assume that sublimation occurs only as h2osoi_liq = 0
elseif t_grnd > tfrz-15 && t_grnd <= tfrz
    vol_ice_proxy = h2osoi_vol*(1 - (t_grnd - (tfrz-15))/15);
    vol_liq_proxy = h2osoi_vol*(t_grnd - (tfrz-15))/15;
    htvp = 2.501*10^6;
else
    vol_ice_proxy = 0; 
    vol_liq_proxy = h2osoi_vol;
    htvp = 2.501*10^6;
end

% Root resistance factors
vol_ice = min(watsat,vol_ice_proxy);
eff_porosity = watsat-vol_ice;
vol_liq = min(eff_porosity,vol_liq_proxy);
rootfr = 0.9949; % since no differentiation of soil layers -> all of soil considered, but due to design total rootfraction not 1
%{
roota_par = 7;      % CLM rooting distribution parameter [1/m], value for NBTs, 6 for DBTs
rootb_par = 2;      % CLM rooting distribution parameter [1/m], value for all boreal trees
%}
if vol_liq <= 0 || t_grnd <= tfrz-15 % tfrz-15 -> as above, but in CLM4.5 tfrz-2 for every layer
    rootr = 0;
else
    s_node = max(vol_liq/eff_porosity,0.01);
    smp_node = max(smpsc,-sucsat*s_node^(-bsw));
    rresis = min((eff_porosity/watsat)*(smp_node - smpsc)/(smpso - smpsc),1);
    rootr = rootfr*rresis;
    btran = btran + rootr;
    smp_node_lf = max(smpsc,-sucsat*(h2osoi_vol/watsat)^(-bsw)) ;
    btran2 = btran2 + rootfr*min((smp_node_lf - smpsc)/(smpso - smpsc),1);
end

%{
% Normalize root resistances to get layer contribution to ET
if btran > 0
    rootr = rootr/btran;
else
    rootr = 0;
end
%}

end