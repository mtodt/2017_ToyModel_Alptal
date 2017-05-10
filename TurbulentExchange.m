function [ch_canopy,ce_canopy,ce_transpiration,ce_interception,ce_condensation,ra,rs]...
    = TurbulentExchange(SW_in_ac,T_air,wind,RelHum,T_can,LAI,z_can,refheight,z_0mg,frac_wet,SD,T_snowsurf)
%{
Turbulent exchange coefficients canopy model
The Shuttelworth and Wallace (1985) approach contains rather complicated
formulas compared to the general low level of understanding of within canopy
turbulent exchange (which was rightly commented by the referees to our paper
Towards an improved...). Also, the old stability correction very easily caused
oscillations in the canopy energy balance.
To meet the referees and to get rid of the oscillations, David introduced a
more simple method, based only on roughness lengths (from Blyth, 1999).
Roughness lengths are estimated from canopy height and leaf area index. A new
stability correction method based on Monin-Obukhov theory was also instead of
the Richardsson number formulation.
Stability correction now also below the canopy (using M-O).
%}

%-------------------------  Physical Parameters  --------------------------

karman = 0.4;       % von Kármán constant

gas_const_air = 287;% specific gas constant for (dry) air [J kg^{-1} K^{-1}]
rho_air = 1.1;      % density of air [kg m^{-3}] (Approximation: use ideal gas law)
C_air = 1004.67;    % specific heat capacity of air [J K^{-1}]

can_ch0 = 3;        % minimum heat exchange at zero wind [W m^{-2} K^{-1}]
rat_displ_can = 0.6667;     % ratio between canopy height and displacement height
rat_rough_can = 0.10;       % ratio between canopy height and roughness length
rat_heat_mom = 0.9999;      % ~=1, but not allowed to be exactly 1
can_rs_mult = 3.0;          % 1+can_rs_mult = max. factor to increase rs below canopy
rs_min = 375;       % minimum canopy surface resistance for transpiration
                    % 500 [s m^{-1}] is for needle leaf trees (van den Hurk et al., 2000),
                    % 75% (Gustafsson et al., 2003)
ra_incr_snow = 10;  % fractional increase of aerodynamic resistance for evaporation of intercepted snow
                    % 10 from Koivusalo and Kokkonen (2002), 8 from calibration with Alptal data
f3_gd = 0.0003;     % parameter for canopy surface resistance response to vapour pressure [Pa^{-1}]
                    % 0.0003 for trees (needle or broadleafs), 0 for crops, grass, tundra, etc.

% phase change constants
lh_subl = 2.838*10^6;        % latent heat of sublimation [J kg^{-1}]
lh_vap = 2.504*10^6;         % latent heat of vaporisation [J kg^{-1}]

% check wind speed to be at least 0.1 m/s
wind_local = max(wind,0.3);

% canopy height above snow surface
%cH = Xdata.cH - Xdata.Ground;       % calculated snow depth
z_can = z_can - SD;%cH;


%----------------------------  Calculations  ------------------------------

%{
1. displacement and roughness (mom) according to Shaw and Perreira (1981)
zdisplcan = 0.803 + 0.108 * CanDensMax - (0.462 - 0.086 * CanDensMax) *->
-> exp(-(0.163 + 0.283 * CanDensMax) * Cdata->lai);
zdisplcan = MAX (0., MIN (refheight - 0.5, zdisplcan * zcan));

1.3 roughness length
const double EQ1 = (0.175 - 0.098 * CanDensMax) + (-0.098 + 0.045 * CanDensMax) * log10(Cdata->lai);
const double EQ2 = (0.150 - 0.025 * CanDensMax) + (0.122 - 0.0135 * CanDensMax) * log10(Cdata->lai);
zomc = MIN(RoughLmax, MAX(zcan * MIN(EQ1, EQ2), RoughLmin)) * CAN_Z0M_COEF;

1. displacement and roughness as simple scaling of canopy height.
Please note:
1) Canopy roughness is not allowed to be smaller than roughness of
snow surface as defined by zomg, otherwise ustar may be negative!
This is now guaranteed by the computation of RoughLmin = MAX(0.01,zomg)
2) refheight is already given as height relative the snow surface
%}

% Shaw Perreira parameters
RoughLmin = 0.01;
RoughLmax = 100;
z_displ_can = max(0,min(refheight-0.5,rat_displ_can*z_can));
z_0mc = max(max(RoughLmin,z_0mg),min(RoughLmax,rat_rough_can*z_can));

% 2. aerodynamic resistances simple approach (Blyth, 1999)
    % 2.1 roughness length for scalars (heat and vapour)
z_0hc = rat_heat_mom*z_0mc;
z_0hg = rat_heat_mom*z_0mg;
	% update Cdata variables
z_0m = z_0mc;
z_0h = z_0hc;
z_displ = z_displ_can;

    % 2.2 Stability correction (adopted from Beljaars and Holtslag, 1991)
% psi_m = 0;
% psi_h = 0;
    % 2.2.1 get Aeta = Monin-Obukhov stability parameter from Richardson number
aeta = RichardsonToAeta(refheight-z_displ_can,T_air,T_air-T_can,wind_local,z_0mc,z_0hc);
[psi_m1,psi_h1] = StabilityFunctions(aeta);
[psi_m2,psi_h2] = StabilityFunctions(aeta*z_0hc/(refheight-z_displ_can));
[psi_m3,psi_h3] = StabilityFunctions(aeta*z_0mc/(refheight-z_displ_can));
psi_h = psi_h2 - psi_h1;
psi_m = psi_m3 - psi_m1;

    % 2.3 friction velocity above canopy
ustar = wind_local*karman/(log((refheight - z_displ_can)/z_0mc) + psi_m);

	% 2.4 transfer coefficient for scalars above canopy
ch_e = ustar*karman/(log((refheight - z_displ_can)/z_0hc) + psi_h);
ch = can_ch0/(rho_air*C_air) + ch_e;

	% 2.5 aerodynamic resistance above canopy
ra = 1/ch;
ra_e = 1/ch_e;

	% 2.6 canopy to canopy level resistance
if log(z_0mc/z_0hc) > 0
	rc = (log(z_0mc/z_0hc))/(karman*ustar);
else
    rc = 0;
end

	% 2.7 surface to canopy level resistance
if log(z_0mc/z_0hg) > 0
    rs = (log(z_0mc/z_0hg))/(karman*ustar) * (1 + can_rs_mult*(1-exp(-LAI)));
else
    rs = 0;
end

    % 2.8 a stability correction is needed for the surface to canopy level resistance
if  rs > 0
    aeta_g = 0;
    i = 0;
	rs_change = 1;
    while i < 100 && abs(rs_change) > 0.0001
        i=i+1;
% 1. estimate ustar and ua(z_displ_can) above surface from ras and z_0mg, z_0hg, and zref=z_displ_can
        [psi_m1,psi_h1] = StabilityFunctions(aeta_g);
        [psi_m2,psi_h2] = StabilityFunctions(aeta_g*z_0hg/z_displ_can);
        [psi_m3,psi_h3] = StabilityFunctions(aeta_g*z_0mg/z_displ_can);
        ustar_below1 = (1/rs)/karman * (log(z_displ_can/z_0hg) - psi_h1 + psi_h2);
        wind_z_displ_can = ustar_below1/karman * (log(z_displ_can/z_0mg) - psi_m1 + psi_m3);
        
% 2. estimate aeta above surface
        if SD > 0
            T_sup = T_snowsurf;
        else
            T_sup = T_air;
        end
        aeta_g = RichardsonToAeta(z_displ_can,T_can,T_can-T_sup,wind_z_displ_can,z_0mg,z_0hg);
% 3. new guess of ustar based on uadisplcan and new aeta_g
        [psi_m1,psi_h1] = StabilityFunctions(aeta_g);
        [psi_m2,psi_h2] = StabilityFunctions(aeta_g*z_0mg/z_displ_can);
        ustar_below2 = wind_z_displ_can*karman/(log((z_displ_can)/z_0mg) - psi_m1 + psi_m2);

% 4. TRANSFER COEFFICIENT FOR SCALARS below CANOPY
        [psi_m3,psi_h3] = StabilityFunctions(aeta_g*z_0hg/z_displ_can);
        ch = ustar_below2*karman/(log((z_displ_can)/z_0hg) - psi_h1 + psi_h3);

% 5. new guess for AERODYNAMIC RESISTANCE below CANOPY
        rs_change = 1/ch - rs;
        rs = 1/ch;
    end
end


%{
Surface resistance for transpiration (van den Hurk et al, 2000)
In case there is no soil data, use the air temperature as guess for the soil
temperature, and skip soil moisture function.
%}
    % multipl. incr. of canopy surface resistance as a function of downward SW (van den Burk et al., 2000)
a = 0.81; b = 0.004; c = 0.05;
f1 = (a*(1 + b*SW_in_ac))/(b*SW_in_ac + c );
f1 = max(f1,1);
    % -"- as a function of atm. vapor pressure deficit (van den Burk et al (2000)
P_sat = WaterSaturationPressure(T_air);
f3 = 1/exp(-f3_gd*((1-RelHum)*P_sat));
    % -"- as a function of soil temperature
F4_A = 1.75; F4_B = 0.5;
	if SD > 0
        Temp = 0;
    else
        Temp = T_air - 273.15;
    end
f4 = 1/(1 - exp(-F4_A*(max(0.00001,Temp)^F4_B)));
%{
if ( useSoilLayers )
    % get_f2f4 = -"- as a function of liquid water content and soil temperature in the root zone
rs_transp = rs_min*f1* get_f2f4(Xdata.SoilNode, &Xdata.Edata[0]) *f3/LAI;
else
%}
rs_transp = rs_min*f1*f4*f3/LAI;
%end

% exchange coefficients sensible heat
ch_canopy = rho_air*C_air/(ra + rc);

% latent heat interception
if T_air < 273.15
    ce_condensation  = 0.622*lh_subl/(gas_const_air*T_air*ra_incr_snow*(ra_e+rc)); %*max(0.1,frac_wet);
    ce_interception  = 0.622*lh_subl/(gas_const_air*T_air*ra_incr_snow*(ra_e+rc)); %*frac_wet;
    ce_transpiration = 0.0;
else
    ce_condensation  = 0.622*lh_vap/(gas_const_air*T_air*(ra_e+rc)); %*max(0.1,frac_wet);
    ce_interception  = 0.622*lh_vap/(gas_const_air*T_air*(ra_e+rc)); %*frac_wet;
    ce_transpiration = 0.622*lh_vap/(gas_const_air*T_air*(ra_e+rc+rs_transp)); %*(1.0-frac_wet);
end

ce_canopy = ce_interception*max(0.001,frac_wet) + ce_transpiration*(1-frac_wet);


end