tic

clear all
close all

%--------------------------------------------------------------------------
%--------------------  location-specific parameters  ----------------------
%{
CLM4.5 is recreated in this Toy Model to compare its parameterisations to
SNOWPACK using stand-scale forcing. Therefore, parameters for vegetation,
soil, etc. have to be prescribed for the particular location and cannot be
extracted from grid-scale forcing data as done for for global simulations.
These parameters include:
- position (latitude, longitude),
- canopy structure (PFT, LAI, SAI, ...), and
- soil type (soil colour; water content of uppermost soil layer is necessary
  but should be calculated).

Meteorological forcing data are necessary for the comparison of the
different parameterisations and an easy choice is Alptal(Switzerland), used
in SnowMIP2, a description of which is given by Rutter et al. (2009).

For the purpose of a more detailed comparison, the individual terms of the
vegetation energy balance are added to the output.
%}
Lat = 47+3/60;
Lon = 8+48/60;
Elevation = 1220;   % open area on same altitude as measurements above forest
PFT = 1;    % needleleaf evergreen forest
LAI = 2.5;  % effective LAI
SAI = 1;
CanCov = 0.96;
BasalArea = 0.004;
can_height = 25;
meteo_height = 35;
alb_soil_dir_open = 0.19;
alb_soil_dif_open = alb_soil_dir_open;
alb_soil_dir_forest = 0.11;
alb_soil_dif_forest = alb_soil_dir_forest;
frac_sand = 0.44;   % "silt" in Rutter et al. (2009)
frac_clay = 0.44;


%--------------------------------------------------------------------------
%-----------------------------  import data  ------------------------------

[Data,Variables] = xlsread('data_Alptal_SNOWMIP2_200304.xls',1);

% time
dates = Variables(4:end,1);
time = nan(length(dates(:,1)),1); date = char(length(time),15);
for t=1:length(time)
    date(t,1:15) = dates{t,:};
end
time(1) = datenum(str2double(date(1,3:6)),str2double(date(1,8:9)),str2double(date(1,11:12))...
    ,str2double(date(1,14:15)),0,0);
time(2) = datenum(str2double(date(2,3:6)),str2double(date(2,8:9)),str2double(date(2,11:12)),...
    str2double(date(2,14:15)),0,0);
dt = time(2)-time(1);
for t=3:length(time)
    time(t) = time(t-1)+dt;
end

dt_CLM = dt*24*3600;    % for CLM timestep required in [s]
dt_SP = dt*24;          % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1


%--------------------------------------------------------------------------
%-----------------------------  observations  -----------------------------

% forcing
SW_in_open = Data(:,2);       % downward shortwave radiation on the open meadow [W m^{-2}]
LW_in_open = Data(:,5);       % downward longwave radiation on the open meadow [W m^{-2}]
T_air_ac = Data(:,8);         % air temperature above the forest canopy [°C]
    T_air_ac = T_air_ac + 273.15;
Wind_ac = Data(:,9);          % wind speed above the forest canopy [m s^{-1}]
RelHum_open = Data(:,4);      % relative humidity on the open meadow [%]
Prec = Data(:,6);             % precipitation, not divided into rain and snow [mm h^{-1}]
    Prec = Prec/3600;         % [mm s^{-1}] required, [mm h^{-1}] measured
load GWL_Proxies_Alptal.mat
    Proxy_GWL = factor_GWL_0304;
z_snow_open = Data(:,7);      % snow depth on the open meadow [cm]
    z_snow_open = z_snow_open/100;
z_snow_forest = 0.37*z_snow_open;       % simple comparison between manual measurements
SnowAge = zeros(size(z_snow_forest));
for t=2:length(time)
    if z_snow_forest(t) > 0 && (z_snow_forest(t)-z_snow_forest(t-1)) <= 0
        SnowAge(t) = SnowAge(t-1) + dt;
    end
end
T_snowsurf1 = Data(:,14);   % snow temperature (10cm above soil) below the canopy
    T_snowsurf1 = T_snowsurf1 + 273.15;
T_snowsurf2 = Data(:,15);   % snow temperature (40cm above soil) below the canopy
    T_snowsurf2 = T_snowsurf2 + 273.15;

% validation
SW_in_ac_val = Data(:,16);      % downward shortwave radiation above the forest canopy [W m^{-2}]
SW_out_ac_val = Data(:,17);     % reflected shortwave radiation above the forest canopy [W m^{-2}]
SW_out_ac_full = Data(:,10);
LW_in_ac_val = Data(:,18);      % downward longwave radiation above the forest canopy [W m^{-2}]
LW_out_ac_val = Data(:,19);     % upward longwave radiation above the forest canopy [W m^{-2}]
SW_in_bc_val = Data(:,21);      % downward shortwave radiation below the forest canopy [W m^{-2}]
SW_out_bc_val = Data(:,22);     % reflected shortwave radiation below the forest canopy [W m^{-2}]
LW_in_bc_val = Data(:,23);      % downward longwave radiation below the forest canopy [W m^{-2}]
LW_out_bc_val = Data(:,24);     % upward longwave radiation below the forest canopy [W m^{-2}]

% differentiating between snow and rain
frac_prec_snow = nan(length(time),1);
forc_snow = nan(length(time),1); forc_rain = nan(length(time),1);
for t=1:length(time)
    frac_prec_snow(t) = RainSnowPartitioningAlptal(T_air_ac(t));
    forc_snow(t) = Prec(t)*frac_prec_snow(t);
    forc_rain(t) = Prec(t)*(1-frac_prec_snow(t));
end

% differentiating between direct/diffuse and VIS/NIR solar radiation
%{
Calculate a diffuse fraction depending on effective emissivity of the sky,
which should be a reliable measure of the cloudiness.
%}
ema_eff = nan(size(T_air_ac));
for t=1:length(time)
ema_eff(t) = LW_in_ac_val(t)/(5.67*10^(-8)*T_air_ac(t)^4);
end
frac_dif = (ema_eff-0.6)*2;
for t=1:length(time)
   if frac_dif(t) > 1
       frac_dif(t) = 1;
   elseif frac_dif(t) < 0
       frac_dif(t) = 0;
   end
   if SW_in_ac_val(t) == 0
       frac_dif(t) = 1;
   end
end
frac_dir = 1 - frac_dif;
forc_sol_dir = zeros(length(time),2); forc_sol_dif = zeros(length(time),2);
forc_sol_dir(:,1) = SW_in_ac_val.*frac_dir;
forc_sol_dif(:,1) = SW_in_ac_val.*frac_dif;
NIR_available = 0;

%{
CLM4.5 features a multiple-layer canopy but only for photosynthesis-related
calculations. This influences the energy balance only via the stomatal
resistance. By setting the canopy layers to 1, the "big leaf"
parameterisation is applied, which should be equal to the previous CLM4
version.
%}
CLM45_canopy_layers = 1;

%{
Variable to switch between CLM4.5 version with or without heat mass
parameterisation from SNOWPACK.
HM_in_CLM = 0   -> no Heat Mass
HM_in_CLM = 1   -> Heat Mass
%}
HM_in_CLM = 0;

%--------------------------------------------------------------------------
%----------------------------  initialisation  ----------------------------
% CLM4.5
IntStor = nan(size(time));
frac_wet = nan(size(time));
SW_net_veg_CLM = nan(size(time));
LW_in_bc_CLM = nan(size(time));
LW_out_ac_CLM = nan(size(time));
T_veg_CLM = nan(size(time));
T_air_wc_CLM = nan(size(time));     % air temperature within canopy space
T_ref2m_CLM = nan(size(time));      % 2 m height surface air temperature
        % CLM4.5 vegetation energy balance terms (TV_old from previous iteration step)
EB_CLM = nan(length(time),17,41);
%{
Terms of Energy Balance matrix:
1)  net direct SW radiation
2)  net diffuse SW radiation
    1) + 2) = net SW radiation
3)  net LW radiation from atmosphere (gain) -> emv*forc_lwrad + (1-emv)*(1-emg)*emv*forc_lwrad
4)  net LW radiation from vegetation (loss) -> -2*emv*sb*TV_old^4 + emv*(1-emg)*emv*sb*TV_old^4
5)  net LW radiation from second vegetation layer (gain) -> 0 for CLM
6)  net LW radiation from ground (gain)     -> emv*emg*sb*Tgrnd^4
7)  net interaction with second vegetation layer -> 0 for CLM
    3) + 4) + [5) - 7)] + 6) = net LW radiation
    1) + 2) + 3) + 4) + [5) - 7)] + 6) + 7) = net radiation
    '-> complicated because of SNOWPACK
8)  net conductive heat flux -> 0 for CLM
9)  net sensible heat flux
10) net latent heat flux
11) dSW / dTveg -> per definition 0
12) dLW / dTveg -> derived from 4)
13) dTT / dTveg -> interaction with second vegetation layer 0 for CLM
14) dHM / dTveg -> conductive heat flux 0 for CLM
15) dSH / dTveg
16) dLH / dTveg
%}

% SNOWPACK
SW_net_can_SP = nan(size(time));
SW_net_trunk_SP = nan(size(time));
LW_in_bc_SP = nan(size(time));
LW_out_ac_SP = nan(size(time));
T_can_SP = nan(size(time));
T_trunk_SP = nan(size(time));
EB_SP_leaf = nan(length(time),16,7);
EB_SP_trunk = nan(length(time),16,7);
%{
Terms of Energy Balance matrix:
1)  net direct SW radiation
2)  net diffuse SW radiation
    1) + 2) = net SW radiation
3)  net LW radiation from atmosphere (gain)
4)  net LW radiation from vegetation (loss) -> emitted by respective layer
5)  net LW radiation from second vegetation layer (gain) -> absorbed from second layer
6)  net LW radiation from ground (gain)
7)  net interaction with second vegetation layer -> due to theory and
    derivation 7) already included in 5), therefore this:
    3) + 4) + [5) - 7)] + 6) = net LW radiation
    1) + 2) + 3) + 4) + [5) - 7)] + 6) + 7) = net radiation
8)  net conductive heat flux
9)  net sensible heat flux
10) net latent heat flux
11) dSW / dTveg -> per definition 0
12) dLW / dTveg -> derived from 4)
13) dTT / dTveg -> 0 for trunk layer because EB derived from this term for leaf layer
14) dHM / dTveg
15) dSH / dTveg
16) dLH / dTveg
%}

% initial conditions
IntStor(4206) = 0.1007;        % value taken from simulation using open field forcing
frac_wet(4206) = 0.4231;       % see -^
T_veg_CLM(4206) = T_air_ac(4206);
T_can_SP(4206) = T_air_ac(4206);
T_trunk_SP(4206) = T_snowsurf2(4206);
%{
Initialising the vegetation temperature is a tough one. Since t=1 is at
midnight (thus no solar radiation) it might be ok to assume the vegetation
temperature, which should be designed similarly to the canopy temperature
in SNOWPACK, equals the temperature above the canopy (10m above) and the
trunk temperature equals the temperature at 40cm above the ground
(T_snowsurf2). Furthermore, the effect of the canopy heat flux, i.e. the
thermal inertia effect of storing absorbed solar radiation during daytime,
is not considered in CLM4.5 anyway. However, for SNOWPACK it would dismiss
the effect of the canopy heat flux.
%}

%--------------------------------------------------------------------------
%-----------------------------  computation  ------------------------------
% roughness length and z0 or bare soil prescribed within Canopy_SNOWPACK2L.m for now
rl = nan;
z0_bs = nan;

alb_snow_fix = 0.8;

for t=4207:length(time)
    if z_snow_forest(t) > 0
        frac_sno = 1;
    else
        frac_sno = 0;
    end
    alb_ground = (1-frac_sno)*alb_soil_dir_forest + frac_sno*...
        ((alb_snow_fix - 0.3)*exp(-SnowAge(t)/7) + 0.3);
    
    date = datestr(time(t),'yyyy-mm-dd-HH-MM-SS');
    JulDay = datenum(0000,str2double(date(6:7)),str2double(date(9:10)),...
        str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));

    [IntStor(t),frac_wet(t),SW_net_veg_CLM(t),LW_in_bc_CLM(t),LW_out_ac_CLM(t),...
        T_veg_CLM(t),T_air_wc_CLM(t),T_ref2m_CLM(t),EB_CLM(t,:,:),EB_SP_leaf(t,:,:),EB_SP_trunk(t,:,:),...
        SW_net_can_SP(t),SW_net_trunk_SP(t),LW_in_bc_SP(t),LW_out_ac_SP(t),...
        T_can_SP(t),T_trunk_SP(t)] = ...
    SNOWPACKwithinCLM45driver(Lat,Lon,Elevation,JulDay,dt_CLM,dt_SP,...
        LAI,SAI,can_height,BasalArea,meteo_height,CLM45_canopy_layers,HM_in_CLM,...
        squeeze(forc_sol_dir(t,:)),squeeze(forc_sol_dif(t,:)),NIR_available,LW_in_ac_val(t),T_air_ac(t),...
        Wind_ac(t),RelHum_open(t),forc_rain(t),forc_snow(t),z_snow_forest(t),frac_sno,alb_ground,...
        T_snowsurf2(t),T_veg_CLM(t-1),T_can_SP(t-1),T_trunk_SP(t-1),...
        frac_sand,frac_clay,1-frac_sand-frac_clay,Proxy_GWL(t),frac_wet(t-1),IntStor(t-1));
end

toc