function [h2ocan,fwet,sabv,dlrad,ulrad,t_veg,taf,t_ref2m,...
    EB_CLM,EB_SP_leaf,EB_SP_trunk,...
    NR_leaf,NR_trunk,LW_in_bc,LW_refl_ac,t_can,t_trunk] = ...
    SNOWPACKwithinCLM45driver(lat,lon,altitude,JulDay,dt_CLM,dt_SP,...
    tlai,tsai,htop,BasalArea,forc_hgt,nlevcan,HM_in_CLM,...
    forc_sol_dir,forc_sol_dif,NIR_logical,forc_lwrad,forc_t,...
    forc_wind,forc_rh,forc_rain,forc_snow,snowdp,frac_sno,alb_ground,...
    t_soisno,t_veg,t_can,t_trunk,sand,clay,om,ProxySoilWater,fwet,h2ocan)
%{
Main function running the CLM4.5 replica module and calling the individual
functions. Similar to clmdriver.f90 within original model code.
Chronological order in CLM:
1) Hydrology1Mod.f90        -> CanopyHydrology45.m (forcing and time step in [s]!!!)
2) SurfaceRadiationMod.f90  -> SurfaceRadiation.m
3) Biogeophysics1Mod.f90    -> Biogeophysics1.m & surrounding parameter calculations
4) CanopyFluxesMod.f90      -> CanopyFluxes.m
5) Biogeophysics2Mod.f90    -> ground/soil/snow calculations, not desired here
6) Hydrology2Mod.f90        -> only soil & snow hydrology - not relevant here
7) SurfaceAlbedoMod.f90     -> TwoStream.m (& SoilAlbedo.m)
Step 7) calculates the surface albedos (incl. two-stream approximation) for
the next timestep. Since this would create additional (difficult)
calculations prior to the first timestep, the calculation of surface albedos
is shifted to the start of the module 0).
    
Modified version to include module for energy fluxes from SNOWPACK.
Included parallel to CanopyFluxes.m so that actual energy balance of
vegetation is calculated twice, by each of the modules, using the same
input (forcing, roughness parameters, canopy coverage, LAI/SAI, results of
mass fluxes, etc.).
%}

%--------------------------------------------------------------------------
%--------------------------  Physical Constants  --------------------------

grav = 9.80616; % acceleration of gravity [m s^{-2}]
rair = 6.02214*10^26 * 1.38065*10^(-23) / 28.966; % dry air gas constant [J K^{-1} kg^{-1}]
cpair = 1.00464*10^3; % specific heat of dry air [J K^{-1} kg^{-1}]

zlnd = 0.01;    % roughness length for soil [m] (tunable)
zsno = 0.0024;  % roughness length for snow [m] (tunable)
z0mr = 0.055;   % ration of momentum roughness length to canopy top height (same for all boreal trees)
displar = 0.67; % ration of displacement height to canopy top height (same for all boreal trees)


%--------------------------------------------------------------------------
%-------------------------------  Forcing  --------------------------------

%{
Forcing data in CLM downscaled depending on elevation difference between
height of atmospheric bottom layer (input elevation) and height of land
grid cell (actual calculation elevation). Since this module is used for
observational sites this is redundant.
%}

% solar elevation angle
    % CLM
[coszen,decl] = CosZenCalculation(JulDay,lat,lon);
    % SNOWPACK
sinelev = SolarElevationAngle(JulDay,lat,lon);
SolAngle = asin(sinelev);

% solar forcing data include nir radiation?
waveband = 1 + NIR_logical; % NIR_logical 1 or 0 depending on available nir forcing

% atmospheric pressure [Pa]
dalr = 0.0065;      % dry adiabatic lapse rate [K m^{-1}]
p0 = 101325.;       % standard pressure at sea level [Pa]
T0 = 288.15;        % standard temperature (25°C) at sea level [K]
R0 = 6356766.0;     % radius of Earth [m]
GasR_da = 287.0423; % gas constant for dry air [J kg^{-1} K^{-1}]
forc_pbot = p0*(1 - (dalr*R0*altitude)/(T0*(R0+altitude)))^(grav/(dalr*GasR_da));

% atmospheric potential temperature [K]
forc_th = forc_t*exp((forc_hgt/(rair*forc_t/grav))*(rair/cpair));

% atmospheric O2 and CO2 concentrations [Pa] -  CO2 calculation neither prognostic nor diagnostic
forc_po2 = 0.209*forc_pbot;          % per constant atmospheric O2 molar ratio [mol mol^{-1}]
forc_pco2 = 355*10^(-6)*forc_pbot;   % per atmospheric CO2 molar ratio (by volume) [umol mol^{-1}]

% atmospheric specific humidity [kg kg^{-1}]
tdc = min(50,max(-50,forc_t-273.15));
if forc_t > 273.15      % saturation vapor pressure over water [Pa]
    e = 100*(6.107799961 + tdc*(4.436518521*10^(-1) + tdc*(1.428945805*10^(-2)...
        + tdc*(2.650648471*10^(-4) + tdc*(3.031240396*10^(-6)...
        + tdc*(2.034080948*10^(-8) + tdc*6.136820929*10^(-11)))))));
else                    % saturation vapor pressure over ice [Pa]
    e = 100*(6.109177956 + tdc*(5.034698970*10^(-1) + tdc*(1.886013408*10^(-2)...
        + tdc*(4.176223716*10^(-4) + tdc*(5.824720280*10^(-6)...
        + tdc*(4.838803174*10^(-8) + tdc*1.838826904*10^(-10)))))));
end
qsat = 0.622*e/(forc_pbot-0.378*e);
forc_q = qsat*forc_rh/100;

% density of air [kg m^{-3}]

forc_vp = forc_q*forc_pbot/(0.622 + 0.378*forc_q);      % atmospheric vapor pressure [Pa]
forc_rho = (forc_pbot - 0.378*forc_vp)/(rair*forc_t);   % density of most air [kg m^{-3}]

% forcing heights
forc_hgt_u = forc_hgt;
forc_hgt_t = forc_hgt;
forc_hgt_q = forc_hgt;


%--------------------------------------------------------------------------
%-------------------------  Physical Parameters  --------------------------

% vertical fraction of vegetation covered by snow
%{
Adjust lai and sai for burying by snow. If exposed lai and sai are less than
0.05, set equal to zero to prevent numerical problems associated with very
small lai and sai. Snow burial fraction for short vegetation (e.g. grasses)
as in Wang and Zeng, 2007.
%}
% ol = min(max(snowdp-hbot,0),htop-hbot);
% fb = 1 - ol/max(10^(-6),htop-hbot);
fb = 1; % hbot for boreal trees at least 7m and snow depth not expected to be that high
elai = max(tlai*fb,0);
if elai < 0.05
    elai = 0;
end
esai = max(tsai*fb,0);
if esai < 0.05
    esai = 0;
end
if elai+esai >= 0.05
    frac_veg_nosno_alb = 1;     % fraction of vegetation free of snow
else
    frac_veg_nosno_alb = 0;
end

% surface emissivities
avmuir = 1;  % infrared inverse optical depth for longwave radiation
emv = 1 - exp(-(elai+esai)/avmuir);
emsoil = 0.96;
emsnow = 0.97;
emg = emsoil*(1-frac_sno) + emsnow*frac_sno;

albgrd = zeros(1,2);
albgrd(1) = alb_ground;
albgri = albgrd;

%--------------------------------------------------------------------------
%-------------------  Computations = Calling Functions  -------------------

%---------------------  0) albedos for this timestep  ---------------------
[fabd,fabi,albd,albi,ftdd,ftid,ftii,omega,gdir,...
    fabd_sun_z,fabd_sha_z,fabi_sun_z,fabi_sha_z,frac_sun_z,nrad,tlai_z,vcmaxcintsun,vcmaxcintsha]...
    = TwoStream_AddOnCLM45(coszen,1,elai,esai,t_veg,fwet,albgrd,albgri,waveband,nlevcan);
                                % 1 -> PFT: boreal ENF

% 1) canopy interception and precipitation on ground (incl. wet fraction)
[IntCapacity,h2ocan,qflx_prec_intr,qflx_through_snow,qflx_through_rain,qflx_candrip_snow,...
    qflx_candrip_rain,qflx_prec_grnd_snow,qflx_prec_grnd_rain,fwet]...
    = CanopyHydrology45(elai,esai,frac_veg_nosno_alb,forc_rain,forc_snow,h2ocan,dt_CLM);


%----------------------  2) surface solar radiation  ----------------------
[sabv,sabvd,sabvi,sabg,fsa,laisun,laisha,parsun_z,parsha_z,laisun_z,laisha_z]...
    = SurfaceRadiation(elai,esai,gdir,coszen,forc_sol_dir,forc_sol_dif,waveband,...
    albgrd,albgri,omega,fabd,fabi,ftdd,ftid,ftii,...
    nrad,tlai_z,frac_sun_z,fabd_sun_z,fabd_sha_z,fabi_sun_z,fabi_sha_z);


%---------------  3) leaf temperature and surface fluxes ?  ---------------
%{
Difficult since several input parameters for soil/ground are necessary and
not available from observations. Therefore, reliable calculations are done
here while the calculations for soil parameters are done in Biogeophysics1,
which would otherwise include all calculations in 3).
Several variables are calculated from soil characteristics such as
fractions of sand, clay or organic matter, but as those are definitely not
expected to be available constant values are prescribed.
Furthermore, the volumetric fraction of water and ice within the soil is
required, which is then replaced by a proxy.
%}
t_grnd = t_soisno;  % calculated out of snow cover fraction and wet fraction in CLM4.5

if isnan(sand)
    sand = 0.3; % boreal value taken from input file for CLM
end
if isnan(clay)
    clay = 0.25; % boreal value taken from input file for CLM
end
if isnan(om)
    om = 0.45;   % ? residual?
end
% saturated vapor pressure, specific humidity
[qg,soilbeta,btran,btran2] = ...
    Biogeophysics1(ProxySoilWater,om,sand,clay,frac_sno,t_grnd,forc_pbot,forc_q);

% ground roughness lengths
if frac_sno > 0
    z0mg = zsno;
else
    z0mg = zlnd;
end
z0hg = z0mg;
z0qg = z0mg;

% potential, virtual temperature, and wind speed at the reference height
beta = 1;
zii = 1000;
thv = forc_th*(1 + 0.61*forc_q);

% roughness lengths over vegetation
z0m = z0mr*htop;
displa = displar*htop;
z0mv = z0m;
z0hv = z0mv;
z0qv = z0mv;
% make forcing height a pft-level quantity that is the atmospheric forcing height plus pft's z0m+displa
if frac_veg_nosno_alb == 0
    forc_hgt_u_pft = forc_hgt_u + z0mg + displa;
    forc_hgt_t_pft = forc_hgt_t + z0mg + displa;
    forc_hgt_q_pft = forc_hgt_q + z0mg + displa;
else
    forc_hgt_u_pft = forc_hgt_u + z0m + displa;
    forc_hgt_t_pft = forc_hgt_t + z0m + displa;
    forc_hgt_q_pft = forc_hgt_q + z0m + displa;
end

% intermediate variable for atmospheric potential temperature [K]
thm = forc_t + 0.0098*forc_hgt_t_pft;


%-----  4) leaf temperature and surface fluxes for vegetated patches  -----
% calculation of soil wetness factor moved from CanopyFluxes to Biogeophysics1
if HM_in_CLM == 0
[dlrad,ulrad,h2ocan,t_veg,taf,t_ref2m,EB_CLM]...
    = CanopyFluxes(lat,decl,dt_CLM,elai,esai,htop,displa,z0mv,z0mg,forc_hgt_u_pft,forc_hgt_t_pft,...
    emv,emg,forc_lwrad,thm,thv,forc_th,forc_q,forc_wind,...
    forc_pbot,forc_rho,forc_po2,forc_pco2,t_grnd,qg,t_veg,sabv,sabvd,sabvi,fwet,h2ocan,snowdp,soilbeta,btran,...
    frac_veg_nosno_alb,nrad,tlai_z,vcmaxcintsun,vcmaxcintsha,parsun_z,parsha_z,laisun_z,laisha_z,laisun,laisha);
elseif HM_in_CLM == 1
[dlrad,ulrad,h2ocan,t_veg,taf,t_ref2m,EB_CLM]...
    = CanopyFluxes_HeatMass(lat,decl,dt_CLM,elai,esai,htop,displa,z0mv,z0mg,forc_hgt_u_pft,forc_hgt_t_pft,...
    emv,emg,forc_lwrad,thm,thv,forc_th,forc_q,forc_wind,...
    forc_pbot,forc_rho,forc_po2,forc_pco2,t_grnd,qg,t_veg,sabv,sabvd,sabvi,fwet,h2ocan,snowdp,soilbeta,btran,...
    frac_veg_nosno_alb,nrad,tlai_z,vcmaxcintsun,vcmaxcintsha,parsun_z,parsha_z,laisun_z,laisha_z,laisun,laisha);
end


[t_can,t_trunk,NR_leaf,NR_trunk,RN_can,H_can,LE_can,IntStor,frac_wet,...
    SW_refl_ac,LW_refl_ac,SW_in_bc,LW_in_bc,SW_refl_bc,LW_refl_bc,SW_net_trunk,LW_net_trunk,em_air_bc,alb_tot,EB_SP_leaf,EB_SP_trunk] ...
    = SNOWPACK2L_EnergyFluxes_withinCLM45(htop,elai,1-emv,BasalArea,SolAngle,...    %= SNOWPACK2L_EnergyFluxes_withinCLM45(htop,elai+esai,1-emv,BasalArea,SolAngle,...
    forc_sol_dir(1)+forc_sol_dif(1),forc_sol_dif(1),forc_lwrad,thm,forc_rh/100,forc_wind,...
    snowdp,t_grnd,emg,albgrd(1),forc_hgt_u,zsno,zlnd,t_can,t_trunk,fwet,h2ocan,IntCapacity,dt_SP);

%-----  5) soil/snow & ground temperature and update surface fluxes  ------
%{
Outgoing LW radiation (and net LW radiation for ground) calculated based on
updated soil/ground temperature. However, soil temperature calculation too
complex and requiring too many parameters for this purpose. (And not
necessary anyway.)
%}
end