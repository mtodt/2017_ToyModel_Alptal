function [T_can,T_trunk,NR_leaf,NR_trunk,RN_can,H_can,LE_can,IntStor,frac_wet,...
    SW_refl_ac,LW_refl_ac,SW_in_bc,LW_in_bc,SW_refl_bc,LW_refl_bc,SW_net_trunk,LW_net_trunk,em_air_bc,alb_tot,EB_leaf,EB_trunk] ...
    = SNOWPACK2L_EnergyFluxes_withinCLM45(z_can,LAI,frac_through,BasalArea,SolAngle,SW_in_ac,SW_in_ac_dif,LW_in_ac,T_air,RelHum,wind,...
    SnowDepth,T_ground,emg,alb_ground,z_wind,roughness_length,z0_BareSoil,T_can,T_trunk,frac_wet,IntStor,IntCapacity,dt)
%{
Canopy surface energy balance (net radiation, sensible and latent heat fluxes)
and final mass balance (evaporation of intercepted water, and transpiration.    
%}
%--------------------------------------------------------------------------
%-------------------------  Physical Parameters  --------------------------
%--------------------------------------------------------------------------

T_freez = 273.15;       % freezing temperature [K]
T_melt = 273.15;        % melting temperature [K]

lh_subl = 2.838*10^6;   % latent heat of sublimation [J kg^{-1}]
lh_vap = 2.504*10^6;    % latent heat of vaporisation [J kg^{-1}]

boltz = 5.67*10^(-8);   % Stefan-Boltzmann constant [W m^{-2} K^{-4}]

%--------------------------  Parameters for EB  ---------------------------

alb_can_dry = 0.11;     % albedo of dry canopy (calibr: 0.09, Alptal)
alb_can_wet = 0.11;     % albedo of wet canopy (calibr: 0.09, Alptal)
alb_can_snow = 0.35;    % albedo of snow covered canopy (calibr. Alptal, but 0.3 in Gouttevin et al. (2015))
alb_trunk = 0.09; 		% trunk albedo
em_leaf = 1;            % canopy emissivity
em_trunk = 1; 			% trunk emissivity

frac_LAI_top = 0.5; 	% fraction of LAI attributed to uppermost layer (default, calibrated for Alptal)
k_LAI = 0.75;           % radiation transmissivity parameter [0.4-0.8] (calibrated on Alptal)

BM_heat_cap = 2800;     % [J kg^{-1} K^{-1}]
BM_density = 900;       % [kg m^{-3}]
Diameter = 1;           % average canopy (tree) diameter [m]
%BasalArea = 0.004;      % Canopy Basal Area [m2/m2], value for closed canopies like Alptal


dT_can_maxperh = 7;     % maximum allowed canopy temperature change [K h^{-1}]

%-----------------  scales for aerodynamical resistances  -----------------

% canopy roughness lengths for heat and momentum (in DataClasses.cc)
z_0m = z_can*0.1;
z_0h = z_0m*0.1;        % ??? z_0h/z_0m = 0.999 in Gouttevin et al. (2015)

% displacement height defining canopy level (in DataClasses.cc)
z_displ = z_can*(2/3);


%--------------------------------------------------------------------------
%-----------------------------  Calculations  -----------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%----------------  2.1 prepare for canopy energy balance  -----------------

% Radiation Transmissivity (Beer-Lambert type of law)
    % first, transmissivity of diffuse (and longwave) radiation
eps = 10^(-3);
if frac_LAI_top < eps || 1-frac_LAI_top < eps
    frac_LAI_top = 1;
end
PAI = 0; % PAI [additional plant area index] could be a function of IntStor
sig_leaf = 1 - exp(-k_LAI*(frac_LAI_top*LAI+PAI)/sin(pi/2));
sig_trunk = 1 - exp(-k_LAI*((1-frac_LAI_top)*LAI+PAI)/sin(pi/2));

    % second, transmissivity of direct solar radiation
sig_leaf_dir = 1 - exp(-k_LAI*(frac_LAI_top*LAI+PAI)/max(sin(SolAngle),0.0001));
sig_trunk_dir = 1 - exp(-k_LAI*((1-frac_LAI_top)*LAI+PAI)/max(sin(SolAngle),0.0001));


% Reference Height [m above snow surface] for meteo input, min. 2 m above z_can above snow surface
ref_height = z_wind;    % (alpine3d)? Xdata.Cdata.height+height_of_wind_val : height_of_wind_val;
z_ref = max(2 + (z_can - SnowDepth),ref_height - SnowDepth);
if SnowDepth > 0.03
    z_0mg = roughness_length;
    z_0 = roughness_length;
else
    z_0mg = z0_BareSoil;
    z_0 = z0_BareSoil;
end

% Turbulent Transport Coefficients
[ch_canopy,ce_canopy,ce_transpiration,ce_interception,ce_condensation,ra,rs]...
    = TurbulentExchange(SW_in_ac,T_air,wind,RelHum,T_can,LAI,z_can,z_ref,z_0mg,frac_wet,SnowDepth,T_ground);

%--------------------------------------------------------------------------
%-------------------  2.2 Energy balance of the canopy  -------------------
%{
The main purpose is to estimate the evapotranspiration loss, and the radiation
balance of the canopy, which influences the snowpack below and the reflection/
emittance to the atmosphere.
Method:
The energy balance of the canopy is assumed to be equal to
                (1)   RnCanopy = HCanopy + LECanopy,
where RnCanopy is a function of
a) incoming shortwave and longwave radiation,
b) albedo, transmissivity, emissivity and temperature of the canopy, and
c) albedo, emissivity and temperature of the ground below, taking multiple
reflection into account (see documentation).
Sensible and latent heat fluxes are computed using the standard bulk
formulations, assuming an logarithmic wind profile above the canopy. The
stomatal control of transpiration is represented by an additional surface
resistance, estimated as a function of incoming solar radiation, atmospheric
vapour pressure deficit and soil water content.
Numerical solution  following SiSPAT manual (Isabelle Braud, 2000):
Equation (1) is linearised around the canopy temperature at time t using the
temperature from the previous timestep t-1, so that:
                (2) RnCanopy  = r0  + r1  * TempCanopy(t)
                (3) HCanopy   = h0  + h1  * TempCanopy(t)
                (4) LECanopy  = le0 + le1 * TempCanopy(t)
TempCanopy(t) is given by inserting eq (2)-(4) into (1). See the functions
for r0,r1, h0,h1, le0, and le1 for details
Alternative (to be implemented shortly):
                (1) RnCanopy = HCanopy + LECanopy + DQCanopy,
where DQCanopy = change of heat content of Canopy
                (5) DQCanopy = dq0 + dq1 * TempCanopy(t)
and dq/dt = HeatCapacity * (TempCanopy(t)-TempCanopy(t-1) with
HeatCapacity = HeatCapacityDry + IntStorage(t-1) * L
INPUT: emissivity eground, surface temperature TGROUND, and albedo ground ALBGROUND
Start of Energy Balance Loop, 3 iterations should be enough in most cases
%}
%{
% local energy flux variables
RNCANOPY = 0; HCANOPY = 0; LECANOPY = 0; LECANOPYCORR = 0;
% double iswrac, rswrac, iswrbc, rswrbc, ilwrac, rlwrac, ilwrbc, rlwrbc, rsnet=IOUtils::nodata;

% local auxiliary variables
double canopyalb=IOUtils::nodata;
double h0, h1, ht0, ht1, le0, le1, let0, let1, HM0, HM1, HMt0, HMt1, TT0, TT1;
double r0, r1, r2, rt0, rt1, rt2;

double canopyclosuredirect=IOUtils::nodata, radfracdirect=IOUtils::nodata, r1p, r2p;
double CanopyEvaporation=0., INTEVAP=0., TRANSPIRATION=0.;
%}

T_can_pre = T_can;
T_trunk_pre = T_trunk;

% calculate Canopy Heat Mass based on canopy basal area and LAI
%{
For leaves: mean leaf thickness assumed to be 1 mm = 0.001 m.
For trunks: trunk volume calculated based on basal area, height and conic
form.
%}
HM_leaf = 0.001*LAI*BM_density*BM_heat_cap;
HM_trunk= 0.5*BasalArea*z_can*BM_density*BM_heat_cap;

i=0;
% Energy Balance output: size stems from SNOWPACK + err  x  iterations
EB_leaf = nan(16,7);
EB_trunk = nan(16,1);

while i < 7
    i=i+1;  % shifted to the start for EB output
    
    T_can_old = T_can;
    T_trunk_old = T_trunk;
    
% update ce_canopy as function of wetfraction
    ce_canopy = max(0.001*ce_interception,ce_interception*frac_wet + ce_transpiration*(1-frac_wet));
    ce_condensation = ce_interception*max(0.1,frac_wet);
    
% canopy albedo - albedo of partly "wet" canopy = weighted average of dry and wet parts
    if T_air > T_melt
        alb_can = frac_wet*alb_can_wet + (1-frac_wet)*alb_can_dry;
    else
        alb_can = frac_wet*alb_can_snow + (1-frac_wet)*alb_can_dry;
    end
    alb_leaf = alb_can;
    
% compute properties r0 and r1 in eq (2) (and downward LW and SW for snowpack model)
    [r0,r1,r2,rt0,rt1,rt2,r1p,r2p,NR_leaf,NR_trunk,CC_dir,frac_dir,...
    r0_SWdir,r0_SWdif,r0_LW,r0_LW_atm,r0_LW_leaf,r0_LW_trunk,r0_LW_ground,r1_LW,r2_LW,...
    rt0_SWdir,rt0_SWdif,rt0_LW,rt0_LW_atm,rt0_LW_trunk,rt0_LW_leaf,rt0_LW_ground,rt1_LW,rt2_LW]...
        = NetRadiation2L(SW_in_ac,SW_in_ac_dif,SolAngle,LW_in_ac,T_air,T_ground,emg,T_can,T_trunk,z_can,LAI,frac_through,frac_wet,alb_ground);
    
% compute properties h0 and h1 in eq (3)
%{
sensible heat flux already a linear function of TC(t)
NOTE: for sparse canopies turbulent fluxes should be scaled in the canopy EB
calculation; for the moment scalingfactor is sigf*(1-direct_throughfall)
%}
    scalingfactor_leaf = sig_leaf*(1-frac_through);
    h1 = scalingfactor_leaf*ch_canopy;
    h0 = -scalingfactor_leaf*ch_canopy*T_air;
    scalingfactor_trunk = sig_trunk*(1-frac_through);
    ht1 = scalingfactor_trunk*ch_canopy;
    ht0 = -scalingfactor_trunk*ch_canopy*T_air;

% compute properties le0 and le1 in eq (4)
%{
LE = ce_canopy*(esat(TC(t)) - eair) is linearised around TC(t) by applying
esat(TC(t)) = esat(TC(t - 1)) + Desat(TC(t - 1)) / DT * (TC(t) - TC(t - 1))
NOTE: for the moment trunks do not exchange latent heat (no interception, no transpiration)
%}
    P_sat = WaterSaturationPressure(T_can);
    if T_can > 273.15
        dpdt = DSaturationPressureDT(lh_vap,T_can);
        le1 = scalingfactor_leaf*ce_canopy*dpdt;
    else
        dpdt = DSaturationPressureDT(lh_subl,T_can);
        le1 = scalingfactor_leaf*ce_canopy*dpdt;
    end
    P_sat_air = WaterSaturationPressure(T_air);
    vp_air = P_sat_air*RelHum;
    le0 = scalingfactor_leaf*ce_canopy *(P_sat-vp_air) - le1*T_can;
    % for the moment trunks do not exchange latent heat (no interc., no transp.)
    let1= 0;
    let0 = 0;

% the conductive heat flux to the leaf layer is already a linear function of TC(t)
    scalingfactor = 1;
    HM0 = -1*scalingfactor*HM_leaf/(dt*3600)*T_can_pre;
    HM1 =  scalingfactor*HM_leaf/(dt*3600);
    HMt0 = -1*scalingfactor*HM_trunk/(dt*3600)*T_trunk_pre;
    HMt1 =  scalingfactor*HM_trunk/(dt*3600);
    
%{
final canopy energy balance
Objective: Derive an analytical expression for Ttrunk = f(TC) from the
Trunk Energy Balance equation, then solve easily the Canopy Energy Balance.
Trunk energy balance equation:
rt0 + rt1 * Ttrunk + rt2 * TC = h0t + h1t *Ttrunk + let0 + let1 * Ttrunk + HMt0 + HMt1 * Ttrunk
=>  Ttrunk = (ht0 + let0 + HMt0 - rt0 -  r2t *TC) / (r1t - ht1 -let1 -HMt1)
rewritten as :  Ttrunk = TT0/r2 + TT1/r2 * TC
so that :       netradcanopy = r0 + r1 * TC + TT0 + TT1 * TC,
to be solved for TC.
%}
    TT0 = r2*(ht0 + let0 + HMt0 - rt0)/(rt1 - ht1 - let1 - HMt1);
    TT1 = -r2*rt2/(rt1 - ht1 - let1 - HMt1);
    T_trunk = (ht0 + let0 + HMt0 - rt0)/(rt1 - ht1 - let1 - HMt1)...
        - rt2/(rt1 - ht1 - let1 - HMt1)*T_can;

% Canopy Energy Balance
    % maximum allowed change of canopy temperature per hour
    T_can_change = (h0 + le0 - r0 + HM0 - TT0)/(r1 - h1 - le1 - HM1 + TT1) - T_can;

    % minimize the rate of change to CANOPYTEMP_MAXCHANGE_PERHOUR [K h^{-1}]
    T_can = T_can + T_can_change;
    T_trunk = T_can*TT1/r2 + TT0/r2 ;

    % re-compute Rn, H, and LE
    RN_can = r0 + r1*T_can + r2*T_trunk;
    H_can = h0 + h1*T_can;
    LE_can = le0 + le1*T_can;
        
    % added for EB output
    le0_out = le0;
    le1_out = le1;
        
    % re-compute in case of condensation/sublimation on canopy
    if LE_can < 0
        T_can = T_can - T_can_change;
        T_can_change = (h0 + le0*ce_condensation/ce_canopy - r0 + HM0 - TT0)/...
            (r1 - h1 - le1*ce_condensation/ce_canopy - HM1 + TT1) - T_can;
        T_can = T_can + T_can_change;

        T_trunk = T_can*TT1/r2 + TT0/r2 ;

        RN_can = r0 +  r1*T_can + r2*T_trunk;
        H_can = h0 + h1*T_can;
        LE_can = le0 + le1*T_can;
        
        % added for EB output
        le0_out = le0*ce_condensation/ce_canopy;
        le1_out = le1*ce_condensation/ce_canopy;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adding individual terms of energy balance to output - LEAF
SWnet_leaf = r0_SWdir + r0_SWdif;
SWnet_leaf_dir = r0_SWdir;
SWnet_leaf_dif = r0_SWdif;
LWnet_leaf = r0_LW + r1_LW*T_can_old;    % = LW_atm + LW_loss + (LW_gain-TTnet) + LW_ground
    LW_leaf_atm = r0_LW_atm;
    LW_leaf_loss = r0_LW_leaf + r1_LW*T_can_old;
    LW_leaf_gain = r0_LW_trunk + r2_LW*T_trunk_old;
    LW_leaf_ground = r0_LW_ground;
TTnet_leaf = TT0 + TT1*T_can_old;        % == r2_LW*T_trunk(_old)
Rnet_leaf = SWnet_leaf + LWnet_leaf + TTnet_leaf;
HMnet_leaf = HM0 + HM1*T_can_old;
SHnet_leaf = h0 + h1*T_can_old;
LHnet_leaf = le0_out + le1_out*T_can_old;
dLWnetdT_leaf = r1_LW;
dTTnetdT_leaf = TT1;
dHMnetdT_leaf = HM1;
dSHnetdT_leaf = h1;
dLHnetdT_leaf = le1_out;
EB_leaf(1,i) = SWnet_leaf_dir;
EB_leaf(2,i) = SWnet_leaf_dif;
EB_leaf(3,i) = LW_leaf_atm;
EB_leaf(4,i) = LW_leaf_loss;
EB_leaf(5,i) = LW_leaf_gain;
EB_leaf(6,i) = LW_leaf_ground;
EB_leaf(7,i) = TTnet_leaf;
EB_leaf(8,i) = HMnet_leaf;
EB_leaf(9,i) = SHnet_leaf;
EB_leaf(10,i) = LHnet_leaf;
EB_leaf(11,i) = 0;    % dSWnetdT_leaf
EB_leaf(12,i) = dLWnetdT_leaf;
EB_leaf(13,i) = dTTnetdT_leaf;
EB_leaf(14,i) = dHMnetdT_leaf;
EB_leaf(15,i) = dSHnetdT_leaf;
EB_leaf(16,i) = dLHnetdT_leaf;

% adding individual terms of energy balance to output - TRUNK
SWnet_trunk = rt0_SWdir + rt0_SWdif;
SWnet_trunk_dir = rt0_SWdir;
SWnet_trunk_dif = rt0_SWdif;
LWnet_trunk = rt0_LW + rt1_LW*T_trunk_old;    % = LW_atm + LW_loss + (LW_gain-TTnet) + LW_ground
    LW_trunk_atm = rt0_LW_atm;
    LW_trunk_loss = rt0_LW_trunk + rt1_LW*T_trunk_old;
    LW_trunk_gain = rt0_LW_leaf + rt2_LW*T_can;
    LW_trunk_ground = rt0_LW_ground;
TTnet_trunk = rt2_LW*T_can;        % uses new temperature of leaf layer
Rnet_trunk = SWnet_trunk + LWnet_trunk + TTnet_trunk;     % A
HMnet_trunk = HMt0 + HMt1*T_trunk_old;                    % |
SHnet_trunk = ht0 + ht1*T_trunk_old;                      % |
LHnet_trunk = let0 + let1*T_trunk_old;                    % |
dLWnetdT_trunk = rt1_LW;                                  % |
% dTTnetdT_trunk = TT1; % 0 according to calculations, see _|
dHMnetdT_trunk = HMt1;
dSHnetdT_trunk = ht1;
dLHnetdT_trunk = let1;
EB_trunk(1,i) = SWnet_trunk_dir;
EB_trunk(2,i) = SWnet_trunk_dif;
EB_trunk(3,i) = LW_trunk_atm;
EB_trunk(4,i) = LW_trunk_loss;
EB_trunk(5,i) = LW_trunk_gain;
EB_trunk(6,i) = LW_trunk_ground;
EB_trunk(7,i) = TTnet_trunk;
EB_trunk(8,i) = HMnet_trunk;
EB_trunk(9,i) = SHnet_trunk;
EB_trunk(10,i) = LHnet_trunk;
EB_trunk(11,i) = 0;   % dSWnetdT_trunk
EB_trunk(12,i) = dLWnetdT_trunk;
EB_trunk(13,i) = 0;   % dTTnetdT_trunk, s.a.
EB_trunk(14,i) = dHMnetdT_trunk;
EB_trunk(15,i) = dSHnetdT_trunk;
EB_trunk(16,i) = dLHnetdT_trunk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% partition LH flux on Int and Trans, correct EB for overestimated interc. evaporation
    if T_air > T_freez
        Evap_can = dt*3600*LE_can/lh_vap;   % [mm]
    else
        Evap_can = dt*3600*LE_can/lh_subl;  % [mm]
    end

    if Evap_can <= 0
        Evap_int = Evap_can;                % [mm]
        Transp = 0;                         % [mm]
        LE_can_corr = LE_can;
    else 
        Transp = Evap_can*ce_transpiration*(1-frac_wet)/ce_canopy;
        Evap_int = Evap_can - Transp;
        if Evap_int > IntStor
            Evap_int = IntStor;
            Evap_can = Evap_int + Transp;
            if T_air > T_freez
                LE_can_corr = Evap_can*lh_vap/(dt*3600);
            else
                LE_can_corr = Evap_can*lh_subl/(dt*3600);
            end	
%{
re-compute T_can from
(R0 + R1 * T_can + TT0 + TT1 * T_can) = (H0 + H1 * T_can + HM0 + HM1 * T_can ) + LE_can_corr
%}
            T_can = (LE_can_corr + h0 - r0 + HM0 - TT0)/(r1 - h1- HM1 + TT1);
            T_trunk = T_can*TT1/r2 + TT0/r2;
% re-compute RN_can, H_can, and LE_can with new temperature
            RN_can = r0 + r1  * T_can + r2 * T_trunk;
            H_can  = h0 + h1 * T_can;
            LE_can = LE_can_corr;
        else
            LE_can_corr = LE_can;
        end
    end

    IntStor_new = IntStor - Evap_int;

% wet surface fraction
    if IntStor_new > 0
    % limit the wet fraction to minimum 0.01 otherwise it will never completely dry
        frac_wet_new = max(0.01,min(1,(IntStor_new/IntCapacity)^(2/3)));
    else
        frac_wet_new = 0;
    end
%{
Changes of temperature induce changes in stability correction. Re-computation
of turbulent exchange coefficients is needed in case of big changes in TC.
%}
    if abs(T_can-T_can_old) > (dT_can_maxperh*dt)
        [ch_canopy,ce_canopy,ce_transpiration,ce_interception,ce_condensation,ra,rs]...
    = TurbulentExchange(SW_in_ac,T_air,wind,RelHum,T_can,LAI,z_can,z_ref,z_0mg,frac_wet_new,SnowDepth,T_ground);
    end

    T_can = (T_can+T_can_old)/2;
    frac_wet = (frac_wet+frac_wet_new)/2;
    
end

% now REDUCE WaterContent in the Soil Elements -> excluded here

% final adjustment of interception storage due to evaporation
IntStor = IntStor - Evap_int;

% radiation above and below canopy
[SW_refl_ac,LW_refl_ac,SW_in_bc,LW_in_bc,SW_refl_bc,LW_refl_bc,SW_net_trunk,LW_net_trunk]...
    = CanopyRadiationOutput(SolAngle,SW_in_ac,frac_dir,LW_in_ac,T_ground,T_can,T_trunk,...
    frac_through,emg,z_can,Diameter,sig_leaf,sig_leaf_dir,sig_trunk,sig_trunk_dir,em_leaf,em_trunk,alb_leaf,alb_trunk,alb_ground);



% atmospheric emissivity as seen by surface below canopy with regard to air temperature
em_air_bc = LW_in_bc/(boltz*T_air^4);   % same variable name as "normal" em_air in SNOWPACK

% adjust friction velocity below canopy using the same reference height as in Meteo.c
z_ref = max(0.5,z_wind - SnowDepth);
ustar = 0.74*log(z_ref/z_0mg)/0.4/(ra+rs);

% function returning the total surface albedo of a canopy covered snow or soil surface
    % total surface albedo (diffuse fraction)
alb_diff = (1-frac_dir)*((sig_leaf*alb_leaf + alb_ground*(1-sig_leaf)^2/(1 - sig_leaf*alb_leaf*alb_ground))...
    *(1-frac_through) + alb_ground*frac_through);
	% total surface albedo (direct fraction)
alb_dir = frac_dir*((sig_leaf_dir*alb_leaf + alb_ground*(1-sig_leaf_dir)^2/(1 - sig_leaf_dir*alb_leaf*alb_ground))...
    *CC_dir + alb_ground*(1-CC_dir));
alb_tot =  alb_diff + alb_dir;
                              

end