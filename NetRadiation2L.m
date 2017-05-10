function [r0,r1,r2,rt0,rt1,rt2,r1p,r2p,NR_leaf,NR_trunk,CC_dir,frac_dir,...
    r0_SWdir,r0_SWdif,r0_LW,r0_LW_atm,r0_LW_leaf,r0_LW_trunk,r0_LW_ground,r1_LW,r2_LW,...
    rt0_SWdir,rt0_SWdif,rt0_LW,rt0_LW_atm,rt0_LW_trunk,rt0_LW_leaf,rt0_LW_ground,rt1_LW,rt2_LW] =...
    NetRadiation2L(SW_in_ac,SW_in_ac_dif,SolAngle,LW_in_ac,T_air,T_surf,emg,T_leaf,T_trunk,z_can,LAI,frac_through,frac_wet,alb_surf)
%{
This routine estimates the radiation balance of a vegetation covered grid.
Inputs:
- incoming shortwave (RG) and longwave (RA) radiation
- vegetation temperature (TV) and surface temperature of the ground below (TG)
- vegetation albedo (AV) and ground albedo (AG)
- vegetation shielding coefficient (SIGF) [shortwave and longwave]
- emissivity of vegetation (EV) and ground (EG)
Outputs:
- Net longwave and shortwave radiation of vegetation (RAV,RGV) and ground (RAG,RGG)
- Total grid albedo (AGRID) and grid surface radiation temperature (TGRID)
Objective:
- netradtrunk = rt0 + rt1 * Ttrunk + rt2 * TC
- netradcanop = r0 + r1*TC + r2* Ttrunk
Method:
- using linearisation of TC^4 and Ttrunk^4
%}

%-------------------------- Physical Constants  ---------------------------

boltz = 5.67*10^(-8);       % Stefan-Boltzmann constant [W m^{-2} K^{-4}]
T_melt = 273.15;            % melting temperature [K]

em_surf = emg;              % emissivity of ground
%{
assumption of 1 for both ground and canopy -> suppression of multiple reflections
%}


%---------------------------  Canopy Parameters  --------------------------

alb_can_dry = 0.11;     % albedo of dry canopy (calibr: 0.09, Alptal)
alb_can_wet = 0.11;     % albedo of wet canopy (calibr: 0.09, Alptal)
alb_can_snow = 0.35;    % albedo of snow covered canopy (calibr. Alptal, but 0.3 in Gouttevin et al. (2015))
alb_trunk = 0.09; 		% trunk albedo
em_leaf = 1;            % canopy emissivity
em_trunk = 1; 			% trunk emissivity

% albedo of partly "wet" canopy = weighted average of dry and wet parts
if T_air > T_melt
    alb_can = frac_wet*alb_can_wet + (1-frac_wet)*alb_can_dry;
else
    alb_can = frac_wet*alb_can_snow + (1-frac_wet)*alb_can_dry;
end
alb_leaf = alb_can;



Diameter = 1.0;             % average canopy (tree) diameter [m]
frac_height_trunk = 0.2;  	% fraction of tree height occupied by trunks, for direct insolation
frac_LAI_top = 0.5;         % fraction of LAI attributed to uppermost layer (calibrated for Alptal)
k_LAI = 0.75;               % radiation transmissivity parameter [0.4-0.8] (calibrated on Alptal)

% canopy absorption fators for direct SW or diffuse SW & LW
sig_leaf = 1 - exp(-k_LAI*frac_LAI_top*LAI);
sig_leaf_dir = 1 - exp(-k_LAI*frac_LAI_top*LAI/max(sin(SolAngle),0.0001));
sig_trunk = 1 - exp(-k_LAI*(1-frac_LAI_top)*LAI);
sig_trunk_dir = 1 - exp(-k_LAI*(1-frac_LAI_top)*LAI/max(sin(SolAngle),0.0001));

% attenuation factors for the radiative impact of the trunk layer
att_SW = 1 - sig_trunk;
att_SW_dir = 1 - sig_trunk_dir;
att_LW = 1 - sig_trunk;



% canopy closure
CC = 1-frac_through;    % for diffuse shortwave and longwave radiation
if SolAngle > 0
    CC_dir = min(1,CC*(1+4*z_can/(pi*Diameter*tan(SolAngle))));
else
    CC_dir = 1;         % for direct shortwave radiation
end
CC_dir_leaf = CC_dir;
CC_dir_trunk = 0;

%{
% optional: allows direct solar insolation of the trunks
if SolAngle > 0
    CC_dir_leaf = min(1,CC*(1+4*z_can*(1-frac_height_trunk)/(pi*Diameter*tan(SolAngle))));
else
    CC_dir_leaf = 1;
end
CC_dir_trunk = CC_dir - CC_dir_leaf;
%}


%--------------------------  Radiation Forcing  ---------------------------

SW_in_dif = SW_in_ac_dif;
SW_in_dir = SW_in_ac-SW_in_dif;
if SW_in_dir > 0
    frac_dir = SW_in_dir/(SW_in_dif+SW_in_dir);
    frac_dif = 1-frac_dir;
else
    frac_dir = 0;
    frac_dif = 1;
end

% LW_in_ac = em_air*boltz*T_air^4;    % longwave radiation above canopy

%--------------------------------------------------------------------------
%---------------------------  Radiation Budget  ---------------------------
%--------------------------------------------------------------------------

%------------------------------  leaf layer  ------------------------------

% diffuse shortwave - different formulation than for total SW in Gouttevin et al. (2015)
SW_net_leaf = frac_dif*SW_in_ac*(1-alb_leaf)*sig_leaf*(1 +...
 alb_surf*(1-sig_leaf)*att_SW/(1-sig_leaf*alb_surf*alb_leaf) +...
 (1-sig_leaf)*sig_trunk*alb_trunk);

% longwave
star = 1 - sig_leaf*(1-em_leaf)*(1-em_surf);
psi = (1-sig_leaf)*(1-em_surf);

% total -> net = net_SW + net_LW
r0p = SW_net_leaf + sig_leaf*em_leaf*((1 + att_LW*psi/star)*LW_in_ac + ...
    (att_LW/star)*em_surf*boltz*T_surf^4);
r1p = -1*sig_leaf*em_leaf*boltz*(2 - em_leaf*sig_leaf*att_LW*(1-em_surf)/star);
r2p = sig_leaf*em_leaf*boltz*em_trunk*sig_trunk*(1+(1-em_surf));

T_leaf_old = T_leaf;
T_trunk_old = T_trunk;

r0 = r0p - 3*r1p*T_leaf_old^4 - 3*r2p*T_trunk_old^4;
r1 = 4*r1p*T_leaf_old^3;
r2 = 4*r2p*T_trunk_old^3;

% additional output for EB analysis
r0_SWdif = SW_net_leaf;
r0_LW = r0 - SW_net_leaf;
r0_LW_atm = sig_leaf*em_leaf*(1 + att_LW*psi/star)*LW_in_ac;
r0_LW_leaf = -3*r1p*T_leaf_old^4;
r0_LW_trunk = -3*r2p*T_trunk_old^4;
r0_LW_ground = sig_leaf*em_leaf*(att_LW/star)*em_surf*boltz*T_surf^4;
r1_LW = r1;
r2_LW = r2;

%{
1) r0p, r1p and r2p correpsond to net(t) = r0p + r1p * T_leaf(t)^4 + r2p * T_trunk^4
2) Linearisation of net around T_leaf and T_trunk by using
T_leaf(t)^4 = T_leaf(t-1)^4 + 4*T_leaf(t-1)^3*(T_leaf(t)-T_leaf(t-1)),
which gives us r0, r1, and r2 correpsonding to net(t) = r0 + r1*T_leaf(t)+ r2*Ttrunk
%}

SW_net_leaf = SW_net_leaf*CC;
r0 = r0*CC;
r1 = r1*CC;
r2 = r2*CC;

% additional output for EB analysis
r0_SWdif = r0_SWdif*CC;
r0_LW = r0_LW*CC;
r0_LW_atm = r0_LW_atm*CC;
r0_LW_leaf = r0_LW_leaf*CC;
r0_LW_trunk = r0_LW_trunk*CC;
r0_LW_ground = r0_LW_ground*CC;
r1_LW = r1_LW*CC;
r2_LW = r2_LW*CC;

% addition of direct SW
SW_net_leaf_dir = CC_dir_leaf*frac_dir*SW_in_ac*(1-alb_leaf)*sig_leaf_dir*...
    (1 + alb_surf*(1-sig_leaf_dir)*att_SW_dir/(1 - sig_leaf_dir*alb_surf*alb_leaf)...
    + (1-sig_leaf_dir)*sig_trunk_dir*alb_trunk)...
    + CC_dir_trunk*frac_dir*SW_in_ac*(1-alb_leaf)*sig_leaf_dir*alb_surf*att_SW_dir/...
    (1 - sig_leaf_dir*alb_surf*alb_leaf);
SW_net_leaf = SW_net_leaf + SW_net_leaf_dir;
r0 = r0 + SW_net_leaf_dir;

r0_SWdir = SW_net_leaf_dir;    % additional output for EB analysis

NR_leaf = SW_net_leaf;

%------------------------------  trunk layer  -----------------------------

% diffuse shortwave
SW_net_trunk = frac_dif*SW_in_ac*(1-sig_leaf)*(1-att_SW)*(1-alb_trunk);

% longwave
rt0p = SW_net_trunk + em_trunk*(1-att_LW)*(em_surf*boltz*T_surf^4 + LW_in_ac*(1-sig_leaf));
rt1p = -2*boltz*em_trunk*(1-att_LW);
rt2p = em_trunk*(1-att_LW)*em_leaf*sig_leaf*boltz;

rt0 = rt0p - 3*rt1p*T_trunk_old^4 - 3*rt2p*T_leaf_old^4;
rt1 = 4*rt1p*T_trunk_old^3;
rt2 = 4*rt2p*T_leaf_old^3;

% additional output for EB analysis
rt0_SWdif = SW_net_trunk;
rt0_LW = rt0 - SW_net_trunk;
rt0_LW_atm = em_trunk*(1-att_LW)*(1-sig_leaf)*LW_in_ac;
rt0_LW_trunk = -3*rt1p*T_trunk_old^4;
rt0_LW_leaf = -3*rt2p*T_leaf_old^4;
rt0_LW_ground = em_trunk*(1-att_LW)*em_surf*boltz*T_surf^4;
rt1_LW = rt1;
rt2_LW = rt2;

SW_net_trunk = SW_net_trunk*CC;
rt0 = rt0*CC;
rt1 = rt1*CC;
rt2 = rt2*CC;

% additional output for EB analysis
rt0_SWdif = rt0_SWdif*CC;
rt0_LW = rt0_LW*CC;
rt0_LW_atm = rt0_LW_atm*CC;
rt0_LW_trunk = rt0_LW_trunk*CC;
rt0_LW_leaf = rt0_LW_leaf*CC;
rt0_LW_ground = rt0_LW_ground*CC;
rt1_LW = rt1_LW*CC;
rt2_LW = rt2_LW*CC;

% total & addition of direct shortwave
SW_net_trunk_dir = CC_dir_leaf*frac_dir*SW_in_ac*(1-sig_leaf_dir)*(1-att_SW_dir)*(1-alb_trunk)...
    + CC_dir_trunk*frac_dir*SW_in_ac*(1-att_SW_dir)*(1-alb_trunk);

SW_net_trunk = SW_net_trunk + SW_net_trunk_dir;
rt0 = rt0 + SW_net_trunk_dir;

rt0_SWdir = SW_net_trunk_dir;    % additional output for EB analysis

NR_trunk = SW_net_trunk;


end