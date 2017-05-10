function [SW_refl_ac,LW_refl_ac,SW_in_bc,LW_in_bc,SW_refl_bc,LW_refl_bc,SW_net_trunk,LW_net_trunk] = ...
    CanopyRadiationOutput(SolAngle,SW_in_ac,frac_dir,LW_in_ac,T_sfc,T_leaf,T_trunk,...
    DirectThroughfall,emg,z_can,Diameter,sig_leaf,sig_leaf_dir,sig_trunk,sig_trunk_dir,em_leaf,em_trunk,alb_leaf,alb_trunk,alb_ground)
%{
Computes upward and downward radiation below and above canopy.
%}

%-------------------------  Physical Parameters  --------------------------

em_ground = emg;          % ground emissivity
boltz = 5.67*10^(-8);     % Stefan-Boltzmann constant [W m^{-2} K^{-4}]

%--------------------------------------------------------------------------
%-----------------------------  Calculations  -----------------------------

% attenuation factors for two-layer canopy

att_SW = 1-sig_trunk;
att_LW = 1-sig_trunk;
att_SWdir = 1-sig_trunk_dir;

%------------------------------  Shortwave  -------------------------------

% diffuse shortwave radiation fluxes above and below canopy
SW_refl_ac_loc = SW_in_ac*(sig_leaf*alb_leaf + ((1-sig_leaf)^2)*sig_trunk*alb_trunk +...
    alb_ground*((1-sig_leaf)^2)/(1 - sig_leaf*alb_leaf*alb_ground)*att_SW);
SW_in_bc_loc = SW_in_ac*(1-sig_leaf)/(1 - sig_leaf*alb_leaf*alb_ground)*att_SW ;
SW_refl_bc_loc = SW_in_bc_loc*alb_ground;

% direct shortwave radiation fluxes above and below canopy
SW_refl_ac_loc2 = SW_in_ac*(sig_leaf_dir*alb_leaf + (1-sig_leaf_dir)*(1-sig_leaf_dir)*sig_trunk_dir*alb_trunk +...
    alb_ground*(1-sig_leaf_dir)*(1-sig_leaf_dir)/(1 - sig_leaf_dir*alb_leaf*alb_ground)*att_SWdir);
SW_in_bc_loc2 = SW_in_ac*(1-sig_leaf_dir)/(1 - sig_leaf_dir*alb_leaf*alb_ground)*att_SWdir ;
SW_refl_bc_loc2 = SW_in_bc_loc2*alb_ground;
        
% additional direct shortwave radiation term due to direct insolation on trunks
SW_refl_ac_loc3 = SW_in_ac*(sig_trunk_dir*alb_trunk +...
    alb_ground*(1-sig_leaf_dir)/(1 - sig_leaf_dir*alb_leaf*alb_ground)*att_SWdir);
SW_in_bc_loc3 = SW_in_ac*att_SWdir/(1 - sig_leaf_dir*alb_leaf*alb_ground);
SW_refl_bc_loc3 = SW_in_bc_loc3*alb_ground;


%-------------------------------  Longwave  -------------------------------

% longwave radiation fluxes above and below canopy
RAG = em_ground*(-boltz*T_sfc^4 +...
    ((1-sig_leaf)*LW_in_ac*att_LW + sig_leaf*em_leaf*boltz*(T_leaf^4)*att_LW +...
    em_ground*sig_leaf*(1-em_leaf)*att_LW*boltz*T_sfc^4)...
    /(1 - sig_leaf*(1-em_leaf)*(1-em_ground))) +...
    em_ground*em_trunk*(1-att_LW)*boltz*(T_trunk^4)*(1 + sig_leaf*(1-em_leaf));
RAV = sig_leaf*em_leaf*(LW_in_ac - 2*boltz*T_leaf^4 +...
    att_LW*(boltz*(em_ground*T_sfc^4 + em_leaf*sig_leaf*(T_leaf^4)*(1-em_ground)) +...
    (1-sig_leaf)*(1-em_ground)*LW_in_ac)...
    /(1 - sig_leaf*(1-em_leaf)*(1-em_ground)) +...
    em_trunk*(1-att_LW)*boltz*(T_trunk^4)*(1 + (1-em_ground)));
RAT = -2*em_trunk*(1-att_LW)*boltz*T_trunk^4 + em_trunk*(1-att_LW)*...
    (em_leaf*sig_leaf*boltz*T_leaf^4 + em_ground*boltz*T_sfc^4 + LW_in_ac*(1-sig_leaf));

LW_in_bc = RAG/em_ground + boltz*T_sfc^4;
LW_refl_bc = (1-em_ground)*LW_in_bc + em_ground*boltz*T_sfc^4;
LW_refl_ac = LW_in_ac - RAG - RAV - RAT;


%---------------------  scaling with Canopy Closure  ----------------------

% canopy closure
CC_dif = 1-DirectThroughfall;       % for diffuse shortwave and longwave radiation
if SolAngle > 0
    CC_dir = min(1,CC_dif*(1+4*z_can/(pi*Diameter*tan(SolAngle))));
else
    CC_dir = 1;                     % for direct shortwave radiation
end

CC_dir_leaf = CC_dir;
CC_dir_trunk = 0;
%{
% optional: allows direct solar insolation of the trunks
if SolAngle > 0
    CC_dir_leaf = min(1,CC_dif*(1+4*z_can*(1-frac_height_trunk)/(pi*Diameter*tan(SolAngle))));
else
    CC_dir_leaf = 1;
end
CC_dir_trunk = CC_dir - CC_dir_leaf;
%}
		 
% shortwave fluxes (diffuse)
SW_refl_ac = (SW_refl_ac_loc*CC_dif + SW_in_ac*alb_ground*(1-CC_dif))*(1-frac_dir);
SW_in_bc = (SW_in_bc_loc*CC_dif + SW_in_ac*(1-CC_dif))*(1-frac_dir);
SW_refl_bc = (SW_refl_bc_loc*CC_dif + SW_in_ac*alb_ground*(1-CC_dif))*(1-frac_dir);
		
% shortwave fluxes (direct)
SW_refl_ac = SW_refl_ac + (SW_refl_ac_loc2*CC_dir_leaf + SW_refl_ac_loc3*CC_dir_trunk +...
    SW_in_ac*alb_ground*(1-CC_dir_trunk-CC_dir_leaf))*frac_dir;
SW_in_bc = SW_in_bc + (SW_in_bc_loc2*CC_dir_leaf + SW_in_bc_loc3*CC_dir_trunk +...
    SW_in_ac*(1-CC_dir_trunk-CC_dir_leaf))*frac_dir;
SW_refl_bc = SW_refl_bc + (SW_refl_bc_loc2*CC_dir_leaf + SW_refl_bc_loc3*CC_dir_trunk +...
    SW_in_ac*alb_ground*(1-CC_dir_trunk-CC_dir_leaf))*frac_dir;
		
% longwave fluxes (treat as diffuse)
LW_refl_ac = LW_refl_ac*CC_dif + boltz*em_ground*T_sfc^4*(1-CC_dif);
LW_in_bc = LW_in_bc*CC_dif + LW_in_ac*(1-CC_dif);
LW_refl_bc = LW_refl_bc*CC_dif + boltz*em_ground*T_sfc^4*(1-CC_dif);
		
% radiations to trunks
SW_net_trunk = (1-frac_dir)*SW_in_ac*(1-sig_leaf)*(1-alb_trunk)*(1-att_SW)*CC_dif +...
    CC_dir_leaf*frac_dir*SW_in_ac*(1-sig_leaf)*(1-alb_trunk)*(1-att_SWdir) +...
    CC_dir_trunk*frac_dir*SW_in_ac*(1-alb_trunk)*(1-att_SWdir);
LW_net_trunk = RAT*CC_dif ;

end