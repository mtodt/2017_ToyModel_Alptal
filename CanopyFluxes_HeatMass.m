function [dlrad,ulrad,h2ocan,t_veg,taf,t_ref2m,EB_veg]...
    = CanopyFluxes_HeatMass(lat,decl,dtime,elai,esai,htop,...
    displa,z0mv,z0mg,forc_hgt_u_pft,forc_hgt_t_pft,emv,emg,forc_lwrad,...
    thm,thv,forc_th,forc_q,forc_wind,forc_pbot,forc_rho,o2,co2,...
    t_grnd,qg,t_veg,sabv,sabvd,sabvi,fwet,h2ocan,snow_depth,soilbeta,btran,frac_veg_nosno,...
    nrad,tlai_z,vcmaxcintsun,vcmaxcintsha,parsun_z,parsha_z,laisun_z,laisha_z,laisun,laisha)
%{
DESCRIPTION:
Calculates the leaf temperature and the leaf fluxes, transpiration,
photosynthesis and  updates the dew accumulation due to evaporation.

Method:
Use the Newton-Raphson iteration to solve for the foliage temperature that
balances the surface energy budget:
    f(t_veg) = Net radiation - Sensible - Latent = 0
    f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)

Note:
(1) In solving for t_veg, t_grnd is given from the previous timestep.
(2) The partial derivatives of aerodynamical resistances, which cannot
    be determined analytically, are ignored for d(H)/dT and d(LE)/dT
(3) The weighted stomatal resistance of sunlit and shaded foliage is used
(4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
                                                         => Ec + Eg = Ea
(5) Energy loss is due to: numerical truncation of energy budget equation
    (*); and "ecidif" (see the code) which is dropped into the sensible
    heat
(6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
    and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
    of water flux from the leaf between the iteration step (n+1) and (n)
    less than 0.1 W/m2; or the iterative steps over 40.
%}

%--------------------------------------------------------------------------
%-------------------------  Physical Parameters  --------------------------

alpha_aero = 1;     % constant for aerodynamic parameter weighting
tlsai_crit = 2;     % critical value of elai+esai for which aerodynamic parameters are maximum
sb = 5.67*10^(-8);  % Stefan-Boltzmann constant [W m^{-2} K^{-4}]
grav = 9.80616;     % acceleration of gravity [m s^{-2}]
csoilc = 0.004;     % drag coefficient for soil under canopy
ria = 0.5;          % free parameter for stable formulation (currently 0.5)
z_dl = 0.05;        % placeholder for (dry) litter layer thickness [m]
lai_dl = 0.5;       % placeholder for (dry) plant litter area index [m^2 m^{-2}]
cpair = 1004.64;    % specific heat of dry air [J kg^{-1} K^{-1}]
hvap = 2.501*10^6;  % latent heat of evaporation [J kg^{-1}]
vkc = 0.4;          % von Karman constant
zii = 1000;         % convective boundary layer height [m]
beta = 1;           % coefficient of convective velocity
delmax = 1;         % maxchange in leaf temperature [K]
dtmin = 0.01;       % max limit for temperature convergence [K]
dlemin = 0.1;       % max limit for energy flux convergence [W m^{-2}]


%--------------------------------------------------------------------------
%----------------------------  Initialisation  ----------------------------

del    = 0;     % change in leaf temperature from previous iteration
efeb   = 0;     % latent head flux from leaf for previous iteration
wtlq0  = 0;
wtalq  = 0;
wtgq   = 0;
wtaq0  = 0;
obuold = 0;


%--------------------------------------------------------------------------
%----------------------------  Calculations  ----------------------------

% calculate daylength [s] control for Vcmax
temp = -(sin(lat*pi/180)*sin(decl))/(cos(lat*pi/180)*cos(decl));
temp = min(1,max(-1,temp));
dayl = 2*13750.9871*acos(temp);
    max_decl = 0.409571;    % +/- 23.4667 degrees = +/- 0.409571 radians
    if lat < 0
        max_decl = -max_decl;
    end
    temp = -(sin(lat*pi/180)*sin(max_decl))/(cos(lat*pi/180)*cos(max_decl));
    temp = min(1,max(-1,temp));
    max_dayl = 2*13750.9871*acos(temp);
dayl_factor = min(1,max(0.01,(dayl^2)/(max_dayl^2)));

rb1 = 0;

% modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
lt = min(elai+esai,tlsai_crit);
egvf = (1 - alpha_aero*exp(-lt))/(1 - alpha_aero*exp(-tlsai_crit));
displa = egvf*displa;   % displacement height [m]
z0mv   = exp(egvf*log(z0mv) + (1-egvf)*log(z0mg));
z0hv   = z0mv;
z0qv   = z0mv;

% net absorbed longwave radiation by canopy and ground = air+bir*t_veg^4+cir*t_grnd^4
air =  emv*(1+(1-emv)*(1-emg))*forc_lwrad;
bir = -(2 - emv*(1-emg))*emv*sb;
cir =  emv*emg*sb;

% saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
[el,deldT,qsatl,qsatldT] = QSat(t_veg,forc_pbot);

% initialise flux profile
nmozsgn = 0;
taf = (t_grnd + thm)/2;
qaf = (forc_q + qg)/2;
ur = max(1,forc_wind);
dth = thm - taf;
dqh = forc_q - qaf;
delq = qg - qaf;
dthv = dth*(1 + 0.61*forc_q) + 0.61*forc_th*dqh;
zldis = forc_hgt_u_pft - displa;    % always >0 since measurements above canopy

% initialisation of Monin-Obukhov length
    % initial values of u* and convective velocity
%    ustar = 0.06;
    wc = 0.5;
    if dthv >= 0
        um = max(ur,0.1);
    else
        um = sqrt(ur^2 + wc^2);
    end
    rib = grav*zldis*dthv/(thv*um*um);
    if rib >= 0     % neutral or stable
        zeta = rib*log(zldis/z0mv)/(1 - 5*min(rib,0.19));
        zeta = min(2,max(zeta,0.01));
    else            % unstable
        zeta = rib*log(zldis/z0mv);
        zeta = max(-100,min(zeta,-0.01));
    end
    obu = zldis/zeta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% include heat mass (values from Alptal site)
BM_heat_cap = 2800;     % [J kg^{-1} K^{-1}]
BM_density = 900;       % [kg m^{-3}]
BasalArea = 0.004;      % Canopy Basal Area [m2/m2], value for closed canopies like Alptal
HM_leaves = 0.001*elai*BM_density*BM_heat_cap;
HM_trunk= 0.5*BasalArea*htop*BM_density*BM_heat_cap;
scalingfactor = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_veg_pre = t_veg;

i = 0;
% Energy Balance output: size stems from SNOWPACK + err  x  iterations
EB_veg = nan(17,41);

% stability iteration
while i <= 40
    % determine friction velocity, and potential temperature and humidity profiles of the surface boundary layer
    [ustar,temp1,temp2,temp12m,temp22m]...
        = FrictionVelocity(displa,z0mv,z0hv,obu,um,forc_hgt_u_pft,forc_hgt_t_pft);
    
    tlbef = t_veg;
    del2 = del;
    
    % determine aerodynamic resistances
    ram1 = 1/(ustar^2/um);
    rah(1) = 1/(temp1*ustar);
    raw(1) = 1/(temp2*ustar);
    
    % determine bulk boundary layer resistances of leaves
    uaf = um*sqrt(1/(ram1*um));
    dleaf = 0.04;   % characteristic leaf dimension [m] (same value for all PFTs)
    cf = 0.01/(sqrt(uaf)*sqrt(dleaf));
    rb = 1/(cf*uaf);
    rb1 = rb;
    
    % parameterisation for variation of csoilc with canopy density
    w = exp(-(elai+esai));
    % transfer coefficient over bare soil is changed to a local variable for readability
    csoilb = vkc/(0.13*(z0mg*uaf/(1.5*10^(-5)))^0.45);
    % compute the stability parameter for ricsoilc
    ri = (grav*htop*(taf-t_grnd))/(taf*uaf^2);
    % modify csoilc value (0.004) if the under-canopy is in stable condition
    if taf-t_grnd > 0
        ricsoilc = csoilc/(1 + ria*min(ri,10));
        csoilcn = csoilb*w + ricsoilc*(1-w);
    else
        csoilcn = csoilb*w + csoilc*(1-w);
    end
    
    rah(2) = 1/(csoilcn*uaf);
    raw(2) = rah(2);
    
    svpts = el;                 % [Pa]
    eah = forc_pbot*qaf/0.622;  % [Pa]
    
    rhaf = eah/svpts;
    rhal = rhaf;
    vpdal = svpts - eah;
    
    rssun = Photosynthesis(forc_pbot,thm,btran,svpts,eah,o2,co2,rb,dayl_factor,t_veg,nrad,tlai_z,vcmaxcintsun,parsun_z,laisun_z);
    rssha = Photosynthesis(forc_pbot,thm,btran,svpts,eah,o2,co2,rb,dayl_factor,t_veg,nrad,tlai_z,vcmaxcintsha,parsha_z,laisha_z);
    
    % sensible heat conductance for air, leaf and ground
    wta = 1/rah(1);
    wtl = (elai+esai)/rb;
    wtg = 1/rah(2);
    wtshi = 1/(wta+wtl+wtg);
    
    wta0 = wta*wtshi;
    wtl0 = wtl*wtshi;
    wtg0 = wtg*wtshi;
    
    wtga = wta0 + wtg0;
    wtal = wta0 + wtl0;

    % fraction of potential evaporation from leaf
    if fwet < 1
        rppdry = (1-fwet)*rb*(laisun/(rb+rssun) + laisha/(rb+rssha))/elai;
    else
        rppdry = 0;
    end
    efpot = forc_rho*wtl*(qsatl - qaf);
    if efpot > 0
        if btran > 0
            qflx_tran_veg = efpot*rppdry;
            rpp = rppdry + fwet;
        else    % no transpiration if btran below 1*10^{-10}
            rpp = fwet;
            qflx_tran_veg = 0;
        end
        % check total transpiration from leaves
        rpp = min(rpp,(qflx_tran_veg+h2ocan/dtime)/efpot);
    else    % no transpiration if potential evaporation less than zero
        rpp = 1;
        qflx_tran_veg = 0;
    end
    
    % update conductances for changes in rpp
    % latent heat conductances for ground & leaf, equal air conductances for SH and LH
    wtaq = frac_veg_nosno * 1/raw(1);
    wtlq = frac_veg_nosno * rpp*(elai+esai)/rb;
    
    % litter layer resistance 
    snow_depth_c = z_dl;            % critical depth for 100% litter burial by snow (=litter thickness)
    fsno_dl = snow_depth/snow_depth_c;  % effective snow cover for (dry) plant litter
    elai_dl = lai_dl*(1 - min(fsno_dl,1));  % exposed (dry) litter area index
    rdl = (1 - exp(-elai_dl))/(0.004*uaf);  % dry litter layer resistance
    % add litter layer resistance
    if delq < 0
        wtgq = frac_veg_nosno * 1/(raw(2)+rdl);
    else
        wtgq = frac_veg_nosno * soilbeta/(raw(2)+rdl);   % soilbeta depending on water content of soil layer
    end
    wtsqi = 1/(wtaq+wtlq+wtgq);
    
    wtaq0 = wtaq*wtsqi;
    wtlq0 = wtlq*wtsqi;
    wtgq0 = wtgq*wtsqi;
    
    wtgaq = wtaq0 + wtgq0;
    wtalq = wtaq0 + wtlq0;
    
    dc1 = forc_rho*cpair*wtl;
    dc2 = hvap*forc_rho*wtlq;
    
    efsh = dc1*(wtga*t_veg - wtg0*t_grnd - wta0*thm);
    efe = dc2*(wtgaq*qsatl - wtgq0*qg - wtaq0*forc_q);
    
    % evaporation flux from foliage
    erre = 0;
    if efe*efeb < 0
        efeold = efe;
        efe = 0.1*efeold;
        erre = efe - efeold;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % include heat mass -> conductive heat flux
        % efhm_leaves = scalingfactor*HM_leaves/dt * (T_can - T_can_pre);
        % efhm_trunk = scalingfactor*HM_trunk/dt * (T_trunk - T_trunk_pre);
    efhm = scalingfactor*(HM_leaves+HM_trunk)/dtime * (t_veg - t_veg_pre);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    Fractioning ground LW emittance not done, because at sites snow fraction
    either 1 or 0 and correct ground/snow temperature passed in anyway.
    %}
    dt_veg = (sabv + air + bir*t_veg^4 + cir*t_grnd^4 - efsh - efe - efhm)/...
        (-4*bir*t_veg^3 + dc1*wtga + dc2*wtgaq*qsatldT + scalingfactor*(HM_leaves+HM_trunk)/dtime);
    t_veg = tlbef + dt_veg;
    dels = dt_veg;
    del = abs(dels);
    err = 0;
    if del > delmax
        dt_veg = delmax*dels/del;
        t_veg = tlbef + dt_veg;
        err = sabv + air + bir*tlbef^3*(tlbef + 4*dt_veg) + cir*t_grnd^4 ...
            - (efsh + dc1*wtga*dt_veg) - (efe + dc2*wtgaq*qsatldT*dt_veg - (efhm + scalingfactor*(HM_leaves+HM_trunk)/dtime * dt_veg));
    end
    
    %{
    fluxes from leaves to canopy space
    efe was limited as its sign changes frequently. This limit may result
    in an imbalance in  hvap*qflx_evap_veg  and  efe + dc2*wtgaq*qsatldT*dt_veg
    %}
    efpot = forc_rho*wtl*(wtgaq*(qsatl + qsatldT*dt_veg) - wtgq0*qg - wtaq0*forc_q);
    qflx_evap_veg = rpp*efpot;
    
    % calculation of evaporation potentials (efpot) and interception losses, flux in [kg m^{-2} s^{-1}]
    %{
    ecidif holds the excess energy if all intercepted water is evaporated
    during the timestep. This energy is later added to the sensible heat
    flux.
    %}
    ecidif = 0;
    if efpot > 0 && btran > 0
        qflx_tran_veg = efpot*rppdry;
    else
        qflx_tran_veg = 0;
    end
    ecidif = max(0,qflx_evap_veg - qflx_tran_veg - h2ocan/dtime);
    qflx_evap_veg = min(qflx_evap_veg,qflx_tran_veg + h2ocan/dtime);
    
    % The energy loss due to above two limits is added to the sensible heat flux.
    eflx_sh_veg = efsh + dc1*wtga*dt_veg + err + erre + hvap*ecidif;
    
    % re-calculate saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
    [el,deldT,qsatl,qsatldT] = QSat(t_veg,forc_pbot);
    
    %{
update vegetation/ground surface temperature, canopy air temperature,
canopy vapor pressure, aerodynamic temperature, and Monin-Obukhov
stability parameter for next iteration.
    %}
    taf = wtg0*t_grnd + wta0*thm + wtl0*t_veg;
    qaf = wtlq0*qsatl + wtgq0*qg + wtaq0*forc_q;
    
    % update Monin-Obukhov length and wind speed including the stability effect
    dth = thm - taf;
    dqh = forc_q - qaf;
    delq = wtalq*qg - wtlq0*qsatl - wtaq0*forc_q;
    
    tstar = temp1*dth;
    qstar = temp2*dqh;
    
    thvstar = tstar*(1 + 0.61*forc_q) + 0.61*forc_th*qstar;
    zeta = zldis*vkc*grav*thvstar/(thv*ustar^2);
    
    if zeta >= 0    % stable
        zeta = min(2,max(zeta,0.01));
        um = max(ur,0.1);
    else            % unstable
        zeta = max(-100,min(zeta,-0.01));
        wc = beta*(-grav*ustar*thvstar*zii/thv)^(1/3);
        um = sqrt(ur^2 + wc^2);
    end
    obu = zldis/zeta;
    if obuold*obu < 0
        nmozsgn = nmozsgn+1;
    end
    if nmozsgn >= 4
        obu = zldis/(-0.01);
        obuold = obu;
    end
    
    % Test For Convergence
    i=i+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adding individual terms of energy balance to output
SWnet = sabv;
SWnet_dir = sabvd;
SWnet_dif = sabvi;
LWnet = air + bir*tlbef^4 + cir*t_grnd^4;
    LW_atm = air;
    LW_loss = bir*tlbef^4;
    LW_ground = cir*t_grnd^4;
SHnet = efsh;
LHnet = efe;
HMnet = efhm;
dLWnetdTveg = 4*bir*tlbef^3;
dHMnetdTveg = scalingfactor*(HM_leaves+HM_trunk)/dtime;
dSHnetdTveg = dc1*wtga;
dLHnetdTveg = dc2*wtgaq*qsatldT;
EB_veg(1,i) = SWnet_dir;
EB_veg(2,i) = SWnet_dif;
EB_veg(3,i) = LW_atm;
EB_veg(4,i) = LW_loss;
EB_veg(5,i) = 0;  % no gain from second layer
EB_veg(6,i) = LW_ground;
EB_veg(7,i) = 0;  % no interaction with another layer
EB_veg(8,i) = HMnet;
EB_veg(9,i) = SHnet;
EB_veg(10,i) = LHnet;
EB_veg(11,i) = 0; % dSWnetdTveg
EB_veg(12,i) = dLWnetdTveg;
EB_veg(13,i) = 0; % no interaction with another layer
EB_veg(14,i) = dHMnetdTveg;
EB_veg(15,i) = dSHnetdTveg;
EB_veg(16,i) = dLHnetdTveg;
EB_veg(17,i) = err;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i > 2
        dele = abs(efe - efeb);
        efeb = efe;
        det = max(del,del2);
        if det < dtmin && dele < dlemin   % convergence criterium
            i=41;   % if understood correctly, PFT in CLM4.5 not re-iterated anymore when above criterium is satisifed
        end
    end
end

% energy balance check
err = sabv + air + bir*tlbef^3 * (tlbef + 4*dt_veg) + cir*t_grnd^4 ...
    - eflx_sh_veg - hvap*qflx_evap_veg -  scalingfactor*(HM_leaves+HM_trunk)/dtime * dt_veg;

% fluxes from ground to canopy space
delt = wtal*t_grnd - wtl0*t_veg - wta0*thm;
tau = -forc_rho*forc_wind/ram1;
eflx_sh_grnd = cpair*forc_rho*wtg*delt;
qflx_evap_soi = forc_rho*wtgq*delq;

% 2m height air temperature
t_ref2m = thm + temp1*dth*(1/temp12m - 1/temp1);

% 2m height specific humidity
q_ref2m = forc_q + temp2*dqh*(1/temp22m - 1/temp2);

% 2m height relative humidity
[e_ref2m,de2mdT,qsat_ref2m,dqsat2mdT] = QSat(t_ref2m,forc_pbot);
rh_ref2m = min(100,q_ref2m/qsat_ref2m * 100);

% downward longwave radiation below the canopy
%dlrad = (1-emv)*emg*forc_lwrad + emv*emg*sb*tlbef^3 * (tlbef + 4*dt_veg);
dlrad = (1-emv)*forc_lwrad + emv*sb*tlbef^3 * (tlbef + 4*dt_veg);
% excluded emg since we want radiation downwards not radiation absorbed by the ground

% upward longwave radiation above the canopy
ulrad = (1-emg)*(1-emv)^2 * forc_lwrad ...
    + emv*(1+(1-emg)*(1-emv))*sb*tlbef^3 * (tlbef + 4*dt_veg) ...
    + emg*(1-emv)*sb*t_grnd^4;

% % derivate of soil energy flux with respect to soil temperature
% cgrnds = cgrnds + cpair*forc_rho*wtg*wtal;
% cgrndl = cgrndl + forc_rho*wtgq*wtalq*dqgdT;
% cgrnd = cgrnds + cgrndl*htvp;

% update dew accumulation [kg m^{-2}]
h2ocan = max(0,h2ocan + (qflx_tran_veg - qflx_evap_veg)*dtime);

% % total photosynthesis
% fpsn = psnsun*laisun + psnsha*laisha;
    
end