function rs = Photosynthesis(forc_pbot,tgcm,btran,esat_tv,eair,oair,cair,rb,dayl_factor,t_veg,nrad,tlai_z,vcmaxcint,par_z,lai_z)
%{
Leaf photosynthesis and stomatal conductance calculation.
All calculations for psn variables are excluded because only the resistances
are necessary for the longwave radiatione/enhancement analysis.

KO option seems suspicious but is used here because
a) these values are stated in the tech notes, and
b) otherwise the 10-day average of the 2m temperature is required for the
calculations.
%}

tfrz = 273.15;  % freezing T of fresh water [K]
rgas = 6.02214*10^(26) * 1.38065*10^(-23);  % universal gas constant [J K^{-1} kmole^{-1}]

% PFT parameters
slatop = 0.008;     % specific leaf area at top of canopy, projected area basis [m^2/gC]
    % slatop = 0.024;   % for NDBTs
    % slatop = 0.03;    % for BDBTs
leafcn = 40;        % leaf C:N [gC/gN]
    % leafcn = 25;      % for NDBTs
    % leafcn = 25;      % for BDBTs
flnr = 0.0466;      % fraction of leaf N in the Rubisco enzyme [gN Rubisco / gN leaf]
    % flnr = 0.0546;    % for NDBTs
    % flnr = 0.1007;    % for BDBTs

%{
===========================================================================
Photosynthesis and stomatal conductance parameters
===========================================================================
%}
% vcmax25 parameters, from CN
fnr = 7.16;
act25 = 3.6;    % umol/mgRubisco/min
    % convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
act25 = act25 * 1000/60;

% activation energy
vcmaxha = 65330;    % KO - vcmaxha = 72000; otherwise
jmaxha  = 43540;    % KO - jmaxha  = 50000; otherwise
tpuha   = 53100;    % KO - tpuha   = 72000; otherwise
    tpuha = vcmaxha;    % same value in tech notes !!!
lmrha   = 46390;

% high temperature deactivation
% The factor "c" scales the deactivation to a value of 1.0 at 25C.
vcmaxhd = 149250;	% KO - vcmaxhd = 200000; otherwise
jmaxhd  = 152040;	% KO - jmaxhd  = 200000; otherwise
tpuhd   = 150650;	% KO - tpuhd   = 200000; otherwise
    tpuhd = vcmaxhd;    % same value in tech notes !!!
lmrhd   = 150650;
vcmaxse = 485;  % KO - otherwise 10-day average of temperature needed later on for calculation
jmaxse  = 495;  % KO - otherwise 10-day average of temperature needed later on for calculation
tpuse   = 490;  % KO - go along with above two (whatever KO means???)
    tpuse = vcmaxse;    % same value in tech notes !!!
lmrse   = 490;
% KO   vcmaxc = fth25(vcmaxhd,vcmaxse);
% KO   jmaxc  = fth25(jmaxhd,jmaxse);
% KO   tpuc   = fth25(tpuhd,tpuse);
lmrc   = fth25(lmrhd, lmrse);

% miscellaneous parameters
fnps = 0.15;
theta_psii = 0.7;

c3flag = 1; % for C3 PFTs (all (boreal) trees)
if c3flag == 1
    bbbopt = 10000;
    mbbopt = 9;
else
    bbbopt = 40000;
    mbbopt = 4;
end

% soil water stress applied to Ball-Berry parameters
bbb = max(bbbopt*btran,1);
mbb = mbbopt;

%{
Multi-layer parameters scaled by leaf nitrogen profile.
Loop through each canopy layer to calculate nitrogen profile using cumulative
lai at the midpoint of the layer.
%}
% leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
lnc = 1/(slatop*leafcn);

% vcmax25 at canopy top, as in CN but using lnc at top of the canopy
vcmax25top = lnc*flnr*fnr*act25*dayl_factor;    % additional scaling depending on CN mode

% parameters derived from vcmax25top, used jmax25 = 1.97 vcmax25
jmax25top = 1.97*vcmax25top;    % KO - otherwise t10 required -v
    % jmax25top = (2.59 - 0.035*min(max((t10-tfrz),11),35))*vcmax25top;
tpu25top = 0.167*vcmax25top;
kp25top = 20000*vcmax25top;

%{
Nitrogen scaling factor. Bonan et al (2011) used kn = 0.11. Here, derive kn
from vcmax25 as in Lloyd et al (2010).
Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
But not used as defined here if using sun/shade big leaf code. Instead,
will use canopy integrated scaling factors from SurfaceAlbedo.
%}
if dayl_factor == 0
    kn =  0;
else
    kn = exp(0.00963*vcmax25top/dayl_factor - 2.43);
end

% leaf maintenance respiration in proportion to vcmax25top
if c3flag == 1
    lmr25top = vcmax25top*0.015;
else
    lmr25top = vcmax25top*0.025;
end

% Loop through canopy layers (above snow). Respiration needs to be calculated
% every timestep. Others are calculated only if daytime.
laican = 0;
for iv = 1:nrad
    % cumulative lai at middle of layer
    if iv == 1
        laican = 0.5*tlai_z(iv);
    else
        laican = laican + 0.5*(tlai_z(iv-1)+tlai_z(iv));
    end
    % Scale for leaf nitrogen profile. If multi-layer code, use explicit
    % profile. If sun/shade big leaf code, use canopy integrated factor.
    if nrad == 1    % nlevcan in FORTRAN code but should be the same
        nscaler = vcmaxcint;
    elseif nrad > 1
        nscaler = exp(-kn*laican);
    end
    % maintenance respiration
    lmr25 = lmr25top*nscaler;
    if c3flag == 1
        lmr_z(iv) = lmr25*ft(t_veg,lmrha)*fth(t_veg,lmrhd,lmrse,lmrc);
    else
        lmr_z(iv) = lmr25 * 2^((t_veg-(tfrz+25))/10);
        lmr_z(iv) = lmr_z(iv)/(1 + exp(1.3*(t_veg-(tfrz+55))));
    end
    if par_z(iv) <= 0   % night time
        vcmax_z(iv) = 0;
        jmax_z(iv) = 0;
        tpu_z(iv) = 0;
        kp_z(iv) = 0;
    else                % day time
        vcmax25 = vcmax25top*nscaler;
        jmax25 = jmax25top*nscaler;
        tpu25 = tpu25top*nscaler;
        kp25 = kp25top*nscaler;
        % adjust for temperature
        % vcmaxse = 668.39 - 1.07 * min(max((t10-tfrz),11),35);
        % jmaxse  = 659.70 - 0.75 * min(max((t10-tfrz),11),35);
        % tpuse = vcmaxse;
        vcmaxc = fth25(vcmaxhd, vcmaxse);
        jmaxc  = fth25(jmaxhd, jmaxse);
        tpuc   = fth25(tpuhd, tpuse);
        vcmax_z(iv) = vcmax25 * ft(t_veg,vcmaxha) * fth(t_veg,vcmaxhd,vcmaxse,vcmaxc);
        jmax_z(iv) = jmax25 * ft(t_veg,jmaxha) * fth(t_veg,jmaxhd,jmaxse,jmaxc);
        tpu_z(iv) = tpu25 * ft(t_veg,tpuha) * fth(t_veg,tpuhd,tpuse,tpuc);
        if c3flag == 0
            vcmax_z(iv) = vcmax25 * 2^((t_veg-(tfrz+25))/10);
            vcmax_z(iv) = vcmax_z(iv)/(1 + exp( 0.2*((tfrz+15)-t_veg)));
            vcmax_z(iv) = vcmax_z(iv)/(1 + exp( 0.3*(t_veg-(tfrz+40))));
        end
        kp_z(iv) = kp25 * 2^((t_veg-(tfrz+25))/10);
    end
    % adjust for soil water
    vcmax_z(iv) = vcmax_z(iv)*btran;
    lmr_z(iv) = lmr_z(iv)*btran;
end

%{
===========================================================================
Leaf-level photosynthesis and stomatal conductance
===========================================================================
%}
rsmax0 = 2*10^4;    % maximum stomatal resistance [s/m]

% leaf boundary layer conductance, umol/m**2/s
cf = forc_pbot/(rgas*1*10^(-3)*tgcm)*1*10^6;
gb = 1/rb;
gb_mol = gb * cf;

% loop through canopy layers (above snow). Only do calculations if daytime
for iv = 1:nrad
    if par_z(iv) <= 0   % night time
        ac(iv) = 0;
        aj(iv) = 0;
        ap(iv) = 0;
        ag(iv) = 0;
        an(iv) = ag(iv) - lmr_z(iv);
        rs_z(iv) = min(rsmax0,cf/bbb);
        ci_z(iv) = 0;
        rh_leaf = 0;
    else                % day time
    % now the constraint is no longer needed, Jinyun Tang
    ceair = min(eair,esat_tv);
    rh_can = ceair/esat_tv;

    % Electron transport rate for C3 plants. Convert par from W/m2 to 
    % umol photons/m**2/s using the factor 4.6
    qabs  = 0.5*(1-fnps)*par_z(iv)*4.6;
    aquad = theta_psii;
    bquad = -(qabs + jmax_z(iv));
    cquad = qabs * jmax_z(iv);
    [r1,r2] = quadratic(aquad,bquad,cquad);
    je = min(r1,r2);

    % iterative loop for ci beginning with initial guess
    if c3flag == 1
        ci_z(iv) = 0.7*cair;
    else
        ci_z(iv) = 0.4*cair;
    end
    niter = 0;

    % Increment iteration counter. Stop if too many iterations
    niter = niter+1;

    % save old ci
    ciold = ci_z(iv);
    
    % find ci and stomatal conductance
    [ciold,gs_mol(iv),an(iv)] = hybrid(ciold,gb_mol,je,cair,oair,lmr_z(iv),par_z(iv),rh_can,tpu_z(iv),vcmax_z(iv),kp_z(iv),bbb,mbb,forc_pbot,t_veg);
                    % ^- correct as output ??? but required below and there's no other location where an is calculated
    
    % End of ci iteration.  Check for an < 0, in which case gs_mol = bbb
    if an(iv) < 0
        gs_mol(iv) = bbb;
    end

    if gs_mol(iv) < 0
        gs_mol(iv) = nan;
    end
    
    % final estimates for cs and ci (needed for early exit of ci iteration when an < 0)
    cs = cair - 1.4/gb_mol * an(iv) * forc_pbot;
    cs = max(cs,1*10^(-6));
    ci_z(iv) = cair - an(iv)*forc_pbot*(1.4*gs_mol(iv)+1.6*gb_mol)/(gb_mol*gs_mol(iv));

    % convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)
    gs = gs_mol(iv)/cf;
    rs_z(iv) = min(1/gs,rsmax0);

    end
end

%{
===========================================================================
Canopy photosynthesis and stomatal conductance
===========================================================================
Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
unit leaf area), which are used in other parts of the model. Here, laican
sums to either laisun or laisha.
%}
gscan = 0;
laican = 0;
for iv = 1:nrad
    gscan = gscan + lai_z(iv)/(rb+rs_z(iv));
    laican = laican + lai_z(iv);
end
if laican > 0
    rs = laican/gscan - rb;
else
    rs = 0;
end