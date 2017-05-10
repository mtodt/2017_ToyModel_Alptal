function [abs_veg,abs_veg_dir,abs_veg_dif,abs_ground,abs_tot,LAI_sun,LAI_sha,parsun_z,parsha_z,LAI_sun_z,LAI_sha_z]...
    = SurfaceRadiation(elai,esai,G_dir,coszen,forc_sol_dir,forc_sol_dif,waveband,...
    alb_ground_dir,alb_ground_dif,omega,fabd,fabi,ftdd,ftid,ftii,...
    nrad,tlai_z,frac_sun_z,fabd_sun_z,fabd_sha_z,fabi_sun_z,fabi_sha_z)
%{
DESCRIPTION:
Calculate solar fluxes absorbed by vegetation and ground surface.

Only considering soil.
%}

% PFT-specific values - differ between evergreen and decidious boreals trees
SLA_top = 0.008;     % specific leaf area at top of canopy, projected area basis [m^2/gC]
dSLAdLAI = 0.001;    % dSLA/dLAI, projected area basis [m^2/gC]
    
% loop over pfts to calculate fsun, etc.
vai = elai+esai;

% initialisation
abs_ground = 0; % solar radiation absorbed by ground (W/m^2)
abs_veg = 0;    % solar radiation absorbed by vegetation (W/m^2)
abs_veg_dir = 0;    % direct solar radiation absorbed by vegetation (W/m^2)
abs_veg_dif = 0;    % diffuse solar radiation absorbed by vegetation (W/m^2)
abs_tot = 0;    % solar radiation absorbed (total) (W/m^2)
parsun_z = zeros(nrad,1);   
parsha_z = zeros(nrad,1);
LAI_sun_z = zeros(nrad,1);
LAI_sha_z = zeros(nrad,1);
%{
if coszen>0 && elai>0 && G_dir>0
    cosz = max(0.001,coszen);
    ext = G_dir/cosz;
    t1 = min(ext*elai,40);
    t2 = exp(-t1);
    frac_sun = (1-t2)/t1;   % sunlit fraction of canopy
    if elai> 0.01
        LAI_sun = elai*frac_sun;      % sunlit leaf area
        LAI_sha = elai*(1-frac_sun);  % shaded leaf area
% specific leaf area for sunlit canopy, projected area basis (m^2/gC)
        SLA_sun = (t2*dSLAdLAI*ext*elai + t2*dSLAdLAI + ...
            t2*SLA_top*ext - dSLAdLAI - SLA_top*ext)/(ext*(t2-1));
% specific leaf area for shaded canopy, projected area basis (m^2/gC)
        SLA_sha = ((SLA_top + (dSLAdLAI*elai/2))*elai - LAI_sun*SLA_sun)/LAI_sha;
    else
        frac_sun = 1;
        LAI_sun = elai;
        LAI_sha = 0;
        SLA_sun = SLA_top;
        SLA_sha = 0;
    end
else
    frac_sun = 0;
    LAI_sun = 0;
    LAI_sha = elai;
    SLA_sun = 0;
    SLA_sha = 0;
end
%}
%{
Loop over pfts to calculate laisun_z and laisha_z for each layer.
Derive canopy laisun, laisha, and fsun from layer sums.
If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
%}
LAI_sun = 0;
LAI_sha = 0;
for iv = 1:nrad
    LAI_sun_z(iv) = tlai_z(iv) * frac_sun_z(iv);
    LAI_sha_z(iv) = tlai_z(iv) * (1-frac_sun_z(iv));
    LAI_sun = LAI_sun + LAI_sun_z(iv);
    LAI_sha = LAI_sha + LAI_sha_z(iv);
end
if elai > 0
    frac_sun = LAI_sun/elai;
else
    frac_sun = 0;
end

for i=1:waveband
% absorbed by canopy
    cad(i) = forc_sol_dir(i)*fabd(i);
    cai(i) = forc_sol_dif(i)*fabi(i);
    abs_veg = abs_veg + cad(i) + cai(i);
    abs_veg_dir = abs_veg_dir + cad(i);     % additional output for EB analysis
    abs_veg_dif = abs_veg_dif + cai(i);     % additional output for EB analysis
    abs_tot  = abs_tot  + cad(i) + cai(i);
    if i == 1
        parveg = cad(i) + cai(i);
    end
%{
absorbed PAR profile through canopy
If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo are canopy
integrated so that layer values equal big leaf values.
%}              
    if i == 1
        for iv = 1:nrad
            parsun_z(iv) = forc_sol_dir(i)*fabd_sun_z(iv) + forc_sol_dif(i)*fabi_sun_z(iv);
            parsha_z(iv) = forc_sol_dir(i)*fabd_sha_z(iv) + forc_sol_dif(i)*fabi_sha_z(iv);
        end
    end
% transmitted = solar fluxes incident on ground
    trd(i) = forc_sol_dir(i)*ftdd(i);
    tri(i) = forc_sol_dir(i)*ftid(i) + forc_sol_dif(i)*ftii(i);
% solar radiation absorbed by ground surface
    Abs  = trd(i)*(1-alb_ground_dir(i)) + tri(i)*(1-alb_ground_dif(i));
    abs_ground = abs_ground + Abs;
    abs_tot  = abs_tot  + Abs;
    %{
% sunlit/shaded canopy algorithm
    if coszen>0 && elai>0 && G_dir>0
    % 1) calculate flux of direct beam radiation absorbed in the sunlit canopy as direct (VegAbs_sun_dir_dir),
    %    and the flux of direct beam radiation absorbed in the total canopy as indirect
        sun_add(i) = forc_sol_dir(i)*(1-ftdd(i))*(1-omega(i));
        tot_aid(i) = forc_sol_dir(i)*fabd(i) - sun_add(i);
        tot_aid(i) = max(tot_aid(i),0);
    % 3) calculate the fraction of indirect radiation being absorbed in the sunlit and shaded canopy fraction.
    %    Some of this indirect originates in the direct beam and some originates in the indirect beam. [ 2) excluded]
        frac_sun_dif_dir(i) = frac_sun;
        frac_sun_dif_dif(i) = frac_sun;
        frac_sha_dif_dir(i) = 1-frac_sun_dif_dir(i);
        frac_sha_dif_dif(i) = 1-frac_sun_dif_dif(i);
    % 4) calculate the total indirect flux absorbed by the sunlit and shaded canopy based on these
    %    fractions and the abs_veg_dir and abs_veg_dif from surface albedo calculations
        sun_aid(i) = tot_aid(i)*frac_sun_dif_dir(i);
        sun_aii(i) = forc_sol_dif(i)*fabi(i)*frac_sun_dif_dif(i);
        sha_aid(i) = tot_aid(i)*frac_sha_dif_dir(i);
        sha_aii(i) = forc_sol_dif(i)*fabi(i)*frac_sha_dif_dif(i);
    % 5) calculate the total flux absorbed in the sunlit and shaded canopy as the sum of these terms
        sun_atot(i) = sun_add(i) + sun_aid(i) + sun_aii(i);
        sha_atot(i) = sha_aid(i) + sha_aii(i);
    % 6) calculate the total flux absorbed by leaves in the sunlit and shaded canopies
        frac_LAI = elai/vai;
        sun_alf(i) = sun_atot(i)*frac_LAI;
        sha_alf(i) = sha_atot(i)*frac_LAI;
    % 7) calculate the fluxes per unit LAI in the sunlit and shaded canopies
        if LAI_sun > 0
            sun_aperlai(i) = sun_alf(i)/LAI_sun;
        else
            sun_aperlai(i) = 0;
        end
        if LAI_sha > 0
            sha_aperlai(i) = sha_alf(i)/LAI_sha;
        else
            sha_aperlai(i) = 0;
        end
        else        % coszen = 0 or elai = 0
            sun_add(i) = 0;
            tot_aid(i) = 0;
            frac_sun_dif_dir(i) = 0;
            frac_sun_dif_dif(i) = 0;
            frac_sha_dif_dir(i) = 0;
            frac_sha_dif_dif(i) = 0;
            sun_aid(i) = 0;
            sun_aii(i) = 0;
            sha_aid(i) = 0;
            sha_aii(i) = 0;
            sun_atot(i) = 0;
            sha_atot(i) = 0;
            sun_aperlai(i) = 0;
            sha_aperlai(i) = 0;
            
            
    end
end
% 8) calculate total and per-unit-lai fluxes for PAR in sunlit and shaded canopy leaf fractions
PAR_sun = sun_aperlai(1);
PAR_sha = sha_aperlai(1);
    %}
end