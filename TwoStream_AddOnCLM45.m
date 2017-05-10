function [abs_veg_dir,abs_veg_dif,alb_veg_dir,alb_veg_dif,trans_veg_dir_dir,trans_veg_dif_dir,trans_veg_dif_dif,omega,G_dir,...
    fabd_sun_z,fabd_sha_z,fabi_sun_z,fabi_sha_z,fsun_z,nrad,tlai_z,vcmaxcintsun,vcmaxcintsha]...
    = TwoStream_AddOnCLM45(coszen,pft,eLAI,eSAI,T_veg,frac_wet,alb_ground_dir,alb_ground_dif,waveband,nlevcan)
%{
DESCRIPTION:
Two-stream fluxes for canopy radiative transfer.
Use two-stream approximation of Dickinson (1983) and Sellers (1985) to
calculate fluxes absorbed by vegetation, reflected by vegetation, and
transmitted through vegetation for unit incoming direct or diffuse flux
given an underlying surface with known albedo.
%}

%------------------  optical parameters for boreal PFTs  ------------------
%---------  1) NET, 2) NDT, 3) BDT, 4) BDS, 5) C3 arctical grass  ---------

% departure of leaf angles from random distribution
orient = [0.01 0.01 0.25 0.25 -0.3];
% reflectances
alpha_leaf(1,:) = [0.07 0.07 0.1 0.1 0.11];     % visual
alpha_leaf(2,:) = [0.35 0.35 0.45 0.45 0.35];   % near-infrared
alpha_stem(1,:) = [0.16 0.16 0.16 0.16 0.31];   % visual
alpha_stem(2,:) = [0.39 0.39 0.39 0.39 0.53];   % near-infrared
% transmittances
tau_leaf(1,:) = [0.05 0.05 0.05 0.05 0.05];     % visual
tau_leaf(2,:) = [0.1 0.1 0.25 0.25 0.34];       % near-infrared
tau_stem(1,:) = [0.001 0.001 0.001 0.001 0.12]; % visual
tau_stem(2,:) = [0.001 0.001 0.001 0.001 0.25]; % near-infrared
% weighted reflectances and transmittances
VAI = eLAI+eSAI;
frac_LAI = eLAI/max(VAI,10^(-6));
frac_SAI = eSAI/max(VAI,10^(-6));

for i=1:2
    alpha(i) = max(alpha_leaf(i,pft)*frac_LAI+alpha_stem(i,pft)*frac_SAI,10^(-6));
    tau(i) = max(tau_leaf(i,pft)*frac_LAI+tau_stem(i,pft)*frac_SAI,10^(-6));
end
% parameters for intercepted snow
omega_snow = [0.8 0.4];
beta_dir_snow = [0.5 0.5];
beta_dif_snow = [0.5 0.5];

% freezing temperature necessary for adjustment for intercepted snow
T_freez = 273.15;

%--------------------------  target definition  ---------------------------

G_dir = 0;
%{
Only calculated in CLM when coszen>0 but needed here later on in any case
within if condition (gdir >0) together with coszen>0 -> gdir not
calculated and therefore still 0 only when coszen<=0 and in that case gdir
in if condition not important anyway.
%}
omega = zeros(1,2); % fraction of intercepted radiation that is scattered (0 to 1)
omega_veg = zeros(1,2);
beta_dir = nan(1,2);    % upscatter parameters for direct radiation
beta_dir_veg = nan(1,2);
beta_dif = nan(1,2);    % upscatter parameters for diffuse radiation
beta_dif_veg = nan(1,2);
abs_veg_dir = zeros(1,2); abs_veg_dif = zeros(1,2);
abs_veg_dir_sun = zeros(1,2); abs_veg_dif_sun = zeros(1,2);
abs_veg_dir_sha = zeros(1,2); abs_veg_dif_sha = zeros(1,2);
alb_veg_dir = zeros(1,2); alb_veg_dif = zeros(1,2);
trans_veg_dir_dir = zeros(1,2); trans_veg_dif_dir = zeros(1,2); trans_veg_dif_dif = zeros(1,2);

%----------------------  CLM4.5 multi-layer canopy  -----------------------
% ncan = 0;
% nrad = 0;
% fabd_sun_z = zeros(length(nlevcan),1);
% fabd_sha_z = zeros(length(nlevcan),1);
% fabi_sun_z = zeros(length(nlevcan),1);
% fabi_sha_z = zeros(length(nlevcan),1);
% fsun_z = zeros(length(nlevcan),1);
% tlai_z = zeros(length(nlevcan),1);
% tsai_z = zeros(length(nlevcan),1);
%{
Diagnose number of canopy layers for radiative transfer, in increments of
dincmax. Add to number of layers so long as cumulative leaf+stem area does
not exceed total leaf+stem area. Then add any remaining leaf+stem area to
next layer and exit the loop. Do this first for elai and esai (not buried
by snow) and then for the part of the canopy that is buried by snow. 
------------------
tlai_z = leaf area increment for a layer
tsai_z = stem area increment for a layer
nrad   = number of canopy layers above snow
ncan   = total number of canopy layers

tlai_z summed from 1 to nrad = elai
tlai_z summed from 1 to ncan = tlai

tsai_z summed from 1 to nrad = esai
tsai_z summed from 1 to ncan = tsai
------------------

Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol" 
because layering is needed for all time steps regardless of radiation

Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered
by nlevcan = 1
%}
dincmax = 0.25;

if nlevcan == 1
    nrad = 1;
    ncan = 1;   % used when considering snow-buried vegetation -> not done here
    tlai_z(1) = eLAI;
    tsai_z(1) = eSAI;
elseif nlevcan > 1
    if eLAI+eSAI == 0
        nrad = 0;
    else
        dincmax_sum = 0;
        for iv = 1:nlevcan
            dincmax_sum = dincmax_sum + dincmax;
            if eLAI+eSAI-dincmax_sum > 10^(-6)
                nrad = iv;
                dinc = dincmax;
                tlai_z(iv) = dinc*eLAI/max(eLAI+eSAI,10^(-6));
                tsai_z(iv) = dinc*eSAI/max(eLAI+eSAI,10^(-6));
            else
                nrad = iv;
                dinc = dincmax - (dincmax_sum - (eLAI+eSAI));
                tlai_z(iv) = dinc*eLAI/max(eLAI+eSAI,10^(-6));
                tsai_z(iv) = dinc*eSAI/max(eLAI+eSAI,10^(-6));
                break
            end
        end
        % mimumum of 4 canopy layers
        if nrad < 4
            nrad = 4;
            for iv = 1:nrad
                tlai_z(iv) = eLAI/nrad;
                tsai_z(iv) = eSAI/nrad;
            end
        end
    end
end

% error check: make sure cumulative of increments does not exceed total
laisum = 0;
saisum = 0;
for iv = 1:nrad
    laisum = laisum + tlai_z(iv);
    saisum = saisum + tsai_z(iv);
end
if abs(laisum-eLAI) > 10^(-6) || abs(saisum-eSAI) > 10^(-6)
    % "endrun" - set everything no nan !!!
end

% now same done for snow-buried part of vegetation - not possible per definition

% zero fluxes for active canopy layers
fabd_sun_z = zeros(nrad,1);
fabd_sha_z = zeros(nrad,1);
fabi_sun_z = zeros(nrad,1);
fabi_sha_z = zeros(nrad,1);
fsun_z = zeros(nrad,1);

%{
Default leaf to canopy scaling coefficients, used when coszen <= 0. This is
the leaf nitrogen profile integrated over the full canopy. Integrate
exp(-kn*x) over x=0 to x=elai and assign to shaded canopy, because sunlit
fraction is 0. Canopy scaling coefficients are set in TwoStream for
coszen > 0. So kn must be set here and in TwoStream.
%}
extkn = 0.30;
if nlevcan == 1
    vcmaxcintsun = 0;
    vcmaxcintsha = (1 - exp(-extkn*eLAI))/extkn;
    if eLAI > 0
        vcmaxcintsha = vcmaxcintsha/eLAI;
    else
        vcmaxcintsha = 0;
    end
elseif nlevcan > 1
    vcmaxcintsun = 0;
    vcmaxcintsha = 0;
end

%--------------------------------------------------------------------------
%-----------------------------  calculations  -----------------------------
%--------------------------------------------------------------------------
if coszen > 0
cosz = max(0.001, coszen);  % cosine of solar zenith angle FOR NEXT TIMESTEP since usually at the end

orient_lim = min(max(orient(pft),-0.4),0.6);    % -0.4 <= orient <= 0.6
if abs(orient_lim) <= 0.01
    orient_lim = 0.01;
end
phi1 = 0.5 - 0.633*orient_lim - 0.33*orient_lim^2;
phi2 = 0.877*(1-2*phi1);
G_dir = phi1 + phi2*cosz; % relative projected area of leaves and stems
twostext = G_dir/cosz;
mu_avg = (1-(phi1/phi2)*log((phi1+phi2)/phi1))/phi2;  % average inverse diffuse optical depth per unit leaf and stem area

temp0 = G_dir + phi2*cosz;
temp1 = phi1*cosz;
temp2 = (1-(temp1/temp0)*log((temp1+temp0)/temp1));

for i=1:waveband
%{
Loop over all wavebands to calculate for the full canopy the scattered fluxes 
reflected upward and transmitted downward by the canopy and the flux absorbed by the 
canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given 
an underlying surface of known albedo.

Output:
------------------
Direct beam fluxes
------------------
albd       - Upward scattered flux above canopy (per unit direct beam flux)
ftid       - Downward scattered flux below canopy (per unit direct beam flux)
ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
fabd       - Flux absorbed by canopy (per unit direct beam flux)
fabd_sun   - Sunlit portion of fabd
fabd_sha   - Shaded portion of fabd
fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
------------------
Diffuse fluxes
------------------
albi       - Upward scattered flux above canopy (per unit diffuse flux)
ftii       - Downward scattered flux below canopy (per unit diffuse flux)
fabi       - Flux absorbed by canopy (per unit diffuse flux)
fabi_sun   - Sunlit portion of fabi
fabi_sha   - Shaded portion of fabi
fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer
%}
    omega_veg(i) = alpha(i)+tau(i);
    asu = 0.5*omega_veg(i)*G_dir/temp0*temp2;
    beta_dir_veg(i) = (1+mu_avg*twostext)/(omega_veg(i)*mu_avg*twostext)*asu;
    beta_dif_veg(i) = ((alpha(i)+tau(i))+(alpha(i)-tau(i))*((1+orient_lim)/2)^2)/(2*omega_veg(i));
    % adjustment for intercepted snow
    if T_veg > T_freez
        tmp0 = omega_veg(i);
        tmp1 = beta_dir_veg(i);
        tmp2 = beta_dif_veg(i);
    else
        tmp0 = (1-frac_wet)*omega_veg(i) + frac_wet*omega_snow(i);
        tmp1 = ((1-frac_wet)*omega_veg(i)*beta_dir_veg(i) + frac_wet*omega_snow(i)*beta_dir_snow(i))/tmp0;
        tmp2 = ((1-frac_wet)*omega_veg(i)*beta_dif_veg(i) + frac_wet*omega_snow(i)*beta_dif_snow(i))/tmp0;
    end
    omega(i) = tmp0;
    beta_dir(i) = tmp1;
    beta_dif(i) = tmp2;

%{
absorbed, reflected, transmitted fluxes per unit incoming radiation
%}
    c1 = omega(i)*beta_dif(i);
    b = 1 - omega(i) + c1;
    tmp0 = mu_avg*twostext;
    d = tmp0 * omega(i)*beta_dir(i);
    f = tmp0 * omega(i)*(1-beta_dir(i));
    tmp1 = b^2 - c1^2;
    h = sqrt(tmp1)/mu_avg;
    sigma = tmp0^2 - tmp1;
    p1 = b + mu_avg*h;
    p2 = b - mu_avg*h;
    p3 = b + tmp0;
    p4 = b - tmp0;

    t1 = min(h*VAI,40); s1 = exp(-t1);          % test to avoid floating...
    t1 = min(twostext*VAI,40); s2 = exp(-t1);   % ... point errors in exp()

%-------------------------  incoming direct flux  -------------------------
    u1 = b - c1/alb_ground_dir(i);
    u2 = b - c1*alb_ground_dir(i);
    u3 = f + c1*alb_ground_dir(i);

    tmp2 = u1 - mu_avg*h;
    tmp3 = u1 + mu_avg*h;
    d1 = p1*tmp2/s1 - p2*tmp3*s1;
    tmp4 = u2 + mu_avg*h;
    tmp5 = u2 - mu_avg*h;
    d2 = tmp4/s1 - tmp5*s1;
    h1 = -d*p4 - c1*f;
    tmp6 = d - h1*p3/sigma;
    tmp7 = (d - c1 - h1/sigma*(u1+tmp0))*s2;
    h2 = (tmp6*tmp2/s1 - p2*tmp7)/d1;
    h3 = -1*(tmp6*tmp3*s1 - p1*tmp7)/d1;
    h4 = -f*p3 - c1*d;
    tmp8 = h4/sigma;
    tmp9 = (u3 - tmp8*(u2-tmp0))*s2;
    h5 = -1*(tmp8*tmp4/s1 + tmp9)/d2;
    h6 = (tmp8*tmp5*s1 + tmp9)/d2;
%     h7 = (c1*tmp2)/(d1*s1);
%     h8 = (-c1*tmp3*s1)/d1;
%     h9 = tmp4/(d2*s1);
%     h10 = (-tmp5*s1)/d2;

% downward direct and diffuse fluxes below vegetation per unit direct flux
    trans_veg_dir_dir(i) = s2;
    trans_veg_dif_dir(i) = h4*s2/sigma + h5*s1 + h6/s1;

% flux reflected by vegetation for direct flux
    alb_veg_dir(i) = h1/sigma + h2 + h3;

% flux absorbed by vegetation per unit direct flux
    abs_veg_dir(i) = 1 - alb_veg_dir(i) - (1-alb_ground_dir(i))*trans_veg_dir_dir(i)...
        - (1-alb_ground_dif(i))*trans_veg_dif_dir(i);

    % additions for sunlit and shaded canopy fractions
    a1 = h1/sigma * (1-s2^2)/(2*twostext) + h2*(1-s2*s1)/(twostext+h) + h3*(1-s2/s1)/(twostext-h);
    a2 = h4/sigma * (1-s2^2)/(2*twostext) + h5*(1-s2*s1)/(twostext+h) + h6*(1-s2/s1)/(twostext-h);
    abs_veg_dir_sun(i) = (1-omega(i)) * (1-s2+1/mu_avg*(a1+a2));
    abs_veg_dir_sha(i) = abs_veg_dir(i) - abs_veg_dir_sun(i);

%-----------------------  incoming diffusive flux  ------------------------
    u1 = b - c1/alb_ground_dif(i);
    u2 = b - c1*alb_ground_dif(i);
%    u3 = f + c1*alb_ground_dif(i);

    tmp2 = u1 - mu_avg*h;
    tmp3 = u1 + mu_avg*h;
    d1 = p1*tmp2/s1 - p2*tmp3*s1;
    tmp4 = u2 + mu_avg*h;
    tmp5 = u2 - mu_avg*h;
    d2 = tmp4/s1 - tmp5*s1;
%    h1 = -d*p4 - c1*f;
%    tmp6 = d - h1*p3/sigma;
%    tmp7 = (d - c1 - h1/sigma*(u1+tmp0))*s2;
%    h2 = (tmp6*tmp2/s1 - p2*tmp7)/d1;
%    h3 = -1*(tmp6*tmp3*s1 - p1*tmp7)/d1;
%    h4 = -f*p3 - c1*d;
%    tmp8 = h4/sigma;
%    tmp9 = (u3 - tmp8*(u2-tmp0))*s2;
%    h5 = -1*(tmp8*tmp4/s1 + tmp9)/d2;
%    h6 = (tmp8*tmp5*s1 + tmp9)/d2;
    h7 = (c1*tmp2)/(d1*s1);
    h8 = (-c1*tmp3*s1)/d1;
    h9 = tmp4/(d2*s1);
    h10 = (-tmp5*s1)/d2;

% downward direct and diffuse fluxes below vegetation per unit diffuse flux
    trans_veg_dif_dif(i) = h9*s1 + h10/s1;

% Flux reflected by vegetation
    alb_veg_dif(i) = h7 + h8;

% Flux absorbed by vegetation
    abs_veg_dif(i) = 1 - alb_veg_dif(i) - (1-alb_ground_dif(i))*trans_veg_dif_dif(i);
    
    % additions for sunlit and shaded canopy fractions
    a1 = h7*(1-s2*s1)/(twostext+h) + h8*(1-s2/s1)/(twostext-h);
    a2 = h9*(1-s2*s1)/(twostext+h) + h10*(1-s2/s1)/(twostext-h);
    abs_veg_dif_sun(i) = (1-omega(i))/mu_avg*(a1+a2);
    abs_veg_dif_sha(i) = abs_veg_dif(i) - abs_veg_dif_sun(i);
  
%{
Repeat two-stream calculations for each canopy layer to calculate derivatives. 
tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives
are calculated at the center of the layer. Derivatives are needed only for
the visible waveband to calculate absorbed PAR (per unit lai+sai) for each
canopy layer. Derivatives are calculated first per unit lai+sai and then
normalized for sunlit or shaded fraction of canopy layer.
Sun/shade big leaf code uses only one layer, with canopy integrated values
from above and also canopy-integrated scaling coefficients.
%}
    if i==1
        if nlevcan == 1
            % sunlit fraction of canopy
            fsun_z(1) = (1-s2)/t1;
            
            % absorbed PAR (per unit sun/shade lai+sai)
            laisum = eLAI+eSAI;
            fabd_sun_z(1) = abs_veg_dir_sun(i)/(fsun_z(1)*laisum);
            fabi_sun_z(1) = abs_veg_dif_sun(i)/(fsun_z(1)*laisum);
            fabd_sha_z(1) = abs_veg_dir_sha(i)/((1 - fsun_z(1))*laisum);
            fabi_sha_z(1) = abs_veg_dif_sha(i)/((1 - fsun_z(1))*laisum);
            
            % leaf to canopy scaling coefficients
            extkn = 0.3;
            extkb = twostext;
            vcmaxcintsun = (1 - exp(-(extkn+extkb)*eLAI))/(extkn + extkb);
            vcmaxcintsha = (1 - exp(-extkn*eLAI))/extkn - vcmaxcintsun;
            if eLAI > 0
                vcmaxcintsun = vcmaxcintsun/(fsun_z(1)*eLAI);
                vcmaxcintsha = vcmaxcintsha/((1 - fsun_z(1))*eLAI);
            else
                vcmaxcintsun = 0;
                vcmaxcintsha = 0;
            end
        elseif nlevcan > 1
            for iv = 1:nrad
                % cumulative lai+sai at center of layer
                if iv == 1
                   laisum = 0.5*(tlai_z(iv)+tsai_z(iv));
                else
                   laisum = laisum + 0.5*((tlai_z(iv-1)+tsai_z(iv-1))+(tlai_z(iv)+tsai_z(iv)));
                end
                
                % coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction
                t1 = min(h*laisum,40); s1 = exp(-t1);
                t1 = min(twostext*laisum,40); s2 = exp(-t1);
                fsun_z(iv) = s2;
                
                % ===============
                % Direct beam
                % ===============
                % coefficients h1-h6 and a1,a2 depend of cumulative lai+sai
                u1 = b - c1/alb_ground_dir(i);
                u2 = b - c1*alb_ground_dir(i);
                u3 = f + c1*alb_ground_dir(i);
                tmp2 = u1 - mu_avg*h;
                tmp3 = u1 + mu_avg*h;
                d1 = p1*tmp2/s1 - p2*tmp3*s1;
                tmp4 = u2 + mu_avg*h;
                tmp5 = u2 - mu_avg*h;
                d2 = tmp4/s1 - tmp5*s1;
                h1 = -d*p4 - c1*f;
                tmp6 = d - h1*p3/sigma;
                tmp7 = (d - c1 - h1/sigma*(u1+tmp0))*s2;
                h2 = (tmp6*tmp2/s1 - p2*tmp7)/d1;
                h3 = -1*(tmp6*tmp3*s1 - p1*tmp7)/d1;
                h4 = -f*p3 - c1*d;
                tmp8 = h4/sigma;
                tmp9 = (u3 - tmp8*(u2-tmp0))*s2;
                h5 = -1*(tmp8*tmp4/s1 + tmp9)/d2;
                h6 = (tmp8*tmp5*s1 + tmp9)/d2;
                
                a1 = h1/sigma*(1-s2^2)/(2*twostext) + ...
                    h2*(1-s2*s1)/(twostext+h) + h3*(1-s2/s1)/(twostext-h);

                a2 = h4/sigma*(1-s2^2)/(2*twostext) + ...
                    h5*(1-s2*s1)/(twostext+h) + h6*(1-s2/s1)/(twostext-h);

                % derivatives for h2, h3, h5, h6 and a1, a2
                v = d1;
                dv = h*p1*tmp2/s1 + h*p2*tmp3*s1;
                u = tmp6*tmp2/s1 - p2*tmp7;
                du = h*tmp6*tmp2/s1 + twostext*p2*tmp7;
                dh2 = (v*du - u*dv)/(v^2);
                u = -tmp6*tmp3*s1 + p1*tmp7;
                du = h*tmp6*tmp3*s1 - twostext*p1*tmp7;
                dh3 = (v*du - u*dv)/(v^2);
                v = d2;
                dv = h*tmp4/s1 + h*tmp5*s1;
                u = -h4/sigma * tmp4/s1 - tmp9;
                du = -h*h4/sigma * tmp4/s1 + twostext*tmp9;
                dh5 = (v*du - u*dv)/(v^2);
                u = h4/sigma * tmp5*s1 + tmp9;
                du = -h*h4/sigma * tmp5*s1 - twostext*tmp9;
                dh6 = (v*du - u*dv)/(v^2);
                da1 = h1/sigma * s2^2 + h2*s2*s1 + h3*s2/s1 ...
                    + (1-s2*s1)/(twostext+h)*dh2 + (1-s2/s1)/(twostext-h)*dh3;
                da2 = h4/sigma * s2^2 + h5*s2*s1 + h6*s2/s1 ...
                    + (1-s2*s1)/(twostext+h)*dh5 + (1-s2/s1)/(twostext-h)*dh6;
                
                % flux derivatives
                d_ftid = -twostext*h4/sigma*s2 - h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1;
                d_fabd = -1*(dh2+dh3) + (1-alb_ground_dir(i))*twostext*s2 - (1-alb_ground_dif(i))*d_ftid;
                d_fabd_sun = (1-omega(i))*(twostext*s2 + 1/mu_avg*(da1+da2));
                d_fabd_sha = d_fabd - d_fabd_sun;
                
                fabd_sun_z(iv) = max(d_fabd_sun,0);
                fabd_sha_z(iv) = max(d_fabd_sha,0);
                
                %{
                Flux derivatives are APARsun and APARsha per unit (LAI+SAI).
                Need to normalize derivatives by sunlit or shaded fraction
                to get APARsun per unit (LAI+SAI)sun and APARsha per unit
                (LAI+SAI)sha.
                %}
                fabd_sun_z(iv) = fabd_sun_z(iv)/fsun_z(iv);
                fabd_sha_z(iv) = fabd_sha_z(iv)/(1 - fsun_z(iv));
                
                % ===============
                % Diffuse
                % ===============
                % coefficients h7-h10 and a1,a2 depend of cumulative lai+sai
                u1 = b - c1/alb_ground_dif(i);
                u2 = b - c1*alb_ground_dif(i);
                tmp2 = u1 - mu_avg*h;
                tmp3 = u1 + mu_avg*h;
                d1 = p1*tmp2/s1 - p2*tmp3*s1;
                tmp4 = u2 + mu_avg*h;
                tmp5 = u2 - mu_avg*h;
                d2 = tmp4/s1 - tmp5*s1;
                h7 = (c1*tmp2)/(d1*s1);
                h8 = (-c1*tmp3*s1)/d1;
                h9 = tmp4/(d2*s1);
                h10 = (-tmp5*s1)/d2;
                
                a1 = h7*(1-s2*s1)/(twostext+h) + h8*(1-s2/s1)/(twostext-h);
                a2 = h9*(1-s2*s1)/(twostext+h) + h10*(1-s2/s1)/(twostext-h);
                
                % derivatives for h7, h8, h9, h10 and a1, a2
                v = d1;
                dv = h*p1*tmp2/s1 + h*p2*tmp3*s1;
                u = c1*tmp2/s1;
                du = h*c1*tmp2/s1;
                dh7 = (v*du - u*dv)/(v^2);
                u = -c1*tmp3*s1;
                du = h*c1*tmp3*s1;
                dh8 = (v*du - u*dv)/(v^2);
                v = d2;
                dv = h*tmp4/s1 + h*tmp5*s1;
                u = tmp4/s1;
                du = h*tmp4/s1;
                dh9 = (v*du - u*dv)/(v^2);
                u = -tmp5*s1;
                du = h*tmp5*s1;
                dh10 = (v*du - u*dv)/(v^2);
                
                da1 = h7*s2*s1 +  h8*s2/s1 + (1-s2*s1)/(twostext+h)*dh7 + (1-s2/s1)/(twostext-h)*dh8;
                da2 = h9*s2*s1 + h10*s2/s1 + (1-s2*s1)/(twostext+h)*dh9 + (1-s2/s1)/(twostext-h)*dh10;
                
                % flux derivatives
                d_ftii = -h*h9*s1 + h*h10/s1 + dh9*s1 + dh10/s1;
                d_fabi = -1*(dh7+dh8) - (1-alb_ground_dif(i))*d_ftii;
                d_fabi_sun = (1-omega(i))/mu_avg * (da1+da2);
                d_fabi_sha = d_fabi - d_fabi_sun;
                
                fabi_sun_z(iv) = max(d_fabi_sun,0);
                fabi_sha_z(iv) = max(d_fabi_sha,0);
                
                %{
                Flux derivatives are APARsun and APARsha per unit (LAI+SAI).
                Need to normalize derivatives by sunlit or shaded fraction
                to get APARsun per unit (LAI+SAI)sun and APARsha per unit
                (LAI+SAI)sha.
                %}
                fabi_sun_z(iv) = fabi_sun_z(iv)/fsun_z(iv);
                fabi_sha_z(iv) = fabi_sha_z(iv)/(1-fsun_z(iv));
            end
        end
    end
end
end