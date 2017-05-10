function [int_stormax,int_stor,int,through_snow,through_rain,drip_snow,drip_rain,ground_snow,ground_rain,frac_wet]...
    = CanopyHydrology45(eLAI,eSAI,frac_veg_nosno,forc_rain,forc_snow,int_stor_pre,dtime)
%{
DESCRIPTION:
Calculation of water storage of intercepted precipitation, direct throughfall
and canopy drainage of precipitation.

Difference between model code and tech notes:
Oleson et al. (2013) -> [kg m^{-2}]
model code -> [mm]
Density of water used to scale, so values are the same.
%}

% definition of parameters
alpha = 0.25;   % fraction of potential interception "scales from point to grid cell" (Oleson et al., 2013)
dewmx = 0.1;    % maximum storage of water [kg m^{-2}] OR "maximum allowed dew [mm]" as in CLM code ??? - equal for density of water
        
fpi = alpha*(1-exp(-0.5*(eLAI+eSAI)));    % coefficient of interception
int_stormax = dewmx*(eLAI+eSAI);     % maximum allowed water on canopy [mm]
%{
The leaf water capacities for solid and liquid are different, generally double for snow,
but these are of somewhat less significance for the water budget because of lower
evaporation rate at lower temperature.  Hence, it is reasonable to assume that vegetation
storage of solid water is the same as liquid water.
%}

if frac_veg_nosno > 0 && (forc_rain+forc_snow) > 0
    frac_snow = forc_snow/(forc_snow + forc_rain);
    frac_rain = forc_rain/(forc_snow + forc_rain);

    through_snow = forc_snow*(1-fpi);
    through_rain = forc_rain*(1-fpi);

    int = (forc_snow+forc_rain)*fpi;   % intercepted precipitation [mm/s]
    int_stor = max(0,int_stor_pre+dtime*int);     % water storage of intercepted precipitation and dew

    xrun = (int_stor-int_stormax)/dtime;    % excess water that exceeds the leaf capacity
    drip = 0;
    if xrun > 0
       drip = xrun;
       int_stor = int_stormax;
    end
    drip_snow = drip*frac_snow;
    drip_rain = drip*frac_rain;
else
    frac_snow = 0; frac_rain = 0;
    through_snow = 0; through_rain = 0;
    int = 0;
    int_stor = int_stor_pre;
    drip_snow = 0; drip_rain = 0;
end
if frac_veg_nosno == 0
    ground_snow = forc_snow;
    ground_rain = forc_rain;
else
    ground_snow = through_snow + drip_snow;
    ground_rain = through_rain + drip_rain;    % plus additional mass balance term in model
end

[frac_wet,frac_dry] = FracWet(eLAI,eSAI,frac_veg_nosno,int_stor);

end