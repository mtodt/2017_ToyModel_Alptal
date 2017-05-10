function P_sat = WaterSaturationPressure(T)
% Standard water vapor saturation pressure [Pa].

t_water_triple_pt = 273.16;     % triple point temperature for water [K]
p_water_triple_pt = 611.73;     % triple point partial vapor pressure for water [Pa]

if T < t_water_triple_pt        % for a flat ice surface
		c2 = 21.88;
		c3 = 7.66;
else                            % for a flat water surface
		c2 = 17.27;
		c3 = 35.86;
end

exp_p_sat = c2*(T - t_water_triple_pt)/(T - c3);

P_sat = p_water_triple_pt*exp(exp_p_sat);




end