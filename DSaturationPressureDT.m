function dpdt = DSaturationPressureDT(L,T)
% Temperature derivative of the saturation pressure function.

lh_subl = 2.838*10^6;        % latent heat of sublimation [J kg^{-1}]

if L ~= lh_subl
    c2 = 17.08085;
    c3 = 234.175;
else
    c2 = 22.44294;
    c3 = 272.44;
end

P_sat = WaterSaturationPressure(T);

dpdt = P_sat*c2*c3/((c3+(T-273.15))^2);

end