function [frac_snow] = RainSnowPartitioningAlptal(T_air)
%{
Since Alptal data only give general precipitation, a function is necessary
to create seperate forcing. This is given by Rutter et al. (2009) as
follows.
!!! In paper, difference (1-T_air)/1.5 results in negative fractions. Thus,
(probably) corrected. !!!
%}

if T_air <= 273.15
    frac_snow = 1;
elseif T_air < 273.15 + 1.5 && T_air > 273.15
    frac_snow = (1.5 - (T_air-273.15))/1.5;
elseif T_air >= 273.15 + 1.5
    frac_snow = 0;
end

end