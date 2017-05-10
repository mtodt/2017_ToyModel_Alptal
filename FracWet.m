function [frac_wet,frac_dry] = FracWet(eLAI,eSAI,frac_veg_nosno,int_stor)
%{
DESCRIPTION:
Determine fraction of vegetated surfaces which are wet and fraction of elai
which is dry. The variable "fwet" is the fraction of all vegetation surfaces
which are wet including stem area which contribute to evaporation. The variable
"fdry" is the fraction of elai which is dry because only leaves can transpire.
Adjusted for stem area which does not transpire.
%}

dewmx = 0.1;    % maximum storage of water [kg m^{-2}] or maximum allowed dew [mm]

if frac_veg_nosno > 0
    if int_stor > 0
        frac_wet = ((int_stor/((eLAI+eSAI)*dewmx)))^(2/3);
        frac_wet = min(frac_wet,1);
    else
        frac_wet = 0;
    end
    frac_dry = (1-frac_wet)*eLAI/(eLAI+eSAI);
elseif frac_veg_nosno == 0
    frac_wet = 0;
    frac_dry = 0;
end

end