function [r1,r2] = quadratic(a,b,c)
%{
DESCRIPTION:
===========================================================================
----------------- Solve quadratic equation for its two roots --------------
===========================================================================
Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
Computing (Cambridge University Press, Cambridge), pp. 145.
%}
if a == 0
    r1 = nan;
    r2 = nan;
else
    if b >= 0
        q = -0.5 * (b + sqrt(b^2 - 4*a*c));
    else
        q = -0.5 * (b - sqrt(b^2 - 4*a*c));
    end
    r1 = q/a;
    if q ~= 0
        r2 = c/q;
    else
        r2 = 1*10^36;
    end
end