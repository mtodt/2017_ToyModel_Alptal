function [x,gs_mol,an] = brent(x1,x2,f1,f2,tol,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg,gs_mol,an)
%{
DESCRIPTION:
Use Brent's method to find the root of a single variable function ci_func,
which is known to exist between x1 and x2. The found root will be updated
until its accuracy is tol.

Since gs_mold and an are required to be defined, but only explicitly
calculated under certain conditions, they are passed in from earlier
calculations in hybrid. I'm not sure if that's correct, though.
%}
ITMAX = 20;         % maximum number of iterations
EPS = 1*10^(-2);    % relative error tolerance

a = x1;
b = x2;
fa = f1;
fb = f2;
if (fa > 0 && fb > 0) || (fa < 0 && fb < 0)
    x = nan;
else

c = b;
fc = fb;
iter = 0;

while iter < ITMAX
    iter = iter+1;
    if (fb > 0 && fc > 0) || (fb < 0 && fc < 0)
        c = a;   % Rename a, b, c and adjust bounding interval d.
        fc = fa;
        d = b-a;
        e = d;
    end
    if abs(fc) < abs(fb)
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
    end
    tol1 = 2*EPS*abs(b) + 0.5*tol;  % Convergence check.   
    xm = 0.5*(c-b);
    if abs(xm) <= tol1 || fb == 0
        % x = b;    % happens after the while loop anyway
        break
    end
    if abs(e) >= tol1 && abs(fa) > abs(fb)
        s = fb/fa;  % Attempt inverse quadratic interpolation.
        if a == c
            p = 2*xm*s;
            q = 1-s;
        else
            q = fa/fc;
            r = fb/fc;
            p = s*(2*xm*q*(q-r) - (b-a)*(r-1));
            q = (q-1)*(r-1)*(s-1);
        end
        if p > 0
            q = -q; % Check whether in bounds.
        end
        p = abs(p);
        if 2*p < min(3*xm*q - abs(tol1*q),abs(e*q))
            e = d;  % Accept interpolation.
            d = p/q;
        else
            d = xm; % Interpolation failed, use bisection.
            e = d;
        end
    else % Bounds decreasing too slowly, use bisection.
        d = xm;
        e = d;
    end
    a = b;  % Move last best guess to a.
    fa = fb;
    if abs(d) > tol1   % Evaluate new trial root.
        b = b+d;
    else
        b = b + abs(tol1)*(xm/abs(xm));    % abs(tol1)*(xm/abs(xm)) == sign(tol1,xm) in FORTRAN -.-
    end
    [fb,gs_mol,an] = ci_func(b,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg);
    if fb == 0
        break
    end
end
x = b;
end