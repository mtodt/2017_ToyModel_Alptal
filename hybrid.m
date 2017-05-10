function [x0,gs_mol,an] = hybrid(x0,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg)
%{
Use a hybrid solver to find the root of equation
f(x) = x - h(x),
we want to find x, s.t. f(x) = 0.
The hybrid approach combines the strength of the newton secant approach
(find the solution domain) and the bisection approach implemented with the
Brent's method to guarrantee convergence.
%}
eps = 1*10^(-2);     % relative accuracy
eps1 = 1*10^(-4);
itmax = 40;          % maximum number of iterations

[f0,gs_mol,an] = ci_func(x0,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg);

if f0 == 0
   x0 = x0;
else
    minx = x0;
    minf = f0;
    x1 = x0*0.99;
    [f1,gs_mol,an] = ci_func(x1,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg);
    
    if f1 == 0
        x0 = x1;
    else
        if f1 < minf
            minx = x1;
            minf = f1;
        end
    
        % first use the secant approach, then use the brent approach as a backup
        iter = 0;
        while iter <= itmax
            iter = iter+1;
            dx = -f1 * (x1-x0)/(f1-f0);
            x = x1 + dx;
            tol = abs(x)*eps;
            if abs(dx) < tol
                x0 = x;
                break
            end
            
            x0 = x1;
            f0 = f1;
            x1 = x;
            [f1,gs_mol,an] = ci_func(x1,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg);
        
            if f1 < minf
                minx = x1;
                minf = f1;
            end
            if abs(f1) <= eps1
                x0 = x1;
                break
            end
        
            % if a root zone is found, use the brent method for a robust backup strategy
            if f1*f0 < 0
                [x,gs_mol,an] = brent(x0,x1,f0,f1,tol,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg,gs_mol,an);
                x0 = x;
                break
            end
        
            if iter > itmax
            %{
            In case of failing to converge within itmax iterations stop at the
            minimum function. This happens because of some other issues besides
            the stomatal conductance calculation and it happens usually in very
            dry places and more likely with c4 plants.
            %}
                [f1,gs_mol,an] = ci_func(minx,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg);
                break
            end
        end
    end
end