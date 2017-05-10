function [rs,psn] = Stomata(forc_pbot,thm,btran,apar,sla,ei,ea,o2,co2,rb,dayl_factor,tl)
%{
Leaf stomatal resistance and leaf photosynthesis. Modifications for CN code.
%}

tgcm = thm; % air temperature at agcm reference height [K]

% constants
mpe = 10^(-6);
kc25 = 30;
akc = 2.1;
ko25 = 30000;
ako = 1.2;
avcmx = 2.4;
bp = 2000;
act25 = 3.6;
q10act = 2.4;
fnr = 7.16;

% PFT-dependent constants
qe25 = 0.06; % quantum efficiency at 25°C [umol CO2 / umol photon] - 0.06 for all trees and shrubs
leafcn = 40; % leaf C:N (gC/gN) - 40 for evergreen boreal trees, 25 for decidious boreal trees
flnr = 0.04; % fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf) !!! values in Tech Notes and PFT file don't match !!!
fnitr = 0.78;% foliage nitrogen limitation factor - values similar for boreal trees
c3psn = 1;   % photosynthetic pathway: 0. = c4, 1. = c3
mp = 6;      % slope of conductance-to-photosynthesis relationship - 6 for needleleaf, 9 for broadleaf

% convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
act25 = act25*1000/60;

% initialsie rs=rsmax and psn=0 because calculations are performed only
% when apar>0, in which case rs<=rsmax and psn>=0
rsmax0 = 20000;
cf = forc_pbot/(8.3145*10^3 * 0.001*tgcm) * 10^6;
if apar <= 0    % night time
    rs = min(rsmax0,cf/bp);
    psn = 0;
    lnc = 0;
    vcmx = 0;
else
    tc = tl - 273.15;
    ppf = 4.6*apar;
    j = ppf*qe25;
    kc = kc25*(akc^((tc-25)/10));
    ko = ko25*(ako^((tc-25)/10));
    awc = kc*(1+o2/ko);
    cp = 0.5*kc/ko*o2*0.21;
    
    lnc = 1/(sla*leafcn);
    act = act25*(q10act^((tc-25)/10));
    f2 = 1 + exp((-2.2*10^5 + 710*(tc+273.15))/(8.3145*10^3 * 0.001*(tc+273.15)));
    vcmx = lnc*flnr*fnr*act/f2*btran*dayl_factor*fnitr;
    
    % first guess ci
    ci = 0.7*co2*c3psn + 0.4*co2*(1-c3psn);
    
    % rb: s m^{-1} -> s m^2 mol^{-1}
    rb = rb/cf;
    
    % constrain ea
    cea = max(0.25*ei*c3psn + 0.4*ei*(1-c3psn),min(ea,ei));
    
    % ci iteration for 'actual' photosynthesis !!! not NEC_SX !!!
    for iter=1:3
       wj = max(ci-cp,0)*j/(ci+2*cp)*c3psn + j*(1-c3psn); 
       wc = max(ci-cp,0)*vcmx/(ci+awc)*c3psn + vcmx*(1-c3psn); 
       we = 0.5*vcmx*c3psn + 4000*vcmx*ci/forc_pbot * (1-c3psn); 
       psn = min(min(wj,wc),we);
       cs = max(co2 - 1.37*rb*forc_pbot*psn,mpe);
       atmp = mp*psn*forc_pbot*cea/(cs*ei) + bp;
       btmp = (mp*psn*forc_pbot/cs + bp)*rb - 1;
       ctmp = -rb;
       if btmp >= 0
           q = -0.5*(btmp + sqrt(btmp^2 - 4*atmp*ctmp));
       else
           q = -0.5*(btmp - sqrt(btmp^2 - 4*atmp*ctmp));
       end
       r1 = q/atmp;
       r2 = ctmp/q;
       rs = max(r1,r2);
       ci = max(cs - psn*forc_pbot*1.65*rs,0);
    end

    % rs, rb: s m^2 umol^{-1} -> s m^{-1}
    rs = min(rsmax0,rs*cf);
    rb = rb*cf;
end
    
end