function [fval,gs_mol,an] = ci_func(ci,gb_mol,je,cair,oair,lmr_z,par_z,rh_can,tpu_z,vcmax_z,kp_z,bbb,mbb,forc_pbot,t_veg)
%{
DESCRIPTION:
evaluate the function
f(ci) = ci - (ca - (1.37rb+1.65rs))*patm*an
%}

theta_ip = 0.95;    % empirical curvature parameter for ap photosynthesis co-limitation

% activation energy
kcha    = 79430;
koha    = 36380;
cpha    = 37830;
%{
kc, ko, cp
kc25 = 404.9 umol/mol
ko25 = 278.4 mmol/mol
cp25 = 42.75 umol/mol
derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
cp to account for variation in O2 using cp = 0.5 O2 / sco
%}
kc25 = (404.9/(1*10^6)) * forc_pbot;
ko25 = (278.4/(1*10^3)) * forc_pbot;
sco  = 0.5 * 0.209/(42.75/(1*10^6));
cp25 = 0.5 * oair/sco;

kc = kc25 * ft(t_veg,kcha);
ko = ko25 * ft(t_veg,koha);
cp = cp25 * ft(t_veg,cpha);

c3flag = 1; % for C3 PFTs (all (boreal) trees)
if c3flag == 1
    % empirical curvature parameter for ac, aj photosynthesis co-limitation
    theta_cj = 0.98;
    % C3: Rubisco-limited photosynthesis
    ac = vcmax_z * max(ci-cp,0)/(ci + kc*(1+oair/ko));
    % C3: RuBP-limited photosynthesis
    aj = je * max(ci-cp,0)/(4*ci + 8*cp);
    % C3: Product-limited photosynthesis 
    ap = 3*tpu_z;
else
    % quantum efficiency, used only for C4 (mol CO2 / mol photons)
    qe = 0.05;
    % empirical curvature parameter for ac, aj photosynthesis co-limitation
    theta_cj = 0.8;
    % C4: Rubisco-limited photosynthesis
    ac = vcmax_z;
    % C4: RuBP-limited photosynthesis
    aj = qe*par_z*4.6;
    % C4: PEP carboxylase-limited (CO2-limited)
    ap = kp_z*max(ci,0)/forc_pbot;
end

% Gross photosynthesis. First co-limit ac and aj. Then co-limit ap.
aquad = theta_cj;
bquad = -1*(ac + aj);
cquad = ac * aj;
[r1,r2] = quadratic (aquad,bquad,cquad);
ai = min(r1,r2);

aquad = theta_ip;
bquad = -1*(ai + ap);
cquad = ai * ap;
[r1,r2] = quadratic (aquad,bquad,cquad);
ag = min(r1,r2);

% Net photosynthesis. Exit iteration if an < 0.
an = ag - lmr_z;
if an < 0
    fval = 0;
    gs_mol = bbb;   % gs_mol not set to bbb for an<0 until in Photosynthesis.f90, but a value already needed here
else
% Quadratic gs_mol calculation with an known. Valid for an >= 0.
% With an <= 0, then gs_mol = bbb
    cs = cair - 1.4/gb_mol * an*forc_pbot;
    cs = max(cs,1*10^(-6));
    aquad = cs;
    bquad = cs*(gb_mol - bbb) - mbb*an*forc_pbot;
    cquad = -gb_mol*(cs*bbb + mbb*an*forc_pbot*rh_can);
    [r1,r2] = quadratic (aquad,bquad,cquad);
    gs_mol = max(r1,r2);
% derive new estimate for ci
    fval = ci - cair + an * forc_pbot * (1.4*gs_mol + 1.6*gb_mol)/(gb_mol*gs_mol);
end