function Eta = RichardsonToAeta(z_a,T_air,DiffTemp,wind,z_0m,z_0h)
%{
Solves equation eta=Ri*Ri2eta(eta) iteratively.
%}

max_itt = 5;    % prescribe maximum amount of iterations

% Calculate Richardson number
Ri = 9.81*DiffTemp*z_a/(T_air*wind^2);

% STEP 1: Compute function Ri2Eta(Eta)
Eta = 0;
[psi_m_Eta,psi_h_Eta] = StabilityFunctions(Eta);
[psi_m_Oma,psi_h_Oma] = StabilityFunctions(Eta*z_0m/z_a);
[psi_m_Oha,psi_h_Oha] = StabilityFunctions(Eta*z_0h/z_a);
Ri2eta = (log(z_a/z_0m) - psi_m_Eta + psi_m_Oma)*(log(z_a/z_0m) - psi_m_Eta + psi_m_Oma)/...
    (log(z_a/z_0h) - psi_h_Eta + psi_h_Oha);

% STEP 2: Compute error in terms of Ri using etaOld and Ri2eta(etaOld)
Error = Eta/Ri2eta - Ri;

% STEP 3: solve iteratively
acc = 0.0001;
itt=1;
while abs(Error) > acc && itt <= max_itt
		% 3.1 new Eta
    Eta = Ri2eta * Ri;
    divider = (log(z_a/z_0h) - psi_h_Eta + psi_h_Oha);
    if divider ~= 0
        Ri2eta = (log(z_a/z_0m) - psi_m_Eta + psi_m_Oma)*(log(z_a/z_0m) - psi_m_Eta + psi_m_Oma)/divider;
    else
        Ri2eta = 10^12;
    end
		% 3.2 error in terms of Richardson number
    Error = Eta/Ri2eta - Ri;
		% 3.3 update number of iterations
    itt = itt + 1;
end

% STEP 4: Return new eta when convergance criteria is fullfilled
Eta = Ri2eta * Ri;
        % check minimum Monin-Obukhov length (there isthing between 2 to some 100 meters in literature)
if Eta > 0
    L = max(z_a,z_a/Eta);
    Eta = z_a/L;
elseif Eta < 0
	L = max(-1000,z_a/Eta);
    Eta = z_a/L;
end
    
end