function res = fth(tl,hd,se,cc)
%{
DESCRIPTION:
photosynthesis temperature inhibition
tl... leaf temperature in photosynthesis temperature function [K]
hd... deactivation energy in photosynthesis temperature function [J mol^{-1}]
se... entropy term in photosynthesis temperature function [J mol^{-1} K^{-1}]
cc... scaling factor for high temperature inhibition (25 C = 1.0)
%}

% universal gas constant [J K^{-1} kmole^{-1}]
rgas = 6.02214*10^(26) * 1.38065*10^(-23);

res = cc/(1 + exp((-hd+se*tl)/(rgas*10^(-3)*tl)));