function res = fth25(hd,se)
%{
DESCRIPTION:
scaling factor for photosynthesis temperature inhibition
hd... deactivation energy in photosynthesis temperature function [J mol^{-1}]
se... entropy term in photosynthesis temperature function [J mol^{-1} K^{-1}]
%}

% universal gas constant [J K^{-1} kmole^{-1}]
rgas = 6.02214*10^(26) * 1.38065*10^(-23);

% freezing T of fresh water [K]
tfrz = 273.15;

res = 1 + exp((-hd+se*(tfrz+25))/(rgas*10^(-3)*(tfrz+25)));