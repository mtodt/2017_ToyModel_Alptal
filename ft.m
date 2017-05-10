function res = ft(tl,ha)
%{
DESCRIPTION:
photosynthesis temperature response
tl... leaf temperature in photosynthesis temperature function [K]
ha... activation energy in photosynthesis temperature function [J mol^{-1}]
%}

% universal gas constant [J K^{-1} kmole^{-1}]
rgas = 6.02214*10^(26) * 1.38065*10^(-23);

% freezing T of fresh water [K]
tfrz = 273.15;

res = exp(ha/(rgas*10^(-3)*(tfrz+25)) * (1 - (tfrz+25)/tl));