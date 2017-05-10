%{
Creating a proxy for the soil humidity from ground water level measurements
at Alptal.

Two time series of ground water level for different soils/soil types are
available. Both of them are used to create a proxy for the ground humidity
necessary for calculations within the CanopyFluxes module.
Therefore, both time series are normalised (using mean and std) and then
treated as normal distributions around 50%.
%}

[Data0304,Variables0304] = xlsread('data_Alptal_SNOWMIP2_200304.xls',1);
GroundWaterLevel10304 = Data0304(:,12);
GroundWaterLevel20304 = Data0304(:,13);

if exist('fig_num','var') == 0
    fig_num = 1;
end

fig=figure(fig_num);fig_num = fig_num+1;
hold on
plot(-GroundWaterLevel10304,'k')
plot(-GroundWaterLevel20304,'r')

% minimum and maximum values determined from time series for both years
Max_GWL1 = 600;
Min_GWL1 = 0;
Max_GWL2 = 1000;
Min_GWL2 = 400;
factor_GWL1_0304 = nan(size(GroundWaterLevel10304)); factor_GWL2_0304 = nan(size(GroundWaterLevel20304));
for t=1:length(GroundWaterLevel10304)
    factor_GWL1_0304(t) = 1 - (GroundWaterLevel10304(t)-Min_GWL1)/(Max_GWL1-Min_GWL1);
    factor_GWL2_0304(t) = 1 - (GroundWaterLevel20304(t)-Min_GWL2)/(Max_GWL2-Min_GWL2);
end

fig=figure(fig_num);fig_num = fig_num+1;
hold on
plot(factor_GWL1_0304,'k')
plot(factor_GWL2_0304,'r')

factor_GWL_0304 = (factor_GWL1_0304 + factor_GWL2_0304)/2;

fig=figure(fig_num);fig_num = fig_num+1;
hold on
plot(factor_GWL_0304,'b')
ylim([0 1])

save('GWL_Proxy_Alptal.mat','')