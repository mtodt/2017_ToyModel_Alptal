# 2017_ToyModel_Alptal

All .m-files in this directory were used to analyse data from the forest site at Alptal, Switzerland, and create the figures displayed in "Simulation of longwave enhancement in a needleleaf forest" by Todt et. al (submitted in May 2017). There is only one file required to reproduce the results within our paper: 'data_Alptal_SNOWMIP2_200304.xls' which can be obtained from Tobias Jonas at SLF, Switzerland, upon request (jonas@slf.ch).

Before running the Toy Model it is necessary to create a proxy for soil moisture, for which the script '???.m' is used. Afterwards running script 'ToyModel_Alptal_ValidationPeriodRun.m' will produce the output. Note that the model evaluation period is limited by the availability of sub-canopy radiation data - hence the file name - although observations required to actually run the model start in late August of 2003. Parameter 'HM_in_CLM' within 'ToyModel_Alptal_ValidationPeriodRun.m' is used to turn the conductive heat flux parameterisation within the CLM4.5 module on or off. 

To reprodcue the figures used within our paper the following steps are necessary:
1) Run '???.m' to create a proxy for soil moisture.
2) Run 'ToyModel_Alptal_ValidationPeriodRun.m' with the conductive heat flux parameterisation turned on.
3) Save variable 'LW_in_bc_CLM' as 'LW_in_bc_CLM_HM' in 'LWsub_HM.mat'.
4) Run 'ToyModel_Alptal_ValidationPeriodRun.m' with the conductive heat flux parameterisation turned off.

To reproduce the results of the sensitivity studies the following steps are necessary:
1) Run 'ToyModel_Alptal_ValidationPeriodRunSensitivity.m' with the conductive heat flux parameterisation turned off multiple times while changing the parameter that is variied in the last loop within the script (index j).
2) Save 'LW_in_bc_CLM' and 'LW_in_bc_SP' for every iteration as '???.mat'.
3) Run '???.m' to produce RMSE and mean bias values as well as the figure of longwave enhancement PDFs.

The following scripts are used for CLM4.5 calculations:
Biogeophysics1.m
CanopyFluxes.m
CanopyFluxes_HeatMass.m
CanopyHydrology45.m
CosZenCalculation.m
DSaturationPressureDT.m
FracWet.m
FrictionVelocity.m
Photosynthesis.m
QSat.m
RichardsonToAeta.m
StabFunc1.m
StabFunc2.m
Stomata.m
SurfaceRadiation.m
TwoStream_AddOnCLM45.m
brent.m
ci_func.m
ft.m
fth.m
fth25.m
hybrid.m
quadratic.m

The following scripts are used for SNOWPACK calculations:
CanopyRadiationOutput.m
NetRadiation2L.m
SNOWPACK2L_EnergyFluxes_withinCLM45.m
SolarElevationAngle.m
StabilityFunctions.m
TurbulentExchange.m
WaterSaturationPressure.m

The following scripts are of general use:
RainSnowPartitioningAlptal.m
SNOWPACKwithinCLM45driver.m
ToyModel_Alptal_ValidationPeriodRun.m
