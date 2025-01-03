function [Cp] = tima_specific_heat_model_hillel(rho_dry,rho)
% TIMA_SPECIFIC_HEAT_MODEL_HILLEL
%   function to calculates specific heat as a funtion of water content and
%   density in J/kgK
%
% Description
%   A Cp model that considers water content and yields similar values as
%   Verhoef model (considering temperature) for mineral soils while saving computation time.
%
% Syntax
%   [Cp] = tima_specific_heat_model_hillel(1300,1600)
%
% Inputs
%   rho_dry: [kg/m^3] dry soil density, for ice set to 950 (scalar)
%   rho: [kg/m^3] current soil density, for ice set to 950 (scalar)
%
% Outputs
%   Cp: [J/kgK] Specific heat capacity  (scalar) 
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Hillel 1980 + Evett in Warrick 2002
%   Pielke 2002; Concrete: 879, Rock: 753, Ice: 2100, Snow: 2093, Stable air: 1005, water; 4186
%   Clay-dry: 890, clay-10%h2o: 1005, clay-20%h2o: 1172, clay-30%h2o: 1340, clay-40%h2o: 1550,  --porosity 40%
%   sand-dry: 800, sand-10%h2o: 1088, sand-20%h2o: 1256, sand-30%h2o: 1423, sand-40%h2o: 1480,  --porosity 40%
%   peat-dry: 1920, peat-10%h2o: 2302, peat-20%h2o: 3098, peat-30%h2o: 3433, peat-40%h2o: 3650, --porosity 80%
%   rooty-soil: 1256;
%   Abu-Hamdeh -- clay: 1170-2250, sand: 830-1670
%   Hanks 1992: Good temperature approximations can be made...even for many nonuniform soils by assuming a uniform thermal diffusivity.

%Assumes Cp_solid == 753;

Cp_H2O = 4186; %J/kgK
Cp = (753.*rho_dry+Cp_H2O.*(rho-rho_dry))./rho; %J/kgK
end