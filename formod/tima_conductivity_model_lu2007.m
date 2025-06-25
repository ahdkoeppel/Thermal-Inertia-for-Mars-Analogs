function [k_eff] = tima_conductivity_model_lu2007(k_dry_std,Soil_Temperature,T_std,VWC,theta_k,m,Soil_RH,material)
%% TIMA_CONDUCTIVITY_MODEL_LU2007
%   function to calculate the effective thermal conductivity of a particulate soil or sediment using the Lu 2007 Method
%
% Description
%   [m = ~0.3 seems best] 0<m<1.33  Best for granular with no hydration regime, Realistic lower limit
%      The lu 2007 fitted values of m are 0.96 and 0.27 for the coarse-textured and fine textured soils, respectively.
%
% Syntax
%   [k_eff] = tima_conductivity_model_lu2007(k_dry_std,Soil_Temperature,T_std,VWC,theta_k,m,Soil_RH,material)
%
% Inputs
%   k_dry_std: [W/mK] stardard bulk dry thermal conductivty at T_std (vector)
%   Soil_Temperature: [K] Soil temperature (vector)
%   T_std: [K] Standard temperature; typically 300 (vector)
%   VWC: [fraction by volume] volumetric water content (vector)
%   k_H2O: [W/mK] Thermal conductivity of soil water (vector)
%   theta_k: [fraction by volume] conductivity soil moisture inflection point - in theory this should be similar to saturation value or porosity (vector)
%   m: [unitless] pore-network-connectivity-parameter (vector)
%   Soil_RH: [fraction] relative humidity of air within the soil (vector)
%   material: ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%
% Outputs
%   k_eff: [W/mK] effective thermal conductivity (vector)
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Zhang and Wang 2017, 
%   Dong 2015
%   Lu and Dong 2015
%   Tsilingiris 2008
%   Piqueux 2009a
%   Piqueux and Christensen 2011
%   Clauser and Huenges [1995]
%   Bristow, 2002: Basalt: 2.2, Granite: 2.0, Quartz: 8.8, Clay: 2.9, organics: 0.25, Ice @ 0C: 2.18
%   Water = 0.552+2.34E-3*Soil_Temperature_C-1.1E-5*Soil_Temperature_C^2
%   Air: 0.0237 + 0.000064*Soil_Temperature_C
%   Campbell 1994: Dry_air = 0.024+7.73E-5.*(Soil_Temperature-273.15) - 2.6E-8*(Soil_Temperature-273.15)^2
%       k_H2O = 0.554+2.23E-3.*(Soil_Temperature - 273.15) - 9.87E-6.*(Soil_Temperature - 273.15).^2; %Campbell 1994
%   Horai 1971
%   Pielke 2002; k_eff Concrete: 4.6, Rock: 2.93, Ice: 2.51, Snow: 0.08-1.67, Stable air: 0.02-0.03, water; 0.57-0.63
%       Clay-dry: 0.25, clay-10%h2o: 0.63, clay-20%h2o: 1.12, clay-30%h2o: 1.33, clay-40%h2o: 1.58, --porosity 40%
%       sand-dry: 0.30, sand-10%h2o: 1.05, sand-20%h2o: 1.95, sand-30%h2o: 2.16, sand-40%h2o: 2.20, --porosity 40%
%       peat-dry: 0.06, peat-10%h2o: 0.10, peat-20%h2o: 0.29, peat-30%h2o: 0.43, peat-40%h2o: 0.50, --porosity 80%
%       rooty-soil: 0.11;
%   Lee & Pielke 1992: Field capacity- Sand: 0.135, Loam: 0.255, Clay: 0.367, Peat: 0.535
%       Saturation- Sand: 0.395, Loam: 0.451, Clay: 0.482, Peat: 0.863

% Air
k_air = 0.024+7.78E-5.*(Soil_Temperature-273.15)-3E-05.*exp(0.0591.*(Soil_Temperature-273.15)).*Soil_RH; %Bristow 2002/Campbell 1994 with Tsilingiris 2008 Fig 3 assuming soil is equilibrated with air
k_air_avg = 0.024+7.78E-5.*(T_std-273.15); %Thermal conductivity of air at zero moisture and standard temperature
k_dry_mod = 2*k_air.^(0.8964.*theta_k + 0.28); % Piquex 2009a Fig 8 - independent of composition and grain size
k_dry_mod_avg = 2*k_air_avg.^(0.8964.*theta_k + 0.28); % Piquex 2009a Fig 8 R > 0.98 within range 0.01-0.035 - independent of composition and grain size
k_dry_mod(k_dry_mod<k_air) = k_air(k_dry_mod<k_air);
k_dry_mod_avg(k_dry_mod_avg<k_air_avg) = k_air_avg(k_dry_mod_avg<k_air_avg);
k_dry = k_dry_std.*k_dry_mod./k_dry_mod_avg; % Since Piquex 2009a Fig 8 was modeled for a specific material type, k_dry_mod does not encessarily apply here, but the ratio of k_dry_mod/k_dry_mod_avg should be similar.

k_H2O = 0.61;%-1.1e-5.*(Soil_Temperature - 273.15).^2 + 0.00234.*(Soil_Temperature - 273.15) + 0.552; %Bristow 2002 or 0.59 from Zhang and Wang 2017 or 0.6096 from Haigh 2012
% k-bulk-dry from theory
if strcmp(material,"basalt")
    k_solid = 1.18 + 474./(350+Soil_Temperature-273.15);%  % Piqueux and Christensen 2011/Clauser and Huenges [1995]
% Roughly = 2.2, Bristow, 2002
elseif strcmp(material,"amorphous")
    k_solid = 0.718 + 0.0015.*(Soil_Temperature); %Piqueux and Christensen 2011/Wood2020/Clauser and Huenges [1995]/Siegler 2012
elseif strcmp(material,"granite")
    k_solid = 2.0; %Wood 2020
    % Roughly = 2.0 (granite), Bristow, 2002
elseif strcmp(material,"sandstone")
    k_solid =  0.18858+2192.2./Soil_Temperature; %Quartz Wood 2020/Birch and Clark 1940
    %Roughly = 7.69; %(horai 1971)
elseif strcmp(material,"clay")
    k_solid = 2.9;
    % Roughly = 2.9 (clay), Bristow, 2002
elseif strcmp(material,"salt")
    k_solid = -2.11 + 2960./(350+Soil_Temperature-273.15); %Piqueux and Christensen 2011/Clauser and Huenges [1995]
elseif strcmp(material,"ice")
    k_solid = 632./Soil_Temperature +0.38-0.00197.*Soil_Temperature; %Wood 2020/Andersson and Inaba 2005
    % Roughly = 2.18, Bristow, 2002
else
    k_solid = 2.2;
end

Sr = VWC./theta_k; %degree of field capacity
k_wet = k_solid.^(1-theta_k).*(k_H2O.^theta_k).*k_dry;%geometric mean to 
% calculate thermal conductivity at field capacity From Johansen 1975 &
% Cote and Konrad 2005

if k_wet(1) < k_dry(1) % In case theoretical k_wet is smaller than k_dry (shouldn't happen because
%water is replacing air)
    k_wet = k_dry-k_air.*theta_k + k_H2O.*theta_k; %Parallel
end
Sr(Sr>1)=1;

%**********K Model************
k_eff = k_dry + exp(m.*(1-(Sr).^(m-1.33))).*(k_wet-k_dry); %Lu 2007

end