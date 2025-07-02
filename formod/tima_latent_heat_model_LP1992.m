function [q_latent,soil_RH,LHCoeff] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_Pa,windspeed_horiz,RH,air_temp_K,T_ground_K,soil_VWC)
%% TIMA_LATENT_HEAT_MODEL_LP1992
%   function to calculate latent heat flux due to evaporation in a given
%   soil layer at a given timestep using a tunable aerodynamic method from
%   Lee & Pielke, 1992
%
% Syntax
%   [q_latent,Soil_RH] = tima_latent_heat_model_LP1992(1000,0.2,1000,5,50,300,290,0.1)
%
% Inputs
%   air_temp_K: [K] near surface air temperature, typically at 3 m AGL.
%           (scalar)
%   CE: [Unitless] resistance to latent heat flux coefficient, similar to the
%       aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
%       (scalar)
%   pressure_air_Pa: [Pa] station pressure, typically at 3m AGL. (scalar)
%   RH: [0-1] array of near surface relative humidity values,
%           typically at 3m AGL. (scalar)
%   soil_VWC: [0-1] soil volumetric water content at depth of interest (scalar)
%   T_ground_K: [K] soil temperature at depth of interest (scalar)
%   theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
%           point (scalar)
%   windspeed_horiz: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL. (scalar)
%
% Outputs
%   q_latent: [W/m^2] latent heat flux from evaporation (scalar)
%   soil_RH: [0-1] Theoretical "relative humidity" within soil (scalar)
%
% Author
%    Ari Koeppel, 2021
%
% Sources
%   Kondo 1990
%   Price 1985
%   http://cires1.colorado.edu/~voemel/vp.html BEST REVIEW OF SATURATION VAPOR PRESSURES
%   Evett in Warrick 2002
%   Sonntag 1994
%   Mahfouf & Noilhan, 1991
%   Tran 2016
%   Lee & Pielke 1992
%   Lee & Pielke 1992: Field capacity- Sand: 0.135, Loam: 0.255, Clay: 0.367, Peat: 0.535
%   Saturation- Sand: 0.395, Loam: 0.451, Clay: 0.482, Peat: 0.863
% Neutral conditions (z = 3 m)
% Smooth water (z0≈10⁻⁴ m): C E ≈ 1.2…1.8×10⁻³
% Flat land (z0≈10⁻² m): C E ≈ 1.0…1.8×10⁻³
% Rough land / short vegetation (z0≈10⁻¹–1 m): C E ≈ 1.6…2.6×10⁻³
% Very rough canopies (z0∼1 m): C E ≈ 2.5…3.5×10⁻³   
% 
% Unstable (convective) stratification
%      C E  typically rises to 2.5…6×10⁻³   (up to ~8×10⁻³ over very rough)   
% Stable stratification
%      C E  typically falls to 0.4…1.2×10⁻³   

    windspeed_horiz(windspeed_horiz<0.5)=0.5; %m/s, lower limit on measureable wind speed
    L_H2O = (2.501-2.37e-3*(T_ground_K-273.15))*10^6;%J/kg, latent heat of evaporation  - Evett in Warrick 2002
    R = 8.31446; %J/(K·mol), gas constant
    Md = 0.0289652; % kg/mol, Dry air
    Mv = 0.018016; % kg/mol, H2O
    Psat_air = tima_psaturationpa_Sonntag1994(air_temp_K); %Pa
    Pv = Psat_air*RH; %Pa, Part Press vapor in air
    if isempty(pressure_air_Pa)
        rho_air = 1.291-0.00418*(air_temp_K-273.15); %kg/m^3, Evett in Warrick 2002
    else
        Pd = pressure_air_Pa - Pv; %Pa, Part Press dry air
        rho_air = (Pd*Md+Pv*Mv)/(R*air_temp_K); %kg/m^3, accounting for humid air
    end
    Psat_surf = tima_psaturationpa_Sonntag1994(T_ground_K); %Pa
    q_air = Mv/Md*Pv/(pressure_air_Pa-(1-Mv/Md)*Pv); %kgh2o/kg dry_air, Specific Humidity Air
    q_sat_surf = Mv/Md*Psat_surf/(pressure_air_Pa-(1-Mv/Md)*Psat_surf); %kgh2o/kg dry_air, Specific Humidity Surf (saturation value at surf temp)
    if soil_VWC <= theta_E % theta_E is max amount a soil will hold if exposed to a constant mist, Mahfouf & Noilhan, 1991, Tran 2016 calls VWC_Sat the evaporation-rate reduction point
        beta = 1/4*(1-cos(soil_VWC/theta_E*pi))^2; % Tran 2016 beta is the resistence to evaporation in top cm resulting from changing water content
    else
        beta = 1; %No resistence when saturated
    end
    if q_sat_surf<=q_air
       soil_RH = RH; 
       Evap=0;
    else
       soil_RH = beta; % Theoretical "relative humidity" within soil, Pct by vol, Lee & Pielke 1992 + Tran, 2016
       Evap = rho_air*beta*windspeed_horiz/CE*(q_sat_surf-q_air); %kg/m^2, Evaporation rate
    end
    q_latent = -L_H2O*Evap; %W/m^2, Latent heat flux
    LHCoeff = -q_latent/(q_sat_surf-q_air);
end
