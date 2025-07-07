function [q_sensible,q_latent,soil_RH] = tima_sensible_and_latent_heat_model(theta_E,windspeed_horiz,air_temp_K,T_surf_K,pressure_air_Pa,RH,soil_VWC,varargin)
% TIMA_SENSIBLE_AND_LATENT_HEAT_MODEL
%   Compute turbulent sensible and latent heat fluxes at the surface, plus
%   an effective soil-relative-humidity metric, using Monin–Obukhov similarity
%   theory (Bonan 2014 §§14.6–14.7) and a soil-moisture–limited evaporation
%   parameterization (Mahfouf & Noilhan 1991; Tran 2016).
%
% Syntax
%   [q_sensible, q_latent, soil_RH] = ...
%       tima_sensible_and_latent_heat_model( ...
%           theta_E, windspeed_horiz, air_temp_K, T_surf_K, ...
%           pressure_air_Pa, RH, soil_VWC, ...
%           'ref_height', h, ...
%           'roughness', z0, ...
%           'displacement', d)
%
% Description
%   Computes:  
%     1. q_sensible — turbulent sensible heat flux [W/m^2]  
%     2. q_latent  — turbulent latent heat flux [W/m^2]  
%     3. soil_RH   — effective near-surface soil relative humidity [0–1]  
%   The method:
%    • Uses aerodynamic resistances from bulk transfer with stability corrections  
%      based on Obukhov length (Bonan 2014).  
%    • Computes latent heat flux limited by soil volumetric water content (VWC).  
%    • Calls tima_psaturationpa_Sonntag1994 for saturation vapor pressures.
%
% Inputs
%   theta_E           [m^3/m^3] — Evaporable soil volumetric water content threshold
%                                  (field capacity); if soil_VWC ≤ theta_E, evaporation
%                                  is reduced by a moisture-dependent factor.
%   windspeed_horiz   [m/s]     — Horizontal wind speed at measurement height.
%   air_temp_K        [K]       — Air (dry‐bulb) temperature at measurement height.
%   T_surf_K          [K]       — Surface (ground) temperature.
%   pressure_air_Pa   [Pa]      — Total air pressure; if empty, a standard
%                                  lapse‐rate estimate is used for density.
%   RH                [0–1]     — Ambient relative humidity (fraction of saturation).
%   soil_VWC          [m^3/m^3] — Soil volumetric water content (top cm).
%
% Name–Value Pair Inputs
%   'ref_height'      [m]       — Reference measurement height; default = 3.
%   'roughness'       [m]       — Surface roughness length z₀; default = 0.01.
%   'displacement'    [m]       — Zero‐plane displacement height; default = 0.
%
% Outputs
%   q_sensible        [W/m^2]   — Sensible heat flux (positive toward the
%                                  atmosphere).
%   q_latent          [W/m^2]   — Latent heat flux (positive toward the
%                                  atmosphere).
%   soil_RH           [0–1]     — Effective soil relative humidity; if
%                                  q_sat_surf ≤ q_air, soil_RH = RH,
%                                  otherwise reduced by soil moisture factor.
%
% Example
%   % compute fluxes for a nearly dry soil under moderate wind
%   [qs, ql, srh] = tima_sensible_and_latent_heat_model( ...
%       0.25,        ... % theta_E
%       2.0,         ... % windspeed_horiz
%       290,         ... % air_temp_K
%       295,         ... % T_surf_K
%       101325,      ... % pressure_air_Pa
%       0.60,        ... % RH
%       0.10,        ... % soil_VWC
%       'ref_height', 5, ...
%       'roughness', 0.02);
%
% References
%   Bonan, G. B. (2014). “Ecological Climatology,” 3rd ed., Chapter 14.
%   Mahfouf, J.-F., & Noilhan, J. (1991). “A Simple Scheme for the
%     Soil–Plant–Atmosphere Continuum,” Q. J. R. Meteorol. Soc.
%   Tran, H. (2016). “Evaporation Parameterization Over Various Soil Moistures.”
%
% See also tima_psaturationpa_Sonntag1994

% Parse name–value pairs
p = inputParser;
p.addParameter('ref_height', 3,   @isnumeric);
p.addParameter('roughness', 0.01, @isnumeric);
p.addParameter('displacement', 0, @isnumeric);
p.parse(varargin{:});
p = p.Results;
displacement=p.displacement; %for unvegetated surfaces
roughness=p.roughness; %Roughness Length (z₀)
ref_height = p.ref_height; % meters, Reference height 

%Some physical constants
vonKarman = 0.4;
gravity = 9.81; % m/s2
Cp_air = 1013;%J/(kg⋅K), Evett in Warrick 2002
R = 8.31446; %J/(K·mol), gas constant
Md = 0.0289652; % kg/mol, Dry air
Mv = 0.018016; % kg/mol, H2O
L_H2O = (2.501-2.37e-3*(T_surf_K-273.15))*10^6;%J/kg, latent heat of evaporation  - Evett in Warrick 2002

% Atmospheric conditions
windspeed_horiz(windspeed_horiz<0.5)=0.5; %m/s, lower limit on measureable wind speed
theta = air_temp_K+0.0098*ref_height; %Potential Temperature of air sensor air at surface
Psat_air = tima_psaturationpa_Sonntag1994(air_temp_K); %Pa
Pv = Psat_air*RH; %Pa, Part Press vapor in air
if isempty(pressure_air_Pa)
    rho_air = 1.291-0.00418*(air_temp_K-273.15); %kg/m^3, Evett in Warrick 2002
else
    Pd = pressure_air_Pa - Pv; %Pa, Part Press dry air
    % Ideal Gas Law
    rho_air = (Pd*Md+Pv*Mv)/(R*air_temp_K); %kg/m^3, accounting for humid air, 
end
Psat_surf = tima_psaturationpa_Sonntag1994(T_surf_K); %Pa
q_air = Mv/Md*Pv/(pressure_air_Pa-(1-Mv/Md)*Pv); %kgh2o/kg dry_air, Specific Humidity Air
q_sat_surf = Mv/Md*Psat_surf/(pressure_air_Pa-(1-Mv/Md)*Psat_surf); %kgh2o/kg dry_air, Specific Humidity Surf (saturation value at surf temp)

% Nondimensional heights from roughness lengths
% Roughness lengths for heat & water ~ 10% of z.0
% (See note in Bonan p 212 following eqn 14.34)
log_zm = log((ref_height-displacement)/roughness);
log_zh = log((ref_height-displacement)/(roughness/10));
log_zw = log_zh;

% Make an initial guess of sensible heat to get the Obukhov Length, then adjust based on stability
CH_guess = vonKarman^2/(log_zm*log_zh); %Neutral flow approximation
q_sensible = -rho_air*Cp_air*windspeed_horiz*CH_guess*(T_surf_K-air_temp_K); %W/m^2

% Friction velocity (from Bonan 14.16)
u_star = vonKarman .* windspeed_horiz ./ log_zm;

% Obukhov Length (Bonan 14.20)
Obukhov_L = -u_star^3 / (vonKarman*gravity/theta * q_sensible / (rho_air*Cp_air));

% Stability parameter (Bonan 14.19)
zeta = (ref_height-displacement) / Obukhov_L; 

%Similarity functions (Bonan 14.26 and 14.27)
if zeta < 0 % unstable case
    x = (1 - 16*zeta)^0.25;
    psi_m = 2 * log((1+x)/2) + log((1+x^2)/2) - 2*atan(x) + pi/2;
    psi_h = 2 * log((1+x^2)/2);
    psi_w = psi_h;
else %stable case
    psi_m = -5 * zeta;
    psi_h = psi_m;
    psi_w = psi_m;
end

% Resistances from Bonan 14.33 
%r_m = (log_zm-psi_m)^2 / (vonKarman^2 * windspeed_horiz);
r_h = (log_zm-psi_m) * (log_zh-psi_h) / (vonKarman^2 * windspeed_horiz);
r_w = (log_zm-psi_m) * (log_zw-psi_w) / (vonKarman^2 * windspeed_horiz);

q_sensible = -rho_air*Cp_air*(T_surf_K-air_temp_K)/r_h; %W/m^2


%% Latent Heat
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
   Evap = rho_air*beta*(q_sat_surf-q_air)/r_w; %kg/m^2, Evaporation rate
end
q_latent = -L_H2O*Evap; %W/m^2, Latent heat flux
end