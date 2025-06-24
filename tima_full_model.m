function [T_surf_C] = tima_full_model(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,T_Start,layer_size,VWC_column,evap_depth,RH,emissivity,...
    pressure_air_pa,T_surf_obs_C,NDAYS,varargin)
%% TIMA_FULL_MODEL
%   This function runs both the subsurface temperature initialization for a
%   given set of coefficents and the subsequent heat transfer forward
%   model.

% Syntax
%   [T_surf_C] = tima_full_model(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
%    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
%    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
%    pressure_air_pa,Observed_Temperatures,NDAYS,'depth_transition',0.5)
%
% varargin options:
%   'IncreaseSampling': Option to increase sampling rate if model does not
%       converge (default=false) (logical)
%   'Initialize': Option to reinitialize subsurface temperatures for run. Increases compute time, 
%       but avoids dramatic shifts with param optimization leading to instability. (default=false) (logical)
%   'Skin': Option to fit surface skin temperature using Newton's method. If false, T_surf is the temperature of the top layer (default=false) (logical)
%
% Input Parameters ('single quotes' means optional varargin)
%   Timeseries data
%       air_Temp_C: [C] near surface air temperature, typically at 3 m AGL.
%           (1D vector)
%       'albedo': [0-1] time variant albedo (e.g., due to wetting) for
%           fitting surface. (1D vector)
%       evap_depth: [m] Depth of the evaporation front (e.g., as function
%           of VWC). (1D vector)
%       'f_diff': [0-1] Fraction of  Global Horizontal Irradiance (GHI)
%           or r_short_upper that is diffuse. (1D vector)
%       RH: [0-1] array of near surface relative humidity values,
%           typically at 3m AGL. (1D vector)
%       pressure_air_pa: [Pa] station pressure, typically at 3m AGL. (1D vector)
%       r_long_upper: [W/m^2] Integrated longwave radiation  (4.5 to 42 μm) incident on flat
%           surface.(1D vector)
%       r_short_lower: [W/m^2] Integrated upwelling shortwave radiation (305 to 2800 nm) from flat
%           surface. (1D vector)
%       r_short_upper: [W/m^2] Integrated shortwave radiation (305 to 2800 nm) incident on flat
%           surface. (1D vector)
%       'shadow_time_ind': Index of shadow_data corresponding to
%           timeseries. (1D vector)
%       'solar_azimuth_cwfromS': [degrees] Solar azimuth in degrees
%           clockwise from South, typically -180:180 or 0:360 (1D vector)
%       'solar_zenith_apparent': [degrees] Solar zenith in degrees,
%           corrected for atmospheric refraction. (1D vector)
%       T_surf_obs_C: [C] Surface temperature values to be used for
%           initialization with any data gaps filled in. (1D vector)
%       VWC_column: [0-1, decimal fraction by volume] array of volumetric water content
%           for each model layer with each time step, typically
%           interpolated. (2D vector)
%       windspeed_horiz: [m/s] Near surface horizontal wind speed,
%           typically at 3m AGL. (1D vector)
% 
%   Constants:
%       'aspect_cwfromS': [degrees] slope aspect clockwise from South. (scalar)
%       CE: [Unitless] resistance to latent heat flux coefficient, similar
%           to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2) (scalar)
%       CH: [Unitless] resistance to sensible heat flux coefficient,
%           similar to the aerodynamic scaling factor rho_air*Cp_air/(log(z1/z0)^2/Kv^2)
%           1/CH should be 0.0028-0.0075 (or CH~100-400) for smooth to
%           roughly open soils on Davenport Scale, CH is larger (more
%           resistence) with rougher topography or larger vegetation (scalar)
%       'depth_transition': [m] Depth of major transition between upper
%           material and lower material. (scalar)
%       dt: [s] Time step (scalar)
%       emissivity: [0-1] Weighted thermal emissivity over wavelength
%           range of sensor. (scalar)
%       'e_fxn': spectral emissivity function that determines emissivity as
%           function of wavelength in um. (function)
%       initial_temps: [K] Initial condition, list of center temperatures
%           for each layer at start of simulation (1D vector)
%       'k_dry_std_lower': [W/mK] bulk dry thermal conductivty of lower layer at T_std (scalar)
%       'k_dry_std_mantle': [W/mK] bulk dry thermal conductivty of topmost
%           mantling layer at T_std (scalar)
%       k_dry_std_upper: [W/mK] bulk dry thermal conductivty of upper layer at T_std (scalar)
%       layer_size: [m] List of vertical thickness for each layer
%           from top to bottom. (1D vector)
%       m: [unitless] pore-network-connectivity-parameter, typically 0-1.3 (scalar)
%       'mantle_thickness': [m] thickness of topmost mantling layer. (scalar)
%       'material': ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at the surface (char)
%       'material_lower':  ['basalt' 'amorphous' 'granite' 'clay' 'salt' 'ice']  primary mineralogy at depth (char)
%       NDAYS: [unitless] # of days to run equilib model for: More = closer to equilib, fewer = faster (scalar)
%       'rho_dry_lower': [kg/m^3] Value for density of lower layer soil.  (scalar)
%       rho_dry_upper: [kg/m^3] Value for density of soil.  (scalar)
%       'shadow_data': [0-1] Fraction of ROI/pixel in shadow corresponding to time. (1D Vector)
%       'slope_angle': [degrees] angle of slope. (scalar)
%       theta_k: [0-1, fraction by volume] conductivity soil moisture inflection
%           point - in theory this should be similar to saturation value or porosity (scalar) 
%       theta_E: [0-1, fraction by volume] latent heat soil moisture inflection
%           point (scalar)
%       T_deep: [K] Lower boundary condition, fixed temperature (scalar)
%       T_std: [K] Standard temperature; typically 300 (scalar)
%       'T_adj1': [index, temperature K] pair used to force column
%           temperature change at a given time point due to wetting
%           (1D vector)
%       'T_adj2': [index, temperature K] pair used to force a second
%           column temperature change at a given time point due to wetting
%           (1D vector)
% 
% Output Parameters:
%   T_Surf_C = [C] Surface temperature time series (1D vector)
%
% Author
%    Ari Koeppel -- Copyright 2025
% ***************
p = inputParser;
p.KeepUnmatched=true;
% option to increase sampling rate if model does not converge
p.addParameter('IncreaseSampling',false,@islogical);
% option to reinitialize subsurface temperatures for run. Increases compute time, 
% but avoids dramatic shifts with param optimization leading to instability
p.addParameter('Initialize',false,@islogical);
p.addParameter('Skin',false,@islogical);
p.parse(varargin{:});
p=p.Results;
if p.Skin
    if p.Initialize
        Subsurface_Temperatures = tima_initialize_skin(k_dry_std_upper,rho_dry_upper,m,...
        theta_k,T_std,T_deep,T_surf_obs_C,dt,layer_size,VWC_column,RH,NDAYS,material,varargin{:});
        T_Start(2:end) = Subsurface_Temperatures(end,end,:);
        T_Start(1) = T_surf_obs_C(1)+273.15;
    end
    [T_surf_C,~,~,~,~,~,~] = tima_heat_transfer_skin(k_dry_std_upper,m,CH,CE,theta_k,theta_E,rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,...
            r_short_lower,r_long_upper,windspeed_horiz,T_deep,T_Start,layer_size,...
            VWC_column,evap_depth,RH,emissivity,...
             pressure_air_pa,varargin{:});
else
    if p.Initialize
        Subsurface_Temperatures = tima_initialize(k_dry_std_upper,rho_dry_upper,m,...
        theta_k,T_std,T_deep,T_surf_obs_C,dt,layer_size,VWC_column,RH,NDAYS,material,varargin{:});
        T_Start(:) = Subsurface_Temperatures(end,end,:);
        T_Start(1) = T_surf_obs_C(1)+273.15;
    end
    [T_surf_C,~,~,~,~,~,~] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,...
            r_short_lower,r_long_upper,windspeed_horiz,T_deep,T_Start,layer_size,...
            VWC_column,evap_depth,RH,emissivity,...
             pressure_air_pa,varargin{:});
end

if any(T_surf_C==-65535) && p.IncreaseSampling
    time_orig = dt.*(0:1:length(air_temp_C)-1);
    dt = dt/2;
    time_new = 0:dt:max(time_orig);
    air_temp_C = interp1(time_orig,air_temp_C,time_new,'linear');
    r_short_upper = interp1(time_orig,r_short_upper,time_new,'linear');
    r_short_lower = interp1(time_orig,r_short_lower,time_new,'linear');
    r_long_upper = interp1(time_orig,r_long_upper,time_new,'linear');
    windspeed_horiz = interp1(time_orig,windspeed_horiz,time_new,'linear');
    VWC_column = interp1(time_orig,VWC_column,time_new,'linear');
    RH = interp1(time_orig,RH,time_new,'linear');
    evap_depth = interp1(time_orig,evap_depth,time_new,'linear');
    pressure_air_pa = interp1(time_orig,pressure_air_pa,time_new,'linear');
    T_surf_C = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,...
        r_short_lower,r_long_upper,windspeed_horiz,T_deep,T_Start,layer_size,...
        VWC_column,evap_depth,RH,emissivity,...
        pressure_air_pa,varargin{:});
    if any(~isreal(T_surf_C)) || any(isnan(T_surf_C))
        T_surf_C(~isreal(T_surf_C)) = -65535;
        T_surf_C(isnan(T_surf_C)) = -65535;
    end
    T_surf_C = interp1(time_new,T_surf_C,time_orig,'linear');
end
end
