function [T_surf_C,T_sub_C,q_latent,k_eff_dt,q_conv,q_rad,q_G] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    pressure_air_pa,varargin)
%% TIMA_HEAT_TRANSFER
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the temperature at the emitting surface of a 1D soil column through time.

% Syntax
%   [T_surf_C,~,~,~,~,~,~,~] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
%    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
%    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
%    pressure_air_pa,'depth_transition',0.5,'k_dry_std_lower',1)
%
%   [T_surf_C,T_sub_C,q_latent,k_eff_dt,q_conv,q_rad,q_G] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
%    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
%    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
%    pressure_air_pa,'depth_transition',0.5,'k_dry_std_lower',1)
%
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
%   T_sub_C = [C] Suburface temperature time series for all layers (2D vector)
%   q_latent = [W/m^2] Latent heat flux time series for all layers (2D vector)
%   k_eff_dt = [W/mK] Effective thermal conductivity time series for all layers (2D vector)
%   q_conv = [W/m^2] Surface sensible heat flux time series (1D vector)
%   q_rad = [W/m^2] Surface radiative heat flux time series (1D vector)
%   q_G = [W/m^2] Heat flux between top and second subsurface layer time series (1D vector)
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% Sources
%   Subsurface heat flux: Kieffer et al., 2013
%   Irradiance Calculation: https://github.com/sandialabs/MATLAB_PV_LIB
%   Sensible heat flux: Cellier 1996
%   Latent heat flux: Daamen and Simmonds 1996 + Mahfouf and Noilhan
%       (1991)+ Kondo1990
%   Thermal conductivity mixing: Zhang and Wang 2017, Dong 2015
%      
% ***************
p = inputParser;
p.KeepUnmatched=true;
p.addRequired('r_short_lower',@(x)all(x>=0) && isnumeric(x));
p.addRequired('r_short_upper',@(x)all(x>=0) && isnumeric(x));
p.addRequired('r_long_upper',@(x)isnumeric(x));
p.addParameter('k_dry_std_lower',k_dry_std_upper,@isnumeric);
p.addParameter('depth_transition',sum(layer_size),@isnumeric);
p.addParameter('rho_dry_lower',rho_dry_upper,@isnumeric);
p.addParameter('mantle_thickness',0,@isnumeric);
p.addParameter('k_dry_std_mantle',k_dry_std_upper,@isnumeric);
p.addParameter('material',"basalt",@ischar);
p.addParameter('material_lower',"basalt",@ischar);
p.addParameter('albedo',[],@isnumeric);
p.addParameter('slope_angle',0,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('aspect_cwfromS',[],@(x)isnumeric(x)&&isscalar(x));
p.addParameter('solar_azimuth_cwfromS',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('solar_zenith_apparent',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('f_diff',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('e_fxn',[],@(x)isa(x,'function_handle'));
p.addParameter('shadow_data',[],@isnumeric);
p.addParameter('shadow_time_ind',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('T_adj1',[]);
p.addParameter('T_adj2',[]);
p.parse(r_short_lower, r_short_upper, r_long_upper, varargin{:});
p=p.Results;

k_dry_std_lower = p.k_dry_std_lower;
depth_transition = p.depth_transition;
rho_dry_lower = p.rho_dry_lower;
albedo = p.albedo;
slope_angle = p.slope_angle;
aspect_cwfromS = p.aspect_cwfromS;
solar_azimuth_cwfromS = p.solar_azimuth_cwfromS;
solar_zenith_apparent = p.solar_zenith_apparent;
f_diff = p.f_diff;
material = p.material;
material_lower = p.material_lower;
shadow_data = p.shadow_data;
shadow_time_ind = p.shadow_time_ind;
T_adj1 = p.T_adj1;if ~isempty(T_adj1), T_adj1(1) = T_adj1(1)/(dt/60);end
T_adj2 = p.T_adj2;if ~isempty(T_adj2), T_adj2(1) = T_adj2(1)/(dt/60);end
e_fxn = p.e_fxn;
mantle_thickness = p.mantle_thickness;
k_dry_std_mantle = p.k_dry_std_mantle;
enable_melting = false; %switch to include ability of ice to melt
% ***************

% ***************
% Inputs and constants:
air_temp_K = air_temp_C+273.15;
rho_H2O = 997; %kg/m^3, does not vary much in T range so assume constant density
omega = (1+cosd(slope_angle))/2; % visible sky fraction ISOtropic sky model
% ***************
% Initial and boundary conditions setup + array initiation:
NLAY = length(layer_size);
T = NaN(length(air_temp_C),NLAY); %Matrix of Temperatures for each layer and time point, T(1,t) is surface temp array that will be fit to obs data
T(1,:) = initial_temps;
T(:,NLAY) = T_deep;
[~,D_z] = min(abs(depth_transition-cumsum(layer_size))); %index of depth for k transition
Mantle_Ind = 0;
if mantle_thickness>0
    [~,Mantle_Ind] = min(abs(mantle_thickness-cumsum(layer_size))); %index of depth for surfce mantle transition
end
%m_melt = zeros(NLAY,1);%area density of melted ice per layer kg/m^2
k = k_dry_std_upper.*ones(1,NLAY);
kay_upper = k_dry_std_upper;
kay_lower = k_dry_std_lower;

q_G = zeros(length(air_temp_C),1);
q_rad = zeros(length(air_temp_C),1);
q_conv = zeros(length(air_temp_C),1);
q_latent = zeros(length(air_temp_C),NLAY);
k_eff_dt = zeros(length(air_temp_C),NLAY);

FC = NaN(NLAY-1,1);
FB = layer_size(2)/layer_size(1);
for z = 1:NLAY-1
    FC(z) = 2*dt/layer_size(z)^2;
end
for t = 2:length(air_temp_C)
    %*********LATENT HEAT & THERMAL CONDUCTIVITY***********
    [q_evap_1,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,1),VWC_column(t-1,1));
    q_latent(t,1) = q_evap_1;
    q_latent(1,1) = q_evap_1; %**CHANGE to q_latent(1,1?)
    if Mantle_Ind >= 1
        k(1) = tima_conductivity_model_lu2007(k_dry_std_mantle,T(t-1,1),T_std,VWC_column(t-1,1),theta_k,m,Soil_RH,material);%Top
    else
        k(1) = tima_conductivity_model_lu2007(kay_upper,T(t-1,1),T_std,VWC_column(t-1,1),theta_k,m,Soil_RH,material);%Top
    end
    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,NLAY),VWC_column(t-1,end));
    k(NLAY) = tima_conductivity_model_lu2007(kay_lower,T(t-1,NLAY),T_std,VWC_column(t-1,end),theta_k,m,Soil_RH,material_lower);%Top;%Bottom
    %****************************************
    
    %*********Specific heat***********
    rho = rho_dry_upper + rho_H2O*VWC_column(t-1,1); %Dry density + water content
    Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,1),material);
    %****************************************

    %*********SENSIBLE HEAT***********
    q_conv(t) = tima_sensible_heat_model(CH,windspeed_horiz(t-1),air_temp_K(t-1),T(t-1,1)); %(W/m^2)
    %*******************************
    %********Surface emission*********
    if isempty(e_fxn) %Treat as black body using Stefan-Boltzmann Law
        sigma = 5.670374419e-8; % Stefan-Boltzmann Const
        r_long_lower = emissivity*sigma*T(t-1,1)^4; %(W/m^2)
    else %Incorperate spectral emission variability via Planck's law
        k_b = 1.380649000000000e-23; %Boltzmann constant (J/K)
        h = 6.626070150000000e-34; %Planck constant (J/Hz)
        c = 299792458; %speed of light in vacuum (m/s)
        M = @(lambda) 2.*pi.*h.*c.^2.*e_fxn(lambda)'./(lambda.^5)./(exp(h.*c./(lambda.*k_b.*T(t-1,1)))-1);
        r_long_lower = integral(M,4.5E-6,42E-6); %Integrate radiance over CNR4 range with wavelength in meters (W/m^2)
    end
    %*********Radiative Flux***********
    if isempty(solar_azimuth_cwfromS) || isempty(solar_zenith_apparent) || isempty(aspect_cwfromS) || isempty(f_diff) && slope_angle == 0 %single ROI data over flat terrain
        if isempty(albedo)
            q_rad(t) = r_short_upper(t-1)-r_short_lower(t-1)+emissivity*r_long_upper(t-1)-r_long_lower;%net heat flux entering surface
        else
            q_rad(t) = (1-albedo(t-1))*r_short_upper(t-1)+emissivity*r_long_upper(t-1)-r_long_lower; %net heat flux entering surface
        end
    else %variably shadowed/sloped terrain
        phi = solar_azimuth_cwfromS(t-1);
        Z = solar_zenith_apparent(t-1);
        Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);
        Incidence(Incidence>1) = 1; Incidence(Incidence<0) = 0;
        if isempty(shadow_data)
            q_rad(t) = (1-albedo(t-1))*Incidence*(1-f_diff(t-1))*r_short_upper(t-1)/cosd(Z)+(1-albedo(t-1))*omega*f_diff(t-1)*r_short_upper(t-1)+(1-albedo(t-1))*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower; %net heat flux entering surface
        else
            q_rad_full = (1-albedo)*Incidence*(1-f_diff(t-1))*r_short_upper(t-1)/cosd(Z)+(1-albedo)*omega*f_diff(t-1)*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower; %net heat flux entering surface
            Lit_Fraction = shadow_data(shadow_time_ind(t-1)); %Shadow = 0, unshadowed = 1;
            q_rad(t) = Lit_Fraction*q_rad_full+(1-Lit_Fraction)*((1-albedo)*omega*f_diff(t-1)*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower); %net heat flux entering surface
        end
    end
    %*******************************
    
    %*********Ground Flux***********
    %ref: Kieffer, 2013
    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,2),VWC_column(t-1,2));
    k(2) = tima_conductivity_model_lu2007(kay_lower,T(t-1,2),T_std,VWC_column(t-1,2),theta_k,m,Soil_RH,material);
    q_G(t) = 2*(T(t-1,2)-T(t-1,1))/(layer_size(1)/k(1)+layer_size(2)/k(2)); %heat flux from conducting with lower layer
    %*******************************

    %*********Combine Heat transfer elements***********
    W = q_G(t)+q_conv(t)+q_rad(t)+q_evap_1; %Heat Flux, W/m^2
    dT = dt/(Cp*rho*layer_size(1))*(W); %Temperature change for given thickness of material with known volumetric Cp
    T(t,1) = T(t-1,1) + dT; %new T at surface
    %*******************************
    %*********Subsurface Flux (from KRC, Kieffer, 2013)***********
    for z = 2:NLAY-1 %layer loop
        if z <= Mantle_Ind && z < D_z %Surface Mantle (low k material that causes high frequency imprint on general diurnal trend
            if sum(layer_size(1:z)) <= evap_depth(t-1)
                [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
            else
                [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                q_evap_z = 0; %Evap_Coeff
            end
            k(z) = tima_conductivity_model_lu2007(k_dry_std_mantle,T(t-1,z),T_std,VWC_column(t-1,z),theta_k,m,Soil_RH,material);
            rho = rho_dry_upper + rho_H2O*VWC_column(t-1,z); %H2O dep
            Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,z),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
        elseif z < D_z % Top Layer
            if sum(layer_size(1:z)) <= evap_depth(t-1)
                [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
            else
                [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                q_evap_z = 0; %Evap_Coeff
            end
            k(z) = tima_conductivity_model_lu2007(kay_upper,T(t-1,z),T_std,VWC_column(t-1,z),theta_k,m,Soil_RH,material);
            rho = rho_dry_upper + rho_H2O*VWC_column(t-1,z); %H2O dep
            Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,z),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
        else %Lower Layer
            if strcmp(material_lower,"ice")
                rho = (917+rho_dry_upper)/2; %kg/m^3 50% ice-rich permafrost by volume
                Cp = tima_specific_heat_model_DV1963(rho,rho,T(t-1,z),material_lower);
                q_evap_z = 0; %Evap_Coeff
                %Soil_RH = 1;
                k(z:end) = 632./T(t-1,z)+0.38-0.00197.*T(t-1,z); %Wood 2020/Andersson and Inaba 2005;
            else %Not ice
                if sum(layer_size(1:z)) <= evap_depth
                    [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                else
                    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                    q_evap_z = 0; %Evap_Coeff
                end
                k(z:end) = tima_conductivity_model_lu2007(kay_lower,T(t-1,z),T_std,VWC_column(t-1,z),theta_k,m,Soil_RH,material_lower);
                rho = rho_dry_lower + rho_H2O*VWC_column(t-1,z); %H2O dep
                Cp = tima_specific_heat_model_DV1963(rho_dry_lower,rho,T(t-1,z),material_lower);
            end
        end
        % if layer_size(z)/sqrt(2*dt*(k(z)/rho/Cp)) < 1
        %     fprintf('Non convergance for k: %0.4f, m: %0.4f, CH: %0.4f, CE: %0.4f, theta_k: %0.4f, theta_E: %0.4f',k_dry_std_upper,m,CH,CE,theta_k,theta_E)
        %     T(t:end,1) = -65535;
        %     break
        % end
        %Subsurface Multip Factors (Kieffer, 2013)
        F1 = FC(z)*k(z)/((Cp*rho)*(1+FB*k(z)/k(z+1)));
        F3 = (1+FB*(k(z)/k(z+1)))/(1+(1/FB)*(k(z)/k(z-1)));
        F2 = -(1 + F3);
        %Temperature response
        dT = F1*(T(t-1,z+1)+F2*T(t-1,z)+F3*T(t-1,z-1))+dt/(Cp*rho*layer_size(z))*q_evap_z;
        T(t,z) = T(t-1,z)+dT;
        if strcmp(material_lower,"ice") && z >= D_z && T(t,z) > 273.15
            if enable_melting==1 %Melt ice if enough energy
                q_ice = (T(t,z)-273.15)*(Cp*rho*layer_size(z));%heat in to ice J/m^2 assuming 50% purity (rho and Cp still as above)
                m_melt(z) = m_melt(z)+q_ice/334000; %mass of ice melted added to prior mass. heat input (J/m^2)/ Latent heat of fusion pure ice (J/kg)=kg/m^2
                if m_melt(z) > 917*0.5*layer_size(z); D_z = D_z+1; end %if mass is greater than total mass of ice in layer at 50%: melting. 917 = Desnity of ice kg/m^3 at 0C
                %To add: elseif freezing
            end
            T(t,z) = 273.15;%set to 0 C if ice and above freezing, heat sink
        end
        q_latent(t,z) = q_evap_z;        
    end
    k_eff_dt(t,:) = k;
end
k_eff_dt(1,:) = k_eff_dt(2,:);
q_G(1)=q_G(2);
q_conv(1)=q_conv(2);
q_rad(1)=q_rad(2);
T_surf_C = T(:,1) - 273.15;
T_sub_C = T(:,2:end)-273.15;
if any(~isreal(T_surf_C)) || any(isnan(T_surf_C))
    fprintf('Non convergance for k: %0.4f, m: %0.4f, CH: %0.4f, CE: %0.4f, theta_k: %0.4f, theta_E: %0.4f \n',k_dry_std_upper,m,CH,CE,theta_k,theta_E)
    T_surf_C(~isreal(T_surf_C)) = -65535;
    T_surf_C(isnan(T_surf_C)) = -65535;
end
end