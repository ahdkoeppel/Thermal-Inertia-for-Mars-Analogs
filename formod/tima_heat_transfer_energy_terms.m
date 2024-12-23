function [T_surf_C,T_sub_C,q_evap,k_eff_dt,q_conv,q_rad,q_G]= tima_heat_transfer_energy_terms(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,RH,emissivity,...
    pressure_air_pa,varargin)

% TIMA_HEAT_TRANSFER
%   This function uses observational data and assigned thermophysical properties 
%   to estimate the emitting surface temperature through time.

% Syntax
%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,'material',material)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,'T_adj1',[index, T_K])

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,k_dry_std_lower,depth_transition)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,k_dry_std_lower,depth_transition,'rho_dry_lower',rho_dry_lower)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,'albedo',albedo)

%   [T_Surf_C] = tima_heat_transfer(k_dry_std_upper,m,CH,CE,theta_k,theta_E,...
    % rho_dry_upper,dt,T_std,air_temp_C,r_short_upper,r_short_lower,r_long_upper,...
    % windspeed_horiz,T_deep,initial_temps,layer_size,VWC_column,evap_depth,RH,emissivity,...
    % pressure_air_pa,'slope_angle',slope_angle,'aspect_cwfromS','solar_azimuth_cwfromS',solar_azimuth_cwfromS,...
    % 'solar_zenith_apparent',solar_zenith_apparent,'f_diff',f_diff,'shadow_data',shadow_data,...
    % 'shadow_time_ind',shadow_time_ind,'MappingMode',true)

% Description
%   This model uses a surrogate optimization solver with 1-7 paramters to fit IR-image-observed
%   surface temperatures using meteorological data and surface parameters derived from
%   aerial/satellite imagery. Fitting Parameters can include any or all of: top layer thermal
%   conductivity at 300K, bottome layer thermal conductivity at 300K, Depth of thermal conductivity
%   transition, pore-network-connectivity, Surf. ex. coef. (sensible), Surf. ex. coef. (latent), Soil Moist. Inflection (conductivity) (%), Soil Moist. Infl. (latent heat) (%)

% Input Parameters
%   Mapping Mode: True or False
%   Observational data
%       dt
%       Air_Temp_C
%       R_Short_Upper: must be positive
%       R_Short_Lower: Must be positive
%       R_Long_Upper
%       WindSpeed_horiz
%       Dug_VWC
%       evap_depth
%       RH
%       rho_dry
%       emissivity
%       albedo
%       pressure_air
%       f_diff
%       material
%    Boundary conditions
%       T_std
%       ST1
%       Initial_Temps
%       T_Deep
%       Cp_Deep
%       B
%       dt
%    Fitting Parameters (each may be fit or assigned)
%       kay_upper
%       kay_lower
%       Depth_transition
%       m
%       Wind_Coeff
%       Evap_Coeff
%       theta_k
%       theta_E

%   OPTIONAL
%       T_adj1: [index, T_K] pair used to adjust temperature during wetting
%       T_adj2: [index, T_K] pair used to adjust temperature during wetting

% Outputs:
%   T_Surf_C = Surface temperature time series (deg C)
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% Sources
%   Optimization tool: https://www.mathworks.com/help/gads/table-for-choosing-a-solver.html
%   Subsurface heat flux: Kieffer et al., 2013
%   Irradiance Calculation: https://github.com/sandialabs/MATLAB_PV_LIB
%   Sensible heat flux:
%   Latent heat flux: Daamen and Simmonds 1996 + Mahfouf and Noilhan (1991)+ Kondo1990
%   Thermal conductivity mixing:
%   
% See also 
%   TIMA_HEAT_TRANSFER TIMA_INITIALIZE TIMA_LATENT_HEAT_MODEL TIMA_LN_PRIOR TIMA_SENSIBLE_HEAT_MODEL TIMA_GWMCMC TIMA_COMBINE_ROWS
% ***************
p = inputParser;
p.addRequired('r_short_lower',@(x)all(x>=0) && isnumeric(x));
p.addRequired('r_short_upper',@(x)all(x>=0) && isnumeric(x));
p.addOptional('k_dry_std_lower',k_dry_std_upper,@isnumeric);
p.addOptional('depth_transition',sum(layer_size),@isnumeric);


p.addParameter('rho_dry_lower',rho_dry_upper,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('albedo',[],@isnumeric);
p.addParameter('slope_angle',0,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('aspect_cwfromS',[],@(x)isnumeric(x)&&isscalar(x));
p.addParameter('solar_azimuth_cwfromS',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('solar_zenith_apparent',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('f_diff',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('e_fxn',[],@(x)isa(x,'cfit'));
p.addParameter('shadow_data',[],@isnumeric);
p.addParameter('shadow_time_ind',[],@isnumeric);
p.addParameter('material',"basalt",@ischar);
p.addParameter('material_lower',"basalt",@ischar);
p.addParameter('MappingMode',false,@islogical);
p.addParameter('T_adj1',[]);
p.addParameter('T_adj2',[]);
p.parse(r_short_lower, r_short_upper, varargin{:});

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
MappingMode = p.MappingMode;
T_adj1 = p.T_adj1;if ~isempty(T_adj1), T_adj1(1) = T_adj1(1)/(dt/60);end
T_adj2 = p.T_adj2;if ~isempty(T_adj2), T_adj2(1) = T_adj2(1)/(dt/60);end
e_fxn = p.e_fxn;
% ***************

% ***************
% Inputs and constants:
air_temp_K = air_temp_C+273.15;
% k_H2O = 0.6096; %from Haigh 2012
rho_H2O = 997; %kg/m^3, does not vary much in T range so assume constant density
sigma = 5.670374419e-8; % Stefan-Boltzmann Const
k_b = 1.380649000000000e-23; 
h = 6.626070150000000e-34;
c = 299792458;
omega = (1+cosd(slope_angle))/2; % visible sky fraction ISOtropic sky model
% T_std = 300; %Standard temperature K
% ***************
% Initial and boundary conditions setup + array initiation:
NLAY = length(layer_size);
T = NaN(length(air_temp_C),NLAY); %Matrix of Temperatures for each layer and time point, T(1,t) is surface temp array that will be fit to obs data
T(1,:) = initial_temps;
T(:,NLAY) = T_deep; %Constant lower boundary temp
k = k_dry_std_upper.*ones(1,NLAY);
kay_upper = k_dry_std_upper;
kay_lower = k_dry_std_lower;
[~,D_z] = min(abs(depth_transition-cumsum(layer_size))); %index of depth for k transition

q_G = zeros(length(air_temp_C),1);
q_rad = zeros(length(air_temp_C),1);
q_conv = zeros(length(air_temp_C),1);
q_evap = zeros(length(air_temp_C),NLAY);
k_eff_dt = zeros(length(air_temp_C),NLAY);

FC = NaN(NLAY-1,1);
FB = layer_size(2)/layer_size(1);
for z = 1:NLAY-1
    FC(z) = 2*dt/layer_size(z)^2;
end

for t = 2:length(air_temp_C)
    %*********LATENT HEAT & THERMAL CONDUCTIVITY***********
    [q_evap_1,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,1),VWC_column(t-1,1));
    q_evap(t,1) = q_evap_1;
    q_evap(1,:) = q_evap_1.*ones(1,NLAY);
    k(1) = tima_conductivity_model_lu2007(kay_upper,T(t-1,1),T_std,VWC_column(t-1,1),theta_k,m,Soil_RH,material);%Top
    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,NLAY),VWC_column(t-1,end));
    k(NLAY) = tima_conductivity_model_lu2007(kay_lower,T(t-1,NLAY),T_std,VWC_column(t-1,end),theta_k,m,Soil_RH,material_lower);%Top;%Bottom
    %****************************************
    
    %*********Specific heat***********
    rho = rho_dry_upper + rho_H2O*VWC_column(t-1,1); %Dry density + water content
    Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,1),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
    %****************************************

    %*********SENSIBLE HEAT***********
    q_conv(t) = tima_sensible_heat_model(CH,windspeed_horiz(t-1),air_temp_K(t-1),T(t-1,1));
    %*******************************
    %********Spectral emissivity*********
    if isempty(e_fxn)
        r_long_lower = emissivity*sigma*T(t-1,1)^4;
    else
        M = @(lambda) 2.*pi.*h.*c.^2.*e_fxn(lambda)'./(lambda.^5)./(exp(h.*c./(lambda.*k_b.*T(t-1,1)))-1);
        r_long_lower = integral(M,4.5E-6,42E-6); %Integrate radiance over CNR4 range.
    end
    %*********Radiative Flux***********
    if MappingMode == false %Tower mode
        if isempty(solar_azimuth_cwfromS) || isempty(solar_zenith_apparent) || isempty(aspect_cwfromS) || isempty(f_diff) && slope_angle == 0
            if isempty(albedo)
                q_rad(t) = r_short_upper(t-1)-r_short_lower(t-1)+emissivity*r_long_upper(t-1)-r_long_lower;% %net heat flux entering surface assuming no transmission (Yes emissivity term in down)
            else
                q_rad(t) = (1-albedo(t-1))*r_short_upper(t-1)+emissivity*r_long_upper(t-1)-r_long_lower; %net heat flux entering surface
            end
        else
            phi = solar_azimuth_cwfromS(t-1);
            Z = solar_zenith_apparent(t-1);
            Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
            q_rad(t) = (1-albedo(t-1))*Incidence*(1-f_diff)*r_short_upper(t-1)/cosd(Z)+(1-albedo(t-1))*omega*f_diff*r_short_upper(t-1)+(1-albedo(t-1))*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower; %net heat flux entering surface
        end
    else %Mapping Mode
        phi = solar_azimuth_cwfromS(t-1);
        Z = solar_zenith_apparent(t-1);
        Incidence = cosd(Z)*cosd(slope_angle)+sind(Z)*sind(slope_angle)*cosd(phi-aspect_cwfromS);Incidence(Incidence>1) = 1; Incidence(Incidence<-1) = -1;
        q_rad_full = (1-albedo)*Incidence*(1-f_diff)*r_short_upper(t-1)/cosd(Z)+(1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower; %net heat flux entering surface
        Lit_Fraction = shadow_data(shadow_time_ind(t-1)); %Assumes shadow = 0, and unshadowed = 1;
        q_rad(t) = Lit_Fraction*q_rad_full+(1-Lit_Fraction)*((1-albedo)*omega*f_diff*r_short_upper(t-1)+(1-albedo)*(1-omega)*r_short_lower(t-1)+omega*emissivity*r_long_upper(t-1)-omega*r_long_lower); %net heat flux entering surface
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
    if ~isempty(T_adj1)
        if VWC_column(t,1) > 0.03 && (t == ceil(T_adj1(1)))
            T(t,1) = T_adj1(2); %Value taken from T109 probe in watering can
        end
    end
    if ~isempty(T_adj2)
        if VWC_column(t,1) > 0.02 && (t == ceil(T_adj2(1)))
            T(t,1) = T_adj2(2); %Value taken from T109 probe in watering can
        end
    end
    %*********Subsurface Flux (from KRC, Kieffer, 2013)***********
    for z = 2:NLAY-1 %layer loop
        if z < D_z
            if sum(layer_size(1:z)) <= evap_depth(t) %Wang 2016 inflection point
                [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
            else
                [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                q_evap_z = 0; %Evap_Coeff
            end
            k(z) = tima_conductivity_model_lu2007(kay_upper,T(t-1,z),T_std,VWC_column(t-1,z),theta_k,m,Soil_RH,material);
            rho = rho_dry_upper + rho_H2O*VWC_column(t-1,z); %H2O dep
            Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,z),material);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
        else
            if material_lower == "ice"
                if T(t-1,z) > 273.15, T(t-1,z) = 273.15; end
                rho = 1500; %kg/m^3
                Cp = tima_specific_heat_model_DV1963(rho,rho,T(t-1,z),material_lower);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
                q_evap_z = 0; %Evap_Coeff
                Soil_RH = 1;
                k(z:end) = 632./T(t-1,z)+0.38-0.00197.*T(t-1,z); %Wood 2020/Andersson and Inaba 2005;
            else
                if sum(layer_size(1:z)) <= evap_depth(t)
                    [q_evap_z,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                else
                    [~,Soil_RH] = tima_latent_heat_model_LP1992(CE,theta_E,pressure_air_pa(t-1),windspeed_horiz(t-1),RH(t-1),air_temp_K(t-1),T(t-1,z),VWC_column(t-1,z));
                    q_evap_z = 0; %Evap_Coeff
                end
                k(z:end) = tima_conductivity_model_lu2007(kay_lower,T(t-1,z),T_std,VWC_column(t-1,z),theta_k,m,Soil_RH,material_lower);
                rho = rho_dry_lower + rho_H2O*VWC_column(t-1,z); %H2O dep
                Cp = tima_specific_heat_model_DV1963(rho_dry_upper,rho,T(t-1,z),material_lower);%tima_specific_heat_model_hillel(rho_dry_upper,rho);%
            end
        end
        % if layer_size(z)/sqrt(2*dt*(k(z)/rho/Cp)) < 1
        %     error('Non convergance for k: %0.4f, m: %0.4f, CH: %0.4f, CE: %0.4f, theta_k: %0.4f, theta_E: %0.4f',k_dry_std_upper,m,CH,CE,theta_k,theta_E)
        % end
        F1 = FC(z)*k(z)/((Cp*rho)*(1+FB*k(z)/k(z+1)));
        F3 = (1+FB*(k(z)/k(z+1)))/(1+(1/FB)*(k(z)/k(z-1)));
        F2 = -(1 + F3);
        %Temperature response
        dT = F1*(T(t-1,z+1)+F2*T(t-1,z)+F3*T(t-1,z-1))+dt/(Cp*rho*layer_size(z))*q_evap_z;
        T(t,z) = T(t-1,z)+dT;
        if ~isempty(T_adj1)
            if VWC_column(t,z) > 0.02 && (t == ceil(T_adj1(1)))
                T(t,z) = T_adj1(2); %Value taken from T109 probe in watering can
            end
        end
        if ~isempty(T_adj2)
            if VWC_column(t,z) > 0.02 && (t == ceil(T_adj2(1)))
                T(t,z) = T_adj2(2); %Value taken from T109 probe in watering can
            end
        end
        q_evap(t,z) = q_evap_z;        
    end
    k_eff_dt(t,:) = k;
end
k_eff_dt(1,:) = k_eff_dt(2,:);
q_G(1)=q_G(2);
q_conv(1)=q_conv(2);
q_rad(1)=q_rad(2);
T_surf_C = T(:,1) - 273.15;
T_sub_C = T(:,2:end)-273.15;
